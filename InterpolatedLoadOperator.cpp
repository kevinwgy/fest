#include<InterpolatedLoadOperator.h>
#include<MathTools/rbf_interp.hpp>
#include<MathTools/weighted_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

// helper function
double ComputeDistance(double *v1, double *v2, int dim)
{
  double dist = 0;
  for(int i=0; i<dim; ++i) {
    try {
      dist += std::pow(v1[i] - v2[i], 2);
    }
    catch (std::exception &e) {
      print_error("*** Error: Parameter size does not match "
                  "the provided dimension size in "
                  "ComputeDistance.\n");
      print(e.what());
      exit_mpi();
    }
  }

  return std::sqrt(dist);
}

//------------------------------------------------------------

InterpolatedLoadOperator::InterpolatedLoadOperator(IoData &iod_, 
                                                   MPI_Comm &comm_)
                        : DynamicLoadOperator(iod_, comm_)
{
   
  int num_points = iod_meta.numPoints;
  // resize solution container
  neighbor_solutions.resize(num_points, nullptr);
  neighbor_surfaces.resize(num_points, nullptr);

}

//------------------------------------------------------------

InterpolatedLoadOperator::~InterpolatedLoadOperator() {
 
  if(!neighbor_solutions.empty())
    for(int i=0; i<(int)neighbor_solutions.size(); ++i) {
      if(neighbor_solutions[i]) delete neighbor_solutions[i];
    }

  if(!neighbor_surfaces.empty())
    for(int i=0; i<(int)neighbor_surfaces.size(); ++i) {
      if(neighbor_surfaces[i]) delete neighbor_surfaces[i];
    }

}

//------------------------------------------------------------

void InterpolatedLoadOperator::LoadExistingSurfaces()
{

  if(!npo->NeedSurface()) return;

  int num_points = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i) {
  
    // create new "empty" triangulated surface
    neighbor_surfaces[i] = new TriangulatedSurface();

    // first get the mesh file name from file handler
    string filename      = file_handler.GetMeshFileForNeighbor(i);

    // read this mesh file and update TriangulatedSurface
    file_handler.ReadMeshFile(filename, neighbor_surfaces[i]->X, 
                              neighbor_surfaces[i]->elems);

    neighbor_surfaces[i]->X0           = neighbor_surfaces[i]->X;
    neighbor_surfaces[i]->active_nodes = neighbor_surfaces[i]->X.size();
    neighbor_surfaces[i]->X0           = neighbor_surfaces[i]->X;
    neighbor_surfaces[i]->BuildConnectivities();
    neighbor_surfaces[i]->CalculateNormalsAndAreas();

  }

}

//------------------------------------------------------------

void InterpolatedLoadOperator::LoadExistingSolutions()
{
  int num_points = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i) {
  
    // create new "empty" solution containers
    neighbor_solutions[i] = new SolutionData3D();

    // first get the solution file name from file handler
    string filename = file_handler.GetSolnFileForNeighbor(i);

    // read this solution file and update the solution container
    file_handler.ReadSolutionFile(filename, *neighbor_solutions[i]);

  } 
}

//------------------------------------------------------------

void
InterpolatedLoadOperator::SetupProjectionMap(TriangulatedSurface &surface)
{

  assert(npo); //cannot be null
  int num_points   = iod_meta.numPoints;
  for(int i=0; i<num_points; ++i)
    npo->SetupProjectionMap(surface, neighbor_surfaces[i]);

}

//------------------------------------------------------------

void
InterpolatedLoadOperator::ComputeForces(TriangulatedSurface &surface, 
                                        std::vector<Vec3D> &force,
                                        std::vector<Vec3D> &force_over_area, 
                                        double t)
{
  int num_points   = iod_meta.numPoints;
  int active_nodes = surface.active_nodes;

  assert((int)surface.X.size() == active_nodes); // cracking not supported.

  // clear existing data
  for(int i=0; i<active_nodes; ++i) {
    force[i] = Vec3D(0.0);
    force_over_area[i] = Vec3D(0.0);
  }

  vector<vector<Vec3D>> neighbor_forces(num_points, 
    vector<Vec3D>(active_nodes, Vec3D(0.0)));

  for(int i=0; i<num_points; ++i) {

    // find the time interval
    double tk, tkp;
    auto time_bounds = neighbor_solutions[i]->GetTimeBounds();

    if(t<time_bounds[0]) {
      tk = tkp = time_bounds[0];
      print_warning("- Warning: Calculating forces at %e by const extrapolation "
                    "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
    }
    else if(t>time_bounds[1]) {
      tk = tkp = time_bounds[1];
      print_warning("- Warning: Calculating forces at %e by const extrapolation "
                    "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
    }
    else {

      auto bracket = neighbor_solutions[i]->GetTimeBracket(t);
      tk  = bracket[0];
      tkp = bracket[1];

    }

    // Get solutions at tk and tkp
    vector<Vec3D> &Sk  = neighbor_solutions[i]->GetSolutionAtTime(tk);
    vector<Vec3D> &Skp = neighbor_solutions[i]->GetSolutionAtTime(tkp);
    vector<Vec3D> &S   = neighbor_forces[i];

    npo->ProjectToTargetSurface(surface, neighbor_surfaces[i], Sk);
    npo->ProjectToTargetSurface(surface, neighbor_surfaces[i], Skp);
    InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                      (double*)S.data(), 3*active_nodes);

  }

  InterpolateInMetaSpace(surface, neighbor_forces, force, force_over_area);

}

//------------------------------------------------------------

void
InterpolatedLoadOperator::ComputePressures(TriangulatedSurface &surface, 
                                           std::vector<double> &pressure,
                                           double t)
{
  int num_points   = iod_meta.numPoints;
  int active_nodes = surface.active_nodes;

  assert((int)surface.X.size() == active_nodes); // cracking not supported.

  // clear existing data
  for(int i=0; i<active_nodes; ++i) {
    pressure[i] = 0.0;
  }

  vector<vector<Vec3D>> neighbor_forces(
    num_points, vector<Vec3D>(active_nodes, Vec3D(0.0)));

  for(int i=0; i<num_points; ++i) {

    // find the time interval
    double tk, tkp;
    auto time_bounds = neighbor_solutions[i]->GetTimeBounds();

    if(t<time_bounds[0]) {
      tk = tkp = time_bounds[0];
      print_warning("- Warning: Calculating forces at %e by const extrapolation "
                    "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
    }
    else if(t>time_bounds[1]) {
      tk = tkp = time_bounds[1];
      print_warning("- Warning: Calculating forces at %e by const extrapolation "
                    "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
    }
    else {

      auto bracket = neighbor_solutions[i]->GetTimeBracket(t);
      tk  = bracket[0];
      tkp = bracket[1];

    }

    // Get solutions at tk and tkp
    vector<Vec3D> &Sk  = neighbor_solutions[i]->GetSolutionAtTime(tk);
    vector<Vec3D> &Skp = neighbor_solutions[i]->GetSolutionAtTime(tkp);
    vector<Vec3D> &S   = neighbor_forces[i];

    npo->ProjectToTargetSurface(surface, neighbor_surfaces[i], Sk);
    npo->ProjectToTargetSurface(surface, neighbor_surfaces[i], Skp);
    InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                      (double*)S.data(), 3*active_nodes);

  }

  InterpolateInMetaSpace(surface, neighbor_forces, pressure);

}
//------------------------------------------------------------

void 
InterpolatedLoadOperator::InterpolateInMetaSpace(
                          TriangulatedSurface &surface, 
                          vector<vector<Vec3D>> &solutions, 
                          vector<double> &pressure) 
{

  int active_nodes = surface.active_nodes;
  int num_points   = iod_meta.numPoints;

  // Get parameters for target and neighbor points --- sorted by distance.
  vector<double> target 
    = file_handler.GetParametersForTarget();
  std::vector<vector<double>> neighbor 
    = file_handler.GetParametersForAllNeighbors();

  vector<double> &N0 = neighbor.front();
  vector<double> &N1 = neighbor.back();
  int var_dim = target.size();
  double rmin = 
    ComputeDistance(N0.data(), target.data(), var_dim);
  double rmax = 
    ComputeDistance(N1.data(), target.data(), var_dim);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
  int nodes_per_rank = active_nodes / mpi_size;
  int remainder      = active_nodes - nodes_per_rank*mpi_size;

  assert(remainder >=0 and remainder < mpi_size);

  vector<int> counts(mpi_size, -1);
  vector<int> start_index(mpi_size, -1);

  for(int i=0; i<mpi_size; ++i) {
    counts[i]      = (i < remainder) ? nodes_per_rank + 1 : nodes_per_rank;
    start_index[i] = (i < remainder) ? (nodes_per_rank + 1)*i : nodes_per_rank*i + remainder;
  }

  assert(start_index.back() + counts.back() == active_nodes);

  int my_block_size  = counts[mpi_rank];
  int my_start_index = start_index[mpi_rank];

  //fprintf(stdout, "rank[%d] block_size %d, start_index %d.\n", mpi_rank, my_block_size, my_start_index);

  //choose a radial basis function for interpolation
  void (*phi)(int, double[], double, double[]); //a function pointer
  switch (iod_meta.basis) {
    case MetaInputData::MULTIQUADRIC :
      phi = MathTools::phi1;  break;
    case MetaInputData::INVERSE_MULTIQUADRIC :
      phi = MathTools::phi2;  break;
    case MetaInputData::THIN_PLATE_SPLINE :
      phi = MathTools::phi3;  break;
    case MetaInputData::GAUSSIAN :
      phi = MathTools::phi4;  break;
    default :
      phi = MathTools::phi4;  break; //Gaussian
  }

  // compute forces
  int index = my_start_index;
  for(int iter=0; iter<my_block_size; ++iter) {

    double xd[var_dim*num_points];
    for(int i=0; i<num_points; ++i)
      for(int j=0; j<var_dim; ++j)
        xd[var_dim*i + j] = neighbor[i][j];


    double r0; //smaller than maximum separation, larger than typical separation
    r0 = rmin + rmax;

    double fd[num_points];
    //interpolate
    for(int i=0; i<num_points; ++i)
      fd[i] = solutions[i][index].norm();

    vector<double> weight(num_points, -1.0);
    vector<double> interp(1, -1.0);
    // MathTools::rbf_weight(var_dim, num_points, xd, r0, phi, fd, weight.data());
    // MathTools::rbf_interp(var_dim, num_points, xd, r0, phi, weight.data(), 1,
    //                       target.data(), interp.data());

    MathTools::weighted_interp(var_dim, num_points, xd, r0, phi, fd, 1, 
                               target.data(), interp.data());

    pressure[index] = interp[0];

    index++;

  }	  

  MPI_Allgatherv(MPI_IN_PLACE, my_block_size, MPI_DOUBLE, pressure.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);

}

//------------------------------------------------------------

void 
InterpolatedLoadOperator::InterpolateInMetaSpace(
                          TriangulatedSurface &surface, 
                          vector<vector<Vec3D>> &solutions, 
                          vector<Vec3D> &force, 
                          vector<Vec3D> &force_over_area) 
{

  int active_nodes = surface.active_nodes;
  int num_points   = iod_meta.numPoints;

  // Get parameters for target and neighbor points --- sorted by distance.
  vector<double> target 
    = file_handler.GetParametersForTarget();
  std::vector<vector<double>> neighbor 
    = file_handler.GetParametersForAllNeighbors();

  vector<double> &N0 = neighbor.front();
  vector<double> &N1 = neighbor.back();
  int var_dim = target.size();
  double rmin = 
    ComputeDistance(N0.data(), target.data(), var_dim);
  double rmax = 
    ComputeDistance(N1.data(), target.data(), var_dim);

  int mpi_rank, mpi_size;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
  int nodes_per_rank = active_nodes / mpi_size;
  int remainder      = active_nodes - nodes_per_rank*mpi_size;

  assert(remainder >=0 and remainder < mpi_size);

  vector<int> counts(mpi_size, -1);
  vector<int> start_index(mpi_size, -1);

  for(int i=0; i<mpi_size; ++i) {
    counts[i]      = (i < remainder) ? nodes_per_rank + 1 : nodes_per_rank;
    start_index[i] = (i < remainder) ? (nodes_per_rank + 1)*i : nodes_per_rank*i + remainder;
  }

  assert(start_index.back() + counts.back() == active_nodes);

  int my_block_size  = counts[mpi_rank];
  int my_start_index = start_index[mpi_rank];

  //fprintf(stdout, "rank[%d] block_size %d, start_index %d.\n", mpi_rank, my_block_size, my_start_index);

  //choose a radial basis function for interpolation
  void (*phi)(int, double[], double, double[]); //a function pointer
  switch (iod_meta.basis) {
    case MetaInputData::MULTIQUADRIC :
      phi = MathTools::phi1;  break;
    case MetaInputData::INVERSE_MULTIQUADRIC :
      phi = MathTools::phi2;  break;
    case MetaInputData::THIN_PLATE_SPLINE :
      phi = MathTools::phi3;  break;
    case MetaInputData::GAUSSIAN :
      phi = MathTools::phi4;  break;
    default :
      phi = MathTools::phi4;  break; //Gaussian
  }

  // compute forces
  int index = my_start_index;
  for(int iter=0; iter<my_block_size; ++iter) {

    // calculate nodal area and current normal.
    Vec3D patch(0.0);
    auto elems  = surface.node2elem[index];
    int  nelems = elems.size();
    for(auto it=elems.begin(); it!=elems.end(); ++it)
      patch += surface.elemArea[*it]*surface.elemNorm[*it];
    
    patch /= nelems;

    double area   = patch.norm();
    Vec3D  normal = patch/area;

    // prepare for interpolation
    double xd[var_dim*num_points];
    for(int i=0; i<num_points; ++i)
      for(int j=0; j<var_dim; ++j)
        xd[var_dim*i + j] = neighbor[i][j];


    double r0; //smaller than maximum separation, larger than typical separation
    r0 = rmin + rmax;

    double fd[num_points];

    //interpolate
    for(int i=0; i<num_points; ++i)
      fd[i] = solutions[i][index].norm();

    vector<double> weight(num_points, -1.0);
    vector<double> interp(1, -1.0);
    // MathTools::rbf_weight(var_dim, num_points, xd, r0, phi, fd, weight.data());
    // MathTools::rbf_interp(var_dim, num_points, xd, r0, phi, weight.data(), 1,
    //                       target.data(), interp.data());

    MathTools::weighted_interp(var_dim, num_points, xd, r0, phi, fd, 1, 
                               target.data(), interp.data());

    force[index]           = interp[0]*area*normal;   
    force_over_area[index] = interp[0]*normal;

    index++;

  }	  

  // communication
  for(int i=0; i<mpi_size; i++) {
    counts[i] *= 3;
    start_index[i] *= 3;
  }
  MPI_Allgatherv(MPI_IN_PLACE, 3*my_block_size, MPI_DOUBLE, (double*)force.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);
  MPI_Allgatherv(MPI_IN_PLACE, 3*my_block_size, MPI_DOUBLE, (double*)force_over_area.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);


}

//------------------------------------------------------------
