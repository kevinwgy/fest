#include<ClosestLoadOperator.h>
#include<MathTools/rbf_interp.hpp>

using std::vector;

extern int verbose;

//------------------------------------------------------------

ClosestLoadOperator::ClosestLoadOperator(IoData &iod_, MPI_Comm &comm_)
                   : DynamicLoadOperator(iod_, comm_)
{
  //
}

//------------------------------------------------------------

ClosestLoadOperator::~ClosestLoadOperator()
{

  if(closest_solution) delete closest_solution;
  if(closest_surface) delete closest_surface;

}

//------------------------------------------------------------

void ClosestLoadOperator::LoadExistingSurfaces()
{

  if(!npo->NeedSurface()) return;
  
  closest_surface = new TriangulatedSurface();
  string filename = file_handler.GetMeshFileForNeighbor(0);
  file_handler.ReadMeshFile(filename, closest_surface->X,
                            closest_surface->elems);
  
  closest_surface->X0           = closest_surface->X;
  closest_surface->active_nodes = closest_surface->X.size();
  closest_surface->X0           = closest_surface->X;
  closest_surface->BuildConnectivities();
  closest_surface->CalculateNormalsAndAreas();  

}

//------------------------------------------------------------

void ClosestLoadOperator::LoadExistingSolutions()
{

  closest_solution = new SolutionData3D();
  string filename  = file_handler.GetSolnFileForNeighbor(0);
  file_handler.ReadSolutionFile(filename, *closest_solution);

}

//------------------------------------------------------------

void
ClosestLoadOperator::SetupProjectionMap(TriangulatedSurface &surface)
{

  assert(npo); //cannot be null
  npo->SetupProjectionMap(surface, closest_surface);

}

//------------------------------------------------------------

void
ClosestLoadOperator::ComputeForces(TriangulatedSurface& surface, 
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

  // find the time interval
  double tk, tkp;
  auto time_bounds = closest_solution->GetTimeBounds();

  if(t<time_bounds[0]) {
    auto bracket = closest_solution->GetTimeBracket(time_bounds[0]);
    tk  = bracket[0];
    tkp = bracket[1];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else if(t>time_bounds[1]) {
    auto bracket = closest_solution->GetTimeBracket(time_bounds[1]);
    tk  = bracket[0];
    tkp = bracket[1];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else {
    auto bracket = closest_solution->GetTimeBracket(t);
    tk  = bracket[0];
    tkp = bracket[1];
  }

  //fprintf(stdout, "Bracket (%e, %e) for time %e.\n", tk, tkp, t);

  // Get solutions at tk and tkp
  vector<Vec3D> &Sk  = closest_solution->GetSolutionAtTime(tk);
  vector<Vec3D> &Skp = closest_solution->GetSolutionAtTime(tkp);
  vector<Vec3D> S(active_nodes, Vec3D(0.0));

  npo->ProjectToTargetSurface(surface, closest_surface, Sk);
  npo->ProjectToTargetSurface(surface, closest_surface, Skp);
  InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                    (double*)S.data(), 3*active_nodes);

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

  int index = my_start_index;
  for(int i=0; i<my_block_size; ++i) {

    // get current area and normal
    Vec3D patch(0.0);
    auto elems  = surface.node2elem[index];
    int  nelems = elems.size();
    for(auto it=elems.begin(); it!=elems.end(); ++it)
      patch += surface.elemArea[*it]*surface.elemNorm[*it];
    
    patch /= nelems;

    double area   = patch.norm();
    Vec3D  normal = patch/area;

    // calculate forces from pressures.
    double pressure = S[index].norm();

    force[index]           = pressure*area*normal;
    force_over_area[index] = pressure*normal; 

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

void
ClosestLoadOperator::ComputePressures(TriangulatedSurface& surface, 
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

  // find the time interval
  double tk, tkp;
  auto time_bounds = closest_solution->GetTimeBounds();

  if(t<time_bounds[0]) {
    auto bracket = true_solution->GetTimeBracket(time_bounds[0]);
    tk = bracket[0];
    tkp = bracket[1];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else if(t>time_bounds[1]) {
    auto bracket = true_solution->GetTimeBracket(time_bounds[1]);
    tk = bracket[0];
    tkp = bracket[1];
    tk = tkp = time_bounds[1];
    print_warning("- Warning: Calculating forces at %e by const extrapolation "
                  "(outside data interval [%e,%e]).\n", t, time_bounds[0], time_bounds[1]);
  }
  else {
    auto bracket = closest_solution->GetTimeBracket(t);
    tk  = bracket[0];
    tkp = bracket[1];
  }

  //fprintf(stdout, "Bracket (%e, %e) for time %e.\n", tk, tkp, t);

  // Get solutions at tk and tkp
  vector<Vec3D> &Sk  = closest_solution->GetSolutionAtTime(tk);
  vector<Vec3D> &Skp = closest_solution->GetSolutionAtTime(tkp);
  vector<Vec3D> S(active_nodes, Vec3D(0.0));

  npo->ProjectToTargetSurface(surface, closest_surface, Sk);
  npo->ProjectToTargetSurface(surface, closest_surface, Skp);
  InterpolateInTime(tk, (double*)Sk.data(), tkp, (double*)Skp.data(), t, 
                    (double*)S.data(), 3*active_nodes);

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

  int index = my_start_index;
  for(int i=0; i<my_block_size; ++i) {

    pressure[index] = S[index].norm(); 
    index++;

  }

  // communication
  MPI_Allgatherv(MPI_IN_PLACE, my_block_size, MPI_DOUBLE, pressure.data(), 
                 counts.data(), start_index.data(), MPI_DOUBLE, comm);

}

//------------------------------------------------------------

