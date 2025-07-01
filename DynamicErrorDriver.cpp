#include<DynamicErrorDriver.h>
#include<cstring>

using std::vector;

extern double start_time;
extern int verbose;

//------------------------------------------------------------

DynamicErrorDriver::DynamicErrorDriver(IoData &iod_, MPI_Comm &comm_)
                  : DynamicDriver(iod_, comm_)
{

}

//------------------------------------------------------------

//! Code reused from M2C's DynamicLoadCalculator::RunForAeroS
void DynamicErrorDriver::Run()
{
  
  // objects for the embedded surface and boundary force.
  // Aero-S will update the surface with computed displacementes.
  // We will populate force and send it to Aero-S.
  TriangulatedSurface surface;
  vector<double>      pressure;

  dlo->InitializeSurface(&surface);
  pressure.assign(surface.active_nodes, 0.0);

  // All meta-level setup, such as nearest neighbor weights and surface maps,
  // will be computed here, before time stepping.
  dlo->LoadExistingSurfaces();
  dlo->LoadExistingSolutions();
  //dlo->SetupProjectionMap(surface);
  double overhead_time = walltime();

  if(verbose>1)
    print("- Time taken for file I/O: %f sec.\n", overhead_time-start_time);

  // Output recieved surface
  lagout.OutputTriangulatedMesh(surface.X0, surface.elems);

  // ---------------------------------------------------
  // Main time loop (mimics the time loop in Main.cpp
  // ---------------------------------------------------
  double error = 0.0;
  double t = 0.0, dt = 0.0, tmax = 0.0;
  int time_step = 0;
  ComputePressures(surface, pressure, t);
  ComputeError(pressure, t, error);

  lagout.OutputResults(t, dt, time_step, surface.X0, surface.X, pressure, true);

  dt   = iod.ts.timestep;
  tmax = iod.ts.maxTime;
  
  while(t<tmax) {

    time_step++;

    if(t+dt >= tmax)
      dt = tmax - t;

    //---------------------------------
    // Move forward by one time-step
    //---------------------------------
    t += dt;
    print("Step %d: t = %e, dt = %e. Computation time: %.4e s.\n", 
          time_step, t, dt, walltime()-start_time);
 
    ComputePressures(surface, pressure, t);

    double error_step = 0;
    ComputeError(pressure, t, error_step);
    error += error_step;
    if(verbose>1) print("  - Error: %e\n", error_step);
    //if(error <= error_step) error = error_step;

    lagout.OutputResults(t, dt, time_step, surface.X0, surface.X, pressure, false);
  }

  lagout.OutputResults(t, dt, time_step, surface.X0, surface.X, pressure, true);

  print("\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("\033[0;32m            NORMAL TERMINATION            \033[0m\n");
  print("\033[0;32m==========================================\033[0m\n");
  print("Avg. Norm. Root Mean Squared Error : %e.\n", error/time_step);
  print("Total File I/O Overhead Time       : %f sec.\n", overhead_time - start_time);
  print("Total Time For Integration         : %f sec.\n", walltime()    - overhead_time);
  print("Total Computation Time             : %f sec.\n", walltime()    - start_time);
  print("\n");

}

//------------------------------------------------------------

void
DynamicErrorDriver::ComputePressures(TriangulatedSurface &surface,
                                     vector<double> &pressure, 
                                     double t)
{

  assert(dlo); // cannot be null
  dlo->ComputePressures(surface, pressure, t);

}

//------------------------------------------------------------

void
DynamicErrorDriver::ComputeError(vector<double> &pressure, double t, double &error)
{
  
  assert(dlo); // cannot be null
  error = dlo->ComputeError(pressure, t);

}

//------------------------------------------------------------

