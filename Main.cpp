#include<ctime>
#include<memory>
#include<ConcurrentProgramsHandler.h>
#include<DynamicDriver.h>
#include<DynamicLoadDriver.h>
#include<DynamicErrorDriver.h>

int verbose;
double start_time;
MPI_Comm fest_comm;

int main(int argc, char* argv[])
{

  //! Initialize MPI 
  MPI_Init(NULL,NULL); //called together with all concurrent programs -> MPI_COMM_WORLD
  start_time = walltime(); //for timing purpose only (calls MPI_Wtime)

  //! Print header (global proc #0, assumed to be a FEST proc)
  fest_comm = MPI_COMM_WORLD;
  printHeader(argc, argv);

  //! Read user's input file (read the parameters)
  IoData iod(argc, argv);
  verbose = iod.calculator.verbose;

  //! Partition MPI, if there are concurrent programs
  MPI_Comm comm; //this is going to be the FEST communicator
  ConcurrentProgramsHandler concurrent(iod, MPI_COMM_WORLD, comm);
  fest_comm = comm; //correct it

  //! Finalize IoData (read additional files and check for errors)
  iod.finalize();

  //! Special tool for nearest neighbor pressure interpolations.
  std::shared_ptr<DynamicDriver> driver;
  if(iod.calculator.mode == DynamicLoadCalculatorData::LOAD)
    driver = std::make_shared<DynamicLoadDriver>(iod, comm, concurrent);
  else
    driver = std::make_shared<DynamicErrorDriver>(iod, comm);

  //! Run time integration. 
  //! Note that FEST does not perform time integration, strictly
  //! speaking. Instead, it will compute fluid-structure interface forces 
  //! based on existing FSI simulation provided by the user.
  driver->Run();
  //MPI_Barrier(fest_comm); 

  //! finalize 
  driver->Destroy();
  concurrent.Destroy();
  MPI_Finalize();
  
  return EXIT_SUCCESS;

}
