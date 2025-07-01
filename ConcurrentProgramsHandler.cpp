/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <ConcurrentProgramsHandler.h>
#include <cassert>

//---------------------------------------------------------

ConcurrentProgramsHandler::ConcurrentProgramsHandler(IoData &iod_, MPI_Comm global_comm_, MPI_Comm &comm_)
                         : iod(iod_), global_comm(global_comm_), 
                           fest_comm(global_comm_), aeros_comm(), aeros(NULL)
{
  coupled = false; //default

  // check if FEST is coupled with Aero-S
  int aeros_color    = -1;

  // Either not coupled with anything, or if coupled with Aero-S, this is the leader FESt instantiation
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    // Within the family... Common codes allow multiple (more than 2) solvers coupled together
    coupled = true;
    //The following parameters are the same as "FLUID_ID" and "MAX_CODES" in AERO-S and AERO-F
    fest_color = 0; // my color
    maxcolor = 4; 
    aeros_color = 1; //"STRUCT_ID" in AERO-S
  }

  // simultaneous operations w/ other programs 
  if(coupled)
    SetupCommunicators();

  // create inter-communicators 
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {
    aeros_comm = c[aeros_color];
    int aeros_size(-1);
    MPI_Comm_size(aeros_comm, &aeros_size);
    assert(aeros_size>0);
  }

  // time-step size suggested by other solvers, will be updated
  dt = -1.0;
  tmax = -1.0;

  // outputs the m2c communicator
  comm_ = fest_comm;
}

//---------------------------------------------------------

ConcurrentProgramsHandler::~ConcurrentProgramsHandler()
{
  if(aeros) delete aeros;
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::InitializeMessengers(TriangulatedSurface *surf_, vector<Vec3D> *F_) //for AERO-S messengers
{
  if(iod.concurrent.aeros.fsi_algo != AerosCouplingData::NONE) {

    assert(surf_); //cannot be NULL
    assert(F_); //cannot be NULL
    aeros = new AerosMessenger(iod.concurrent.aeros, fest_comm, aeros_comm, *surf_, *F_); 

    dt = aeros->GetTimeStepSize();
    tmax = aeros->GetMaxTime();

  }
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::SetupCommunicators()
{

  MPI_Comm_rank(global_comm, &global_rank);
  MPI_Comm_size(global_comm, &global_size);

  //! split the global communicator into groups. This split call
  //! gives us the processes assigned to fest. A similar call is
  //! undertaken in aero-s as well.
  MPI_Comm_split(global_comm, fest_color + 1, global_rank, &fest_comm);
  MPI_Comm_rank(fest_comm, &fest_rank);
  MPI_Comm_size(fest_comm, &fest_size);
  assert(fest_rank<fest_size); //rank must be 0 -- (size-1)

  c.resize(maxcolor);

  c[fest_color] = fest_comm;

  vector<int> leaders(maxcolor, -1);
  vector<int> newleaders(maxcolor, -1);

  if(fest_rank == 0) {
    leaders[fest_color] = global_rank;
  }

  //! collects leaders from all programs.
  MPI_Allreduce(leaders.data(), newleaders.data(), maxcolor, MPI_INTEGER, MPI_MAX, global_comm);

  for(int i=0; i<maxcolor; i++) {
    if(i != fest_color && newleaders[i] >= 0) {
      // create a communicator between fest and program i
      int tag;
      if(fest_color < i)
        tag = maxcolor * (fest_color + 1) + i + 1;
      else
        tag = maxcolor * (i + 1) + fest_color + 1;

      MPI_Intercomm_create(fest_comm, 0, global_comm, newleaders[i], tag, &c[i]);
    }
  }

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Destroy()
{
  if(aeros)
    aeros->Destroy();

  for(int i=0; i<(int)c.size(); i++)
    MPI_Comm_free(&c[i]);

}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::CommunicateBeforeTimeStepping(/*SpaceVariable3D *coordinates_,
                               DataManagers3D *dms_,
                               std::vector<GhostPoint> *ghost_nodes_inner_,
                               std::vector<GhostPoint> *ghost_nodes_outer_,
                               GlobalMeshInfo *global_mesh_, SpaceVariable3D *V,
                               SpaceVariable3D *ID, std::set<Int3> *spo_frozen_nodes*/)
{
  aeros->CommunicateBeforeTimeStepping();
  dt = aeros->GetTimeStepSize();
  tmax = aeros->GetMaxTime();
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FirstExchange(/*SpaceVariable3D *V, double dt0, double tmax0*/)
{
  aeros->FirstExchange();
  dt = aeros->GetTimeStepSize();
  tmax = aeros->GetMaxTime();
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::Exchange(/*SpaceVariable3D *V, double dt0, double tmax0*/)
{
  aeros->Exchange();
  dt = aeros->GetTimeStepSize();
  tmax = aeros->GetMaxTime();
}

//---------------------------------------------------------

void
ConcurrentProgramsHandler::FinalExchange(/*SpaceVariable3D *V*/)
{
  aeros->FinalExchange();
  dt = aeros->GetTimeStepSize();
  tmax = aeros->GetMaxTime();
}

//---------------------------------------------------------


//---------------------------------------------------------



