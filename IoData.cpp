/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#include <Utils.h>
#include <IoData.h>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cmath>
#include <unistd.h>
#include <bits/stdc++.h> //INT_MAX
//#include <dlfcn.h>
using namespace std;

RootClassAssigner *nullAssigner = new RootClassAssigner;

//------------------------------------------------------------------------------

// TODO: Refactor these data structures to lagrangian probes.
ProbeNode::ProbeNode() {
  locationX = locationY = locationZ = -1.0e20;
}

//------------------------------------------------------------------------------

Assigner* ProbeNode::getAssigner()
{
  ClassAssigner *ca = new ClassAssigner("normal", 3, nullAssigner);

  new ClassDouble<ProbeNode>(ca, "X",this,&ProbeNode::locationX);
  new ClassDouble<ProbeNode>(ca, "Y",this,&ProbeNode::locationY);
  new ClassDouble<ProbeNode>(ca, "Z",this,&ProbeNode::locationZ);

  return ca;
}

//------------------------------------------------------------------------------

Probes::Probes() {

  frequency = -100;
  frequency_dt = -1;

  density = "";
  pressure = "";
  temperature = "";
  delta_temperature = "";
  velocity_x = "";
  velocity_y = "";
  velocity_z = "";
  sound_speed = "";
  materialid = "";
  laser_radiance = "";
  levelset0 = "";
  levelset1 = "";
  levelset2 = "";
  levelset3 = "";
  levelset4 = "";
  ionization_result = "";
  reference_map = "";
  principal_elastic_stresses = "";

}

//------------------------------------------------------------------------------

void Probes::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 21, father);

  new ClassInt<Probes>(ca, "Frequency", this, &Probes::frequency);
  new ClassDouble<Probes>(ca, "TimeInterval", this, &Probes::frequency_dt);
  new ClassStr<Probes>(ca, "Density", this, &Probes::density);
  new ClassStr<Probes>(ca, "Pressure", this, &Probes::pressure);
  new ClassStr<Probes>(ca, "Temperature", this, &Probes::temperature);
  new ClassStr<Probes>(ca, "DeltaTemperature", this, &Probes::delta_temperature);
  new ClassStr<Probes>(ca, "VelocityX", this, &Probes::velocity_x);
  new ClassStr<Probes>(ca, "VelocityY", this, &Probes::velocity_y);
  new ClassStr<Probes>(ca, "VelocityZ", this, &Probes::velocity_z);
  new ClassStr<Probes>(ca, "SoundSpeed", this, &Probes::sound_speed);
  new ClassStr<Probes>(ca, "MaterialID", this, &Probes::materialid);
  new ClassStr<Probes>(ca, "LaserRadiance", this, &Probes::laser_radiance);
  //KW: laser_radiance is actually irradiance. (A misnomer in the past)
  new ClassStr<Probes>(ca, "LaserIrradiance", this, &Probes::laser_radiance);
  new ClassStr<Probes>(ca, "LevelSet0", this, &Probes::levelset0);
  new ClassStr<Probes>(ca, "LevelSet1", this, &Probes::levelset1);
  new ClassStr<Probes>(ca, "LevelSet2", this, &Probes::levelset2);
  new ClassStr<Probes>(ca, "LevelSet3", this, &Probes::levelset3);
  new ClassStr<Probes>(ca, "LevelSet4", this, &Probes::levelset4);
  new ClassStr<Probes>(ca, "IonizationResult", this, &Probes::ionization_result);
  new ClassStr<Probes>(ca, "ReferenceMap", this, &Probes::reference_map);
  new ClassStr<Probes>(ca, "PrincipalElasticStresses", this, &Probes::principal_elastic_stresses);

  myNodes.setup("Node", ca);

}

//------------------------------------------------------------------------------

AerosCouplingData::AerosCouplingData()
{
  fsi_algo = NONE;
}

//------------------------------------------------------------------------------

void AerosCouplingData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 1, father);

  new ClassToken<AerosCouplingData> (ca, "FSIAlgorithm", this,
     reinterpret_cast<int AerosCouplingData::*>(&AerosCouplingData::fsi_algo), 4,
     "None", 0, "ByAeroS", 1, "C0", 2, "A6", 3);
}

//------------------------------------------------------------------------------


ConcurrentProgramsData::ConcurrentProgramsData()
{

}

//------------------------------------------------------------------------------

void ConcurrentProgramsData::setup(const char *name, ClassAssigner *father)
{
  //ClassAssigner *ca = new ClassAssigner(name, 3, father);
  new ClassAssigner(name, 1, father);
  aeros.setup("AeroS");
} 

//------------------------------------------------------------------------------

LagrangianMeshOutputData::LagrangianMeshOutputData()
{
  frequency = 0;
  frequency_dt = -1.0;

  prefix = "";

  orig_config = "";
  disp = "";
  force  = "";
  force_over_area = "";
  surface_pressure = "";

  wetting_output_filename = "";
}

//------------------------------------------------------------------------------

void LagrangianMeshOutputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 9, father);
  
  new ClassInt<LagrangianMeshOutputData>(ca, "Frequency", this, &LagrangianMeshOutputData::frequency);
  new ClassDouble<LagrangianMeshOutputData>(ca, "TimeInterval", this, &LagrangianMeshOutputData::frequency_dt);

  new ClassStr<LagrangianMeshOutputData>(ca, "Prefix", this, &LagrangianMeshOutputData::prefix);
  new ClassStr<LagrangianMeshOutputData>(ca, "Mesh", this, &LagrangianMeshOutputData::orig_config);

  new ClassStr<LagrangianMeshOutputData>(ca, "Displacement", this, &LagrangianMeshOutputData::disp);
  new ClassStr<LagrangianMeshOutputData>(ca, "Force", this, &LagrangianMeshOutputData::force);
  new ClassStr<LagrangianMeshOutputData>(ca, "ForceOverArea", this, &LagrangianMeshOutputData::force_over_area);
  new ClassStr<LagrangianMeshOutputData>(ca, "SurfacePressure", this, &LagrangianMeshOutputData::surface_pressure);

  new ClassStr<LagrangianMeshOutputData>(ca, "ContactSurfaceOutput", this,
                                         &LagrangianMeshOutputData::wetting_output_filename);

}

//------------------------------------------------------------------------------

MetaInputData::MetaInputData()
{
  metafile = "";
  //snapshot_file_prefix = "";
  //snapshot_file_suffix = "";

  basis = GAUSSIAN;
  numPoints = 3;
}

//------------------------------------------------------------------------------

void MetaInputData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 3, father);
  
  new ClassStr<MetaInputData>(ca, "MetaFile", this, &MetaInputData::metafile);

  //new ClassStr<MetaInputData>(ca, "SnapshotFilePrefix", this, 
  //        &MetaInputData::snapshot_file_prefix);

  //new ClassStr<MetaInputData>(ca, "SnapshotFileSuffix", this, 
  //        &MetaInputData::snapshot_file_suffix);

  new ClassToken<MetaInputData> (ca, "MetaInterpolationBasis", this,
     reinterpret_cast<int MetaInputData::*>(&MetaInputData::basis), 4,
     "Multiquadric", 0, "InverseMultiquadric", 1, "ThinPlateSpline", 2, "Gaussian", 3);

  new ClassInt<MetaInputData>(ca, "NumberOfBasisPoints", this, &MetaInputData::numPoints);
}

//------------------------------------------------------------------------------

SpatialInterpolationData::SpatialInterpolationData()
{
  type = EXACT;
}

//------------------------------------------------------------------------------

void SpatialInterpolationData::setup(const char *name, ClassAssigner *father)
{

  //ClassAssigner *ca = new ClassAssigner(name, 1, father);
  //
  //new ClassToken<SpatialInterpolationData> (ca, "ProjectionType", this,
  //   reinterpret_cast<int SpatialInterpolationData::*>(&SpatialInterpolationData::projection), 2,
  //   "Matched", 0, "Interpolate", 1);

}

//------------------------------------------------------------------------------

DynamicLoadCalculatorData::DynamicLoadCalculatorData()
{
  mode     = LOAD;
  type     = NONE;
  verbose  = LOW;
  pressure = 1e5;
}

//------------------------------------------------------------------------------

void DynamicLoadCalculatorData::setup(const char *name, ClassAssigner *father)
{
  ClassAssigner *ca = new ClassAssigner(name, 6, father);

  new ClassToken<DynamicLoadCalculatorData>(ca, "Mode", this,
      reinterpret_cast<int DynamicLoadCalculatorData::*>(&DynamicLoadCalculatorData::mode), 2,
      "Load", 0, "Error", 1);

  new ClassToken<DynamicLoadCalculatorData>(ca, "Type", this,
     reinterpret_cast<int DynamicLoadCalculatorData::*>(&DynamicLoadCalculatorData::type), 4,
     "None", 0, "Constant", 1, "Closest", 2, "Interp", 3);

  new ClassToken<DynamicLoadCalculatorData>(ca, "VerboseScreenOutput", this,
      reinterpret_cast<int DynamicLoadCalculatorData::*>(&DynamicLoadCalculatorData::verbose), 3,
      "Low", 0, "Medium", 1, "High", 2);

  new ClassDouble<DynamicLoadCalculatorData>(ca, "Pressure", this, &DynamicLoadCalculatorData::pressure);

  meta_input.setup("MetaInputData");
  spatial_interp.setup("SpatialInterpolationData");
} 

//------------------------------------------------------------------------------

TsData::TsData()
{

  maxIts = INT_MAX;
  timestep = -1.0;
  maxTime = 1e10;

}

//------------------------------------------------------------------------------

void TsData::setup(const char *name, ClassAssigner *father)
{

  ClassAssigner *ca = new ClassAssigner(name, 3, father);

  new ClassInt<TsData>(ca, "MaxIts", this, &TsData::maxIts);
  new ClassDouble<TsData>(ca, "TimeStep", this, &TsData::timestep);
  new ClassDouble<TsData>(ca, "MaxTime", this, &TsData::maxTime);

}

//------------------------------------------------------------------------------

IoData::IoData(int argc, char** argv)
{
  //Should NOT call functions in Utils (e.g., print(), exit_mpi()) because the
  //M2C communicator may have not been properly set up.
  readCmdLine(argc, argv);
  readCmdFile();
}

//------------------------------------------------------------------------------

void IoData::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stdout,"\033[0;31m*** Error: Input file not provided!\n\033[0m");
    exit(-1);
  }
  cmdFileName = argv[1];
}

//------------------------------------------------------------------------------

void IoData::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stdout,"\033[0;31m*** Error: could not open \'%s\'\n\033[0m", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stdout,"\033[0;31m*** Error: command file contained parsing errors.\n\033[0m");
    exit(error);
  }
  fclose(cmdFilePtr);
}

//------------------------------------------------------------------------------
// This function is supposed to be called after creating M2C communicator. So, 
// functions in Utils can be used.
void IoData::finalize()
{
  //
}

//------------------------------------------------------------------------------

void IoData::setupCmdFileVariables()
{

  concurrent.setup("ConcurrentPrograms");
  calculator.setup("DynamicsCalculator");
  output.setup("Output");
  ts.setup("Time");

}

//------------------------------------------------------------------------------
