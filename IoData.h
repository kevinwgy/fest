/************************************************************************
 * Copyright © 2020 The Multiphysics Modeling and Computation (M2C) Lab
 * <kevin.wgy@gmail.com> <kevinw3@vt.edu>
 ************************************************************************/

#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <cstdio>
#include <map>
#include <parser/Assigner.h>
#include <parser/Dictionary.h>
#include <Utils.h>

using std::map;

/*********************************************************************
 * class IoData reads and processes the input data provided by the user
 *********************************************************************
*/
//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:
  map<int, DataType *> dataMap;

  void setup(const char *name, ClassAssigner *p) {
    SysMapObj<DataType> *smo = new SysMapObj<DataType>(name, &dataMap);
    if (p) p->addSmb(name, smo);
    else addSysSymbol(name, smo);
  }

  ~ObjectMap() {
    for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      delete it->second;
  }
};

//------------------------------------------------------------------------------

// TODO: (AN) update these structures to handle lagrangian probes.
struct ProbeNode {

  double locationX;
  double locationY;
  double locationZ;

  ProbeNode();
  ~ProbeNode() {}

  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct Probes {

  ObjectMap<ProbeNode> myNodes;

  int frequency;
  double frequency_dt;

  enum Vars  {DENSITY = 0, VELOCITY_X = 1, VELOCITY_Y = 2, VELOCITY_Z = 3, PRESSURE = 4, TEMPERATURE = 5, 
              DELTA_TEMPERATURE = 6, MATERIALID = 7, LASERRADIANCE = 8, LEVELSET0 = 9, LEVELSET1 = 10, 
              LEVELSET2 = 11, LEVELSET3 = 12, LEVELSET4 = 13, IONIZATION = 14, REFERENCE_MAP = 15,
              PRINCIPAL_ELASTIC_STRESSES = 16, SOUND_SPEED = 17, SIZE = 18};

  const char *density;
  const char *velocity_x;
  const char *velocity_y;
  const char *velocity_z;
  const char *sound_speed;
  const char *pressure;
  const char *temperature;
  const char *delta_temperature;
  const char *materialid;
  const char *laser_radiance; //!< note: this is actually irradiance
  const char *levelset0;
  const char *levelset1;
  const char *levelset2;
  const char *levelset3;
  const char *levelset4;
  const char *ionization_result;
  const char *reference_map;
  const char *principal_elastic_stresses;

  Probes();
  ~Probes() {}

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct LagrangianMeshOutputData {

  int frequency;
  double frequency_dt; //!< -1 by default. To activate it, set it to a positive number

  const char* prefix; //!< path

  const char* orig_config; //!< original mesh
  const char* disp; //!< displacement
  const char* force; //!< force
  const char* force_over_area; //!< force over area
  const char* surface_pressure; //!< norm of interpolated values

  const char *wetting_output_filename; //!< optional output file that shows the detected wetted side(s)

  LagrangianMeshOutputData();
  ~LagrangianMeshOutputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------


struct AerosCouplingData {

  enum FSICouplingAlgorithm {NONE = 0, BY_AEROS = 1, C0 = 2, A6 = 3} fsi_algo;

  AerosCouplingData();
  ~AerosCouplingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ConcurrentProgramsData {

  AerosCouplingData aeros;

  ConcurrentProgramsData();
  ~ConcurrentProgramsData() {} 

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

//! Optimization level data is encapsulated here. For example, paths to existing
//! FSI simulations, number of simulations used for approximation, and radial 
//! basis function used for nearest neighbor interpolations.
struct MetaInputData {

  const char* metafile;
  //const char* snapshot_file_prefix; 
  //const char* snapshot_file_suffix; 
  
  enum BasisFunction {MULTIQUADRIC = 0, INVERSE_MULTIQUADRIC = 1, 
                      THIN_PLATE_SPLINE = 2, GAUSSIAN = 3, SIZE = 4} basis; //basis function for interpolation
  int numPoints; //number of points for (unstructured) interpolation

  MetaInputData();
  ~MetaInputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

//! This data structure contains user inputs for spacial reconstruction of
//! fluid-structure interface pressure based on exisiting FSI simulations.
//! Barycentric/radial basis interpolations are used for the reconstruction.
//! Currently, not implemented.
struct SpatialInterpolationData {

  enum Type {EXACT=0, INITIAL, UPDATED} type;
  
  SpatialInterpolationData();
  ~SpatialInterpolationData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DynamicLoadCalculatorData {

  enum Mode {LOAD = 0, ERROR = 1} mode;
  enum Type {NONE = 0, CONSTANT, CLOSEST, INTERP} type;

  enum VerbosityLevel {LOW = 0, MEDIUM = 1, HIGH = 2} verbose;

  double pressure;
  MetaInputData meta_input;
  SpatialInterpolationData spatial_interp;

  DynamicLoadCalculatorData();
  ~DynamicLoadCalculatorData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData {

  int maxIts;
  double timestep;
  double maxTime;

  TsData();
  ~TsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  ConcurrentProgramsData concurrent;

  DynamicLoadCalculatorData calculator;

  LagrangianMeshOutputData output;

  TsData ts;

public:

  IoData() {}
  IoData(int, char**);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void finalize();

};
#endif
