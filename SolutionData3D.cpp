#include<SolutionData3D.h>
#include<iterator>

//------------------------------------------------------------

SolutionData3D::SolutionData3D()
{
  data_ptr      = nullptr;
  num_stamps    = -1;
  num_data_rows = -1;
}

//------------------------------------------------------------
/*
SolutionData3D::SolutionData3D(std::map<double, std::vector<Vec3D>> &&input_map)
{
  data_ptr      = std::make_shared<std::map<double, std::vector<Vec3D>>>(
                    std::move(input_map));
  num_stamps    = data_ptr->size();
  num_data_rows = data_ptr->begin()->second.size();
}
*/

//------------------------------------------------------------

SolutionData3D::~SolutionData3D()
{
  // smart pointer is automatically deleted.
}

//------------------------------------------------------------

bool SolutionData3D::CheckSizes()
{

  assert(data_ptr); // cannot be null

  int first_size = data_ptr->begin()->second.size();
  for(const auto &snaps : *data_ptr)
    if(snaps.second.size() != first_size) return false;

  return true;

}

//------------------------------------------------------------

//! Takes ownership of the input map. It will be destroyed in the
//! calling function.
void
SolutionData3D::MoveMap(std::map<double, std::vector<Vec3D>> &&input_map)
{
  data_ptr       = std::make_shared<std::map<double, std::vector<Vec3D>>>(
                     std::move(input_map));
  if(!CheckSizes()) {
    print_error("*** Error: Solution snapshots must contain the same number of "
                "rows.\n");
    exit_mpi();
  }
  num_stamps     = data_ptr->size();
  num_data_rows  = data_ptr->begin()->second.size();
}

//------------------------------------------------------------

std::array<double,2>
SolutionData3D::GetTimeBounds()
{

  assert(data_ptr); // should not be null
  double tmin = data_ptr->begin()->first;
  double tmax = data_ptr->rbegin()->first;
  return std::array<double,2>{tmin, tmax};

}

//------------------------------------------------------------

std::array<double,2>
SolutionData3D::GetTimeBracket(double time)
{

  assert(data_ptr); // should not be null

  auto first = data_ptr->begin();
  auto last  = data_ptr->rbegin();

  //-----------------------------------
  // time out-of-bounds (left)
  //-----------------------------------
  if(time <= first->first) {
    auto other = std::next(first);
    return {first->first, other->first};
  }

  //-----------------------------------
  // time out-of-bounds (right)
  //-----------------------------------
  if(time >= last->first) {
    auto other = std::prev(last);
    return {other->first, last->first};
  }

  //-----------------------------------
  // Within bounds
  //----------------------------------
  auto range = data_ptr->equal_range(time);

  double t0 = range.first->first;
  double t1 = range.second->first;

  auto low   = range.first;
  auto upp   = range.second;

  // find correct brackets
  if(low == data_ptr->begin() && t0 != time) {
    print_error("*** Error: Cannot find a valid time bracket "
                "for time %e.\n", time);
    exit_mpi();
  }
  else if(t0 != time)
    low = std::prev(range.first);

  t0 = low->first;
  t1 = upp->first;

  return {t0, t1};

}

//------------------------------------------------------------

void
SolutionData3D::Insert(double time, std::vector<Vec3D> &&data)
{

  if(!data_ptr)
    data_ptr = std::make_shared<std::map<double, std::vector<Vec3D>>>();

  (*data_ptr)[time] = std::move(data);

  if(!CheckSizes()) {
    print_error("*** Error: Solution snapshots must contain the same number of "
                "rows.\n");
    exit_mpi();
  }

  num_stamps        = (int)data_ptr->size();
  num_data_rows     = (int)data_ptr->begin()->second.size();

}

//------------------------------------------------------------

std::vector<Vec3D>&
SolutionData3D::GetSolutionAtTime(double t)
{
  assert(data_ptr); // should not be null
  try {
    return data_ptr->at(t); // will throw an error if not found.
  } catch (std::out_of_range &e) {
    print_error("*** Error: Time stamp for t = %e not found.\n", t);
    exit_mpi();
  }
}

//------------------------------------------------------------

void
SolutionData3D::Flatten(std::vector<double> &flat)
{

  assert(data_ptr); // should not be null

  int container_size = GetSize();
  if(flat.empty())
    flat.resize(container_size, 0.0);
  
  // additional check
  assert((int)flat.size() == container_size);
   
  int index = 0;
  for(const auto& snaps : *data_ptr) {

    const double &time             = snaps.first;
    const std::vector<Vec3D> &data = snaps.second;

    // might not be needed.
    assert(index < container_size);

    // store time stamp
    flat[index] = time;
    index++;

    // store Vec3D data
    for(const auto& vec : data) {
      for(int i=0; i<3; ++i) {
        flat[index] = vec[i];
	index++;
      }
    }
    
    //fprintf(stdout, "time = %e, soln size = %d, stamps = %d\n", time, (int)data.size(), num_stamps);
    //fprintf(stdout, "index: %d, rows: %d, size: %d.\n", index, GetRows(), container_size);

  }

}

//------------------------------------------------------------

//SolutionData3D&
void
SolutionData3D::Rebuild(const std::vector<double> &other, int rows)
{

  // assign an empty map
  if(!data_ptr)
    data_ptr = std::make_shared<std::map<double, std::vector<Vec3D>>>();

  int size = other.size();
  
  for(int index=0; index<size;) {

    double time;
    std::vector<Vec3D> data(rows, Vec3D(0.0));

    time = other[index];
    index++;

    for(int i=0; i<rows; ++i) {
      for(int d=0; d<3; ++d) {
        data[i][d] = other[index];
	index++;
      }
    }

    //fprintf(stdout, "time = %e, soln size = %d, index = %d\n", time, (int)data.size(), index);

    (*data_ptr)[time] = std::move(data);

  }

  if(!CheckSizes()) {
    print_error("*** Error: Solution snapshots must contain the same number of "
                "rows.\n");
    exit_mpi();
  }

  num_stamps    = data_ptr->size();
  num_data_rows = data_ptr->begin()->second.size();

  //return *this;

}

//------------------------------------------------------------

