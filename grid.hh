#ifndef GRID_HH
#define GRID_HH

/* -------------------------------------------------------------------------- */
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <mpi.h>
/* -------------------------------------------------------------------------- */
using std::string;
using std::ifstream;
using std::ios;

class Grid {
public:
  Grid(int m, int n, MPI_Comm communicator_);

  /// read values from binary file
  void read_file(string filename);

  /// write values to binary file
  void write_file(string filename);

  /// access the value [i][j] of the grid
  inline double & operator()(int i, int j) { return m_storage[i * nx + j]; }
  inline const double & operator()(int i, int j) const {
    return m_storage[i * nx + j];
  }

  void resize(int m, int n);

  /// set the grid to 0
  void clear();

  int get_nx() const;
  int get_ny() const;

private:
  int ny, nx;
  std::vector<double> m_storage;
  MPI_Comm communicator;
};

#endif /* GRID_HH */
