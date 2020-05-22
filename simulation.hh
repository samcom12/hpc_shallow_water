#ifndef SIMULATION_HH
#define SIMULATION_HH

/* -------------------------------------------------------------------------- */
#include "double_buffer.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
#include "math.h"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
class Simulation {
public:
  Simulation(int nx_, double size_, double tend_);

  /// set the initial conditions
  virtual void set_initial_conditions();

  /// set the initial conditions
  virtual void save_results();

  /// perform the simulation
  int compute();

  /// access the precision
  void set_treshold(double epsilon_);
  double treshold() const;

protected:
  /// compute one step and return an error
  virtual void compute_step();

  /// enforce zero gradient boundary condition
  virtual void boundary_conditions();

private:
  /// Global problem size
  int nx;

  /// Size of map, Size*Size [km]
  int size;

  /// Grid spacing
  double dx;

  /// Simulation time in hours [hr]
  double T = 0, Tend;

  /// Inactivity treshold
  double epsilon;

  /// Hight and velocity grids storage
  DoubleBuffer h, hu, hv;

  /// Topography storage
  Grid zdx, zdy;

};

#endif /* SIMULATION_HH */
