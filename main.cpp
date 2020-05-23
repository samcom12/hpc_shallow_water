#include <iostream>
#include <fstream>
#include <string>
#include "grid.hh"
#include "simulation.hh"
#include "double_buffer.hh"

using namespace std;

int main() {
    double g = 127267.2; // Gravity, 9.82*(3.6)^2*1000 in [km / hr^2]
    int Size = 500;      // Size of map, Size*Size [km]
    int nx = 2001;       // Number of cells in each direction on the grid
    double Tend = 0.20;  // Simulation time in hours [hr]
    double dx = Size/nx; // Grid spacening

    int T = 0, nt = 0;
    double max_hu, max_hv, mu = 0, dt;


    Simulation simu(nx, Size, Tend);
    simu.set_initial_conditions();
    const clock_t begin_time = clock();
    simu.compute();
    double time = (double)(clock() - begin_time) / CLOCKS_PER_SEC;

    // Communicate time-to-compute
    int ops = nt * (15 + 2 + 11 + 30 + 30 + 1) * nx^2;
    double flops = ops / time;
    cout << "Time to compute solution : " << time << " seconds\n";
    cout << "Average performance      : " << flops / 1.0e9 << " gflops\n";

    return 0;
}
