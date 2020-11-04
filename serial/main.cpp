#include <iostream>
#include <fstream>
#include <string>
#include "grid.hh"
#include "simulation.hh"
#include "double_buffer.hh"

using namespace std;

int main(int argc, char * argv[]) {
    int Size = 500;      // Size of map, Size*Size [km]
    int nx = 2001;       // Number of cells in each direction on the grid
    double Tend = 0.20;  // Simulation time in hours [hr]

    int nt = 0;
    
    if (argc == 2) {
        nx = atoi(argv[1]);
    }

    Simulation simu(nx, Size, Tend);
    simu.set_initial_conditions();
    const clock_t begin_time = clock();
    nt = simu.compute();
    double time = (double)(clock() - begin_time) / CLOCKS_PER_SEC;

    // Communicate time-to-compute
    int ops = nt * (15 + 2 + 11 + 30 + 30 + 1) * nx^2;
    double flops = ops / time;
    cout << "Time to compute solution : " << time << " seconds\n";
    cout << "Average performance      : " << flops / 1.0e9 << " gflops\n";

    return 0;
}
