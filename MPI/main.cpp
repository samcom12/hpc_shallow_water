#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <omp.h>
#include <mpi.h>
#include "simulation.hh"

using namespace std;
typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

int main(int argc, char * argv[]) {
    int Size = 500;      // Size of map, Size*Size [km]
    int nx = 2001;       // Number of cells in each direction on the grid
    double Tend = 0.20;  // Simulation time in hours [hr]
    int nt = 0;

    if (argc == 2) {
        nx = atoi(argv[1]);
    }

    auto filename = to_string(Tend);
    filename = filename.substr(0, filename.find(".") + 2);
    filename = "../data/Data_nx" + to_string(nx) + "_" + to_string(Size) + "km_T" + filename;

    MPI_Init(&argc, &argv);
    int prank, psize;

    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

    Simulation simu(nx, Size, Tend, MPI_COMM_WORLD);
    simu.set_initial_conditions(filename);

    auto start = clk::now();
    nt = simu.compute();
    auto end = clk::now();

    second time = end - start;

    // Communicate time-to-compute
    if(prank == 0) {
        int ops = nt * (15 + 2 + 11 + 30 + 30 + 1) * nx^2;
        double flops = ops / time.count();
        cout << "Time to compute solution : " << time.count() << " seconds\n";
        cout << "Average performance      : " << flops / 1.0e9 << " gflops\n";
    }

    return 0;
}
