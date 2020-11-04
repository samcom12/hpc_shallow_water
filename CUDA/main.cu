#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

#include "kernels.cuh" //parallel kernels
#include "utils.h"     //sequential functions

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) getchar();
    }
}

int main(int argc, char *argv[])
{
    int Size = 500;      // Size of map, Size*Size [km]
    int nx = 2001;       // Number of cells in each direction on the grid
    double Tend = 0.20;  // Simulation time in hours [hr]
    int Nthreads = 256;  // Number of threads per block
    int nt = 0;
    double T = 0.0;

    if (argc == 3) {
        nx = atoi(argv[1]);
        Nthreads = atoi(argv[2]);
    }

    cout << "nx = " << nx << endl;
    cout << "Nthreads = " << Nthreads << endl;

    auto filename = to_string(Tend);
    filename = filename.substr(0, filename.find(".") + 2);
    filename = "../data/Data_nx" + to_string(nx) + "_" + to_string(Size) + "km_T" + filename;

    double dx = ((double)Size)/((double)nx);           // Grid spacening
    int numElements = nx * nx;                     // Total number of elements
    size_t memsize = numElements * sizeof(double); // Memory size of one array

    // Host data
    double GPU_dt, C = 0.0;
    double *H = (double*)malloc(memsize);
    double *HU = (double*)malloc(memsize);
    double *HV = (double*)malloc(memsize);
    double *Ht = (double*)malloc(memsize);
    double *HUt = (double*)malloc(memsize);
    double *HVt = (double*)malloc(memsize);
    double *Zdx = (double*)malloc(memsize);
    double *Zdy = (double*)malloc(memsize);
    double *dt = (double *)malloc(sizeof(double));
    double *GPU_mu = (double *)malloc(sizeof(int));
    read_data(filename, H, HU, HV, Zdx, Zdy, nx);;

    // Device data
    double *d_H, *d_HU, *d_HV, *d_Ht, *d_HUt, *d_HVt;
    double *d_Zdx, *d_Zdy;
    double *d_mu;
    int *d_mutex;

    gpuErrchk(cudaMalloc((void **) &d_H, memsize));
    gpuErrchk(cudaMalloc((void **) &d_HU, memsize));
    gpuErrchk(cudaMalloc((void **) &d_HV, memsize));
    gpuErrchk(cudaMalloc((void **) &d_Ht, memsize));
    gpuErrchk(cudaMalloc((void **) &d_HUt, memsize));
    gpuErrchk(cudaMalloc((void **) &d_HVt, memsize));
    gpuErrchk(cudaMalloc((void **) &d_Zdx, memsize));
    gpuErrchk(cudaMalloc((void **) &d_Zdy, memsize));
    gpuErrchk(cudaMalloc((void **) &d_mu, sizeof(double)));
    gpuErrchk(cudaMalloc((void **) &d_mutex, sizeof(int)));
    gpuErrchk(cudaMemset(d_mu, 0.0, sizeof(double)));
    gpuErrchk(cudaMemset(d_mutex, 0, sizeof(float)));

    // Copy data from host to device
    gpuErrchk(cudaMemcpy(d_H, H, memsize, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HU, HU, memsize, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HV, HV, memsize, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_Zdx, Zdx, memsize, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_Zdy, Zdy, memsize, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_Ht,  Ht, memsize,  cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HUt, HUt, memsize,  cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HVt, HVt, memsize,  cudaMemcpyHostToDevice));

    cpy_to(Ht,H,numElements);
    cpy_to(HUt,HU,numElements);
    cpy_to(HVt,HV,numElements);
    copy_host2device(d_H, d_HU, d_HV, H, HU, HV, memsize);
    copy_host2device(d_Ht, d_HUt, d_HVt, Ht, HUt, HVt, memsize);

    // Set grid and block dimensions
    // For finite volume kernel
    int NblocksFV = ((nx-2)*(nx-2) + Nthreads -1) / Nthreads;
    // For tolerances kernel
    int NblocksBC = pow(2, ceil(log2(nx))) / Nthreads;
    dim3 GridDimBC(NblocksBC,4);
    // For Enforce BC kernel
    int NblocksTol = (nx*nx + Nthreads -1) / Nthreads;

    const clock_t begin_time = clock();

    while (T < Tend) {
        find_mumax_kernel<<<Nthreads, Nthreads>>>(d_H, d_HU, d_HV, d_mutex, d_mu, numElements);
        cudaMemcpy(GPU_mu, d_mu, sizeof(double), cudaMemcpyDeviceToHost);
        GPU_dt = dx / (sqrt(2.0) * GPU_mu[0]);
        if(T + GPU_dt > Tend) {
            GPU_dt = Tend - T;
        }

        // Print status
        cout << "Computing T: " << T + GPU_dt << ". " << 100 * (T + GPU_dt) / Tend << "%\n";

        swap(d_H, d_HU, d_HV, d_Ht, d_HUt, d_HVt);
        cudaDeviceSynchronize();

        enforce_BC_kernel<<<GridDimBC,Nthreads>>>(d_Ht, d_HUt, d_HVt, nx);
        cudaDeviceSynchronize();

        C = (.5 * GPU_dt / dx);
        FV_iterator_kernel<<<NblocksFV,Nthreads>>>(d_H, d_HU, d_HV, d_Zdx, d_Zdy, d_Ht, d_HUt, d_HVt, C, GPU_dt, nx);
        cudaDeviceSynchronize();

        impose_tolerances_kernel<<<NblocksTol,Nthreads>>>(d_H, d_HU, d_HV, numElements);
        T = T + GPU_dt;
        nt++;
        cudaDeviceSynchronize();
    }

    cudaMemcpy(H, d_H, memsize, cudaMemcpyDeviceToHost);

    filename = to_string(Tend);
    filename = filename.substr(0, filename.find(".") + 2);
    filename = "Solution_nx" + to_string(nx) + "_" + to_string(Size) + "km_T" + filename + "_h.bin";
    write_file(filename, H, nx);

    cudaFree(d_H);
    cudaFree(d_HU);
    cudaFree(d_HV);
    cudaFree(d_Zdx);
    cudaFree(d_Zdy);
    cudaFree(d_Ht);
    cudaFree(d_HUt);
    cudaFree(d_HVt);
    cudaFree(d_mu);
    cudaFree(d_mutex);

    free(H);
    free(HU);
    free(HV);
    free(Zdx);
    free(Zdy);
    free(Ht);
    free(HUt);
    free(HVt);
    free(GPU_mu);

    double time = (double)(clock() - begin_time) / CLOCKS_PER_SEC;

    // Communicate time-to-compute
    int ops = nt * (15 + 2 + 11 + 30 + 30 + 1) * nx^2;
    double flops = (double)ops / time;
    cout << "Time to compute solution : " << time << " seconds\n";
    cout << "Average performance      : " << flops / 1.0e9 << " gflops\n";

    return 0;
}