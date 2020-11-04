#ifndef g
#define g 127267.2 //Gravity, 9.82*(3.6)^2*1000 in [km / hr^2]
#endif
#ifndef __KERNELS_CUH__
#define __KERNELS_CUH__

// Global Functions (call assignated device function)

__global__ void enforce_BC_kernel(double *d_Ht, double *d_HUt, double *d_HVt,
 	int nx);

__global__ void FV_iterator_kernel(double *d_H, double *d_HU, double *d_HV,
	const double *d_Zdx, const double *d_Zdy, double *d_Ht, double *d_HUt,
	double *d_HVt, double C, double dt, int nx);

__global__ void impose_tolerances_kernel(double *d_Ht, double *d_HUt,
		double *d_HVt, int nx);

__global__ void find_mumax_kernel(const double *H, const double *HU, const double *HV, int *mutex,
  double* mu,int numElements);

// Device Functions

__device__ void enforce_BC_device(double *d_Ht, double *d_HUt, double *d_HVt,
 	int nx);

__device__ void FV_iterator_device(double *d_H, double *d_HU, double *d_HV,
	const double *d_Zdx, const double *d_Zdy, double *d_Ht, double *d_HUt,
	double *d_HVt, double C, double dt, int nx);

__device__ void impose_tolerances_device(double *d_Ht, double *d_HUt,
		double *d_HVt, int nx);

__device__ void find_mumax_device(const double *H, const double *HU, const double *HV, int *mutex,
  double* mu,int numElements);

// Additionnal Functions

void copy_host2device(double * d_Ht,double * d_HUt,double * d_HVt,double * Ht,
                      double * HUt,double * HVt,size_t memsize);

void copy_device2host(double * H,double * HU,double * HV,double * d_H,
                      double * d_HU,double * d_HV,size_t memsize);

#endif
