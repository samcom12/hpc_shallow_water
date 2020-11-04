#include "kernels.cuh"
#include <stdio.h>
#include <assert.h>
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) getchar();
   }
}
// Global Functions (call assignated device function)

__global__ void enforce_BC_kernel(double *d_Ht, double *d_HUt, double *d_HVt,
	int nx){
  	enforce_BC_device(d_Ht,d_HUt,d_HVt,nx);
	}

__global__ void
__launch_bounds__(256, 6)
 FV_iterator_kernel(double *d_H, double *d_HU, double *d_HV,
	const double *d_Zdx, const double *d_Zdy, double *d_Ht, double *d_HUt,
	double *d_HVt, double C, double dt, int nx){
		FV_iterator_device(d_H,d_HU,d_HV,d_Zdx,d_Zdy,d_Ht,d_HUt,d_HVt,C,dt,nx);
	}


__global__ void impose_tolerances_kernel(double *d_Ht, double *d_HUt,
	double *d_HVt, int nx){
			impose_tolerances_device(d_Ht,d_HUt,d_HVt,nx);
		}

__global__ void find_mumax_kernel(const double *H, const double *HU, const double *HV, int *mutex,
	double* mu,int numElements){
			find_mumax_device(H, HU, HV, mutex, mu, numElements);
		}

// Device Functions

__device__ void enforce_BC_device(double *d_Ht, double *d_HUt, double *d_HVt,
		int nx){
			unsigned int idx = threadIdx.x + blockIdx.x*blockDim.x;
			int offset;
			int stride = gridDim.x * blockDim.x;

			if(blockIdx.y == 0){ // Upper side
				offset = 0;
				while((idx + offset)< nx){
						d_Ht [idx]  = d_Ht [idx + nx];
						d_HUt[idx]  = d_HUt[idx + nx];
						d_HVt[idx]  = d_HVt[idx + nx];
						offset += stride;
					}
		  }
			if(blockIdx.y == 1){ // Right side
				offset = 0;
				while((idx + offset)< nx){
						d_Ht [nx-1 + idx * nx]  = d_Ht [nx-1 + idx * nx - 1];
						d_HUt[nx-1 + idx * nx]  = d_HUt[nx-1 + idx * nx - 1];
						d_HVt[nx-1 + idx * nx]  = d_HVt[nx-1 + idx * nx - 1];
						offset += stride;
					}
		  }
			if(blockIdx.y == 2){ // Downer side
				offset = 0;
				while((idx + offset)< nx){
						d_Ht [(nx-1)*(nx)+idx]  = d_Ht [(nx-1)*(nx)+idx - nx];
						d_HUt[(nx-1)*(nx)+idx]  = d_HUt[(nx-1)*(nx)+idx - nx];
						d_HVt[(nx-1)*(nx)+idx]  = d_HVt[(nx-1)*(nx)+idx - nx];
						offset += stride;
					}
		  }
			if(blockIdx.y == 3){ // Left side
				offset = 0;
				while((idx + offset)< nx){
						d_Ht [idx * nx]  = d_Ht [idx * nx + 1];
						d_HUt[idx * nx]  = d_HUt[idx * nx + 1];
						d_HVt[idx * nx]  = d_HVt[idx * nx + 1];
						offset += stride;
					}
		  }
    }

__device__ void FV_iterator_device(double *d_H, double *d_HU, double *d_HV,
		const double *d_Zdx, const double *d_Zdy, double *d_Ht, double *d_HUt,
		double *d_HVt, double C, double dt, int nx){
			unsigned int idx		= threadIdx.x + blockIdx.x*blockDim.x;
			unsigned int y;
			unsigned int x;
			unsigned int offset = 0;
      if(idx<(nx-2)*(nx-2)){
				y = (idx+offset)/(nx-2)+1;
				x = (idx+offset)%(nx-2)+1;
						d_H[y * (nx) + x] =
							0.25*( d_Ht [y * (nx) + (x+1)] + d_Ht [y 		 * (nx) + (x-1)]
										+d_Ht [(y+1) * (nx) + x] + d_Ht [(y-1) * (nx) + x])
						 +C   *( d_HUt[(y-1) * (nx) + x] - d_HUt[(y+1) * (nx) + x]
										+d_HVt[y * (nx) + (x-1)] - d_HVt[y 		 * (nx) + (x+1)]);

						d_HU[y * (nx) + x] =
							0.25*( d_HUt [y * (nx) + (x+1)] + d_HUt[y * (nx) + (x-1)]
										+d_HUt [(y+1) * (nx) + x] + d_HUt[(y-1) * (nx) + x])
						 -dt * g * d_H[y * (nx) + x]*d_Zdx[y * (nx) + x]
						 +C   *( d_HUt[y * (nx) + (x-1)]*d_HVt[y * (nx) + (x-1)]/d_Ht[y * (nx) + (x-1)]
										-d_HUt[y * (nx) + (x+1)]*d_HVt[y * (nx) + (x+1)]/d_Ht[y * (nx) + (x+1)])
						 +C   *( pow(d_HUt[(y-1) * (nx) + x],2)/d_Ht[(y-1) * (nx) + x]
										+0.5 * g * pow(d_Ht[(y-1) * (nx) + x],2)
										-pow(d_HUt[(y+1) * (nx) + x],2)/d_Ht[(y+1) * (nx) + x]
										-0.5 * g * pow(d_Ht[(y+1) * (nx) + x],2));

						d_HV[y * (nx) + x]  =
							0.25*( d_HVt[y * (nx) + (x+1)] + d_HVt[y * (nx) + (x-1)]
										+d_HVt[(y+1) * (nx) + x] + d_HVt[(y-1) * (nx) + x])
						 -dt * g * d_H[y * (nx) + x]*d_Zdy[y * (nx) + x]
						 +C   *( d_HUt[(y-1) * (nx) + x]*d_HVt[(y-1) * (nx) + x]/d_Ht[(y-1) * (nx) + x]
										-d_HUt[(y+1) * (nx) + x]*d_HVt[(y+1) * (nx) + x]/d_Ht[(y+1) * (nx) + x])
						 +C   *( pow(d_HVt[y * (nx) + (x-1)],2)/d_Ht[y * (nx) + (x-1)]
										+0.5 * g * pow(d_Ht[y * (nx) + (x-1)],2)
										-pow(d_HVt[y * (nx) + (x+1)],2)/d_Ht[y * (nx) + (x+1)]
									  -0.5 * g * pow(d_Ht[y * (nx) + (x+1)],2));
      }
		}

__device__ void impose_tolerances_device(double *d_Ht, double *d_HUt,
	double *d_HVt, int nx){
		unsigned int idx = threadIdx.x + blockIdx.x*blockDim.x;
		if(idx < nx * nx){
			if(d_Ht[idx]<0){
				d_Ht[idx] = 1e-5;
			}
			if(d_Ht[idx] <= 1e-5){
				d_HUt[idx] = 0;
				d_HVt[idx] = 0;
			}
		}
	}

__device__ void find_mumax_device(const double *H, const double *HU,
  const double *HV, int *mutex,double* mu,int numElements){

 	 	unsigned int index 						= threadIdx.x + blockIdx.x*blockDim.x;
 	 	unsigned int stride 					= gridDim.x*blockDim.x;
 	 	unsigned int offset						= 0;
		double 			 mu_ij 						= 0.0;
		double 			 mu_max 					= -1.0;
		__shared__ double cache[256];
    if(index == 0) mu[0] = -1.0;
		while((index+offset)<numElements){
	    mu_ij   = sqrt(pow(fmaxf(abs(HU[index+offset]/H[index+offset]
                                    -sqrt(H[index+offset]*g)
                                  ),
															 abs(HU[index+offset]/H[index+offset]
                                    +sqrt(H[index+offset]*g)
                                  )
                              )
                        ,2)
			              +pow(fmaxf(abs(HV[index+offset]/H[index+offset]
                                    -sqrt(H[index+offset]*g)
                                  ),
										           abs(HV[index+offset]/H[index+offset]
                                    +sqrt(H[index+offset]*g)
                                  )
                                )
                        ,2)
                      );
			mu_max 	= fmaxf(mu_max,mu_ij);
			offset += stride;
		}
	 	cache[threadIdx.x] = mu_max;
	 	__syncthreads();
	 	// reduction
	 	unsigned int i = blockDim.x/2;
	 	while(i != 0){
	 		if(threadIdx.x < i){
	 			cache[threadIdx.x] = fmaxf(cache[threadIdx.x], cache[threadIdx.x + i]);
	 		}
	 		__syncthreads();
	 		i /= 2;
	 	}
	 	if(threadIdx.x == 0){
	 		while(atomicCAS(mutex,0,1) != 0);  //lock
			mu[0] = max(mu[0],cache[0]);
	 		atomicExch(mutex, 0);  //unlock
	 		}
	}

// Additionnal Functions

void copy_host2device(double * d_Ht,double * d_HUt,double * d_HVt,double * Ht,
                      double * HUt,double * HVt,size_t memsize){
    gpuErrchk(cudaMemcpy(d_Ht,   Ht,    memsize,  cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HUt,  HUt,   memsize,  cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_HVt,  HVt,   memsize,  cudaMemcpyHostToDevice));
}

void copy_device2host(double * H,double * HU,double * HV,double * d_H,
                      double * d_HU,double * d_HV,size_t memsize){
    gpuErrchk(cudaMemcpy(H,     d_H,  memsize,  cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(HU,    d_HU, memsize,  cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(HV,    d_HV, memsize,  cudaMemcpyDeviceToHost));
}
