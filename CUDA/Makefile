NVCC= nvcc               #cuda c compiler
opt= -O3
#-O2 -g -G           #optimization flag
ARCH= -arch=sm_35       #cuda compute capability
LIBS=  
execname=CUDA

objects= main.o kernels.o utils.o

#compile
$(execname): $(objects)
	$(NVCC) $(opt) -o $(execname) $(objects) $(LIBS) 

kernels.o: kernels.cu
	$(NVCC) $(opt) $(ARCH) -c kernels.cu
main.o: main.cu
	$(NVCC) $(opt) $(ARCH) -std=c++11 -c main.cu


#clean Makefile
clean:
	rm -rf hello CUDA *.o *~

