# Shallow water parallel solver 
Shallow Water Wave Equation with a finite volume solver

Requirements
====

- for serial version: `icc`
- for MPI version: `mpiicpc`
- for GPU version: CUDA® compiler

Data
====

The data must be downloaded from [here](https://drive.switch.ch/index.php/s/7QFOGwphJun6mln). Data folder should be in a project folder.

Build
====
Compile serial version on `fidis`:

```
cd serial
module load intel
make
```

Compile MPI version on `fidis`:
	
```
cd MPI
module load intel intel-mpi
make 
```
	
Compile CUDA version on `deneb{1,2}`:

```
cd CUDA
slmodules -s x86_E5v2_Mellanox_GPU
module load gcc cuda       
make
```

Run
====
Serial version:

```
./serial <nx>

```

MPI version:

```
./MPI <nx>

```

CUDA version:

```
./CUDA <nx> <block_size>

```

