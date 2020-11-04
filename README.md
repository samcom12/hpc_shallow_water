# Shallow water parallel solver 
Shallow Water Wave Equation with a finite volume solver

Requirements
====

- for serial version: `icc`
- for MPI version: `mpiicpc`
- for GPU version: CUDA® compiler

Data
====

The data files can be downloaded from [here](https://drive.switch.ch/index.php/s/7QFOGwphJun6mln). Data folder should be located in a project folder.

Build
====
To compile serial version:

```
cd serial
make
```

To compile MPI version:

```
cd MPI
make
```
	
To compile CUDA version:

```
cd CUDA      
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

