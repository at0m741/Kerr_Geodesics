# Geodesic Simulation in Kerr Spacetime

![singularity](https://github.com/at0m741/Kerr_Geodesics/assets/20189027/77c3b1ed-d70e-4319-81d4-32004044a585)

This project contains a C code to simulate geodesic trajectories in Kerr spacetime, which describe a rotating black hole.

# Theory
For more details on the operation and underlying theory, please see my [paper](Geodesics_trajectories_simulation (4).pdf).

# Usage
Make sure you have GCC-1X or MPICC installed on your system if you want more accurate and laborious calculations.
```bash
make x86
```
if you are using a Intel Knight Landing system using AVX512 :
```bash
make knl
```
MPI code is not ready at the moment but.. soon
```bash
make mpicc
```
Paraview is needed to plot the results

# Results and benchmarks
The results are stored in the results folder. The code will generate a .vtk file with the results of the simulation.\
The code is able to simulate the geodesic trajectory of a particle. All the code is written in C and optimized for the Intel Knight Landing architecture.\
Intel AVX2 and AVX512 intrinsics are used to optimize the code.\

For 4 millions of iterations, the code takes 2.1 seconds to compute Christoffel symbols and geodesic equations (Runge-Kutta 4th order method have been used).\
At the moment, I'm planing to rewrite a custom fprintf function to write the results in a file (writing on VTK file is the slowest part of the execution).\
```bash

# TODO
-Add the possibility to change the metric \
-Add cacheline optimizations \
-Make it work with bigger Mass and more metrics \
-Translate paper, correct potential mistakes and add bibliography\
-Add MPI code to simulate multiple particles at the same time\
-Add a custom fprintf function to write the results in a file\
-Add an adaptive step size method to optimize the code\


