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


# TODO

-Clean and optimize \
-Add mpi parallel code and Intel AVX2/AVX512 Intrinsics \
-Make it work with bigger Mass and more metrics \
-Translate paper, correct potential mistakes and add bibliography\


