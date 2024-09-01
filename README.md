# Geodesic Simulation in Kerr Spacetime

![singularity](https://github.com/at0m741/Kerr_Geodesics/assets/20189027/77c3b1ed-d70e-4319-81d4-32004044a585)

This project contains a C code to simulate geodesic trajectories in Kerr spacetime, which describe a rotating black hole.

# Theory
For more details on the operation and underlying theory, please see my [paper]([Paper](https://github.com/at0m741/Kerr_Geodesics/blob/main/Geodesics_trajectories_simulation%20(4).pdf)).

# Usage
Make sure you have GCC-1X installed on your system.
```bash
make avx2
```
If you are using a Intel Knight Landing or anyother system that can use AVX512 :
```bash
make avx512
```
Paraview is needed to plot the results


# TODO

-Make it work with bigger Mass and more metrics \
-Paper correct potential mistakes and add bibliography\


