#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>
#include <xmmintrin.h>
// #include <mpi.h>
#include <sys/time.h>
#include "/nfs/homes/ltouzali/.local/include/hdf5/include/hdf5.h"

#define MAX_POINTS 100000
#define c 299792458.0
#define G 6.67430e-11
#define M 1.0
#define a 0.935
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define TT 0
#define DT 0.0000005
#define max_dt 2.1
#define ALIGNMENT 32

#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
#define ARCH "AVX2"

#if defined(__INTEL_COMPILER) || defined(__ICC) || \
    defined(__INTEL_LLVM_COMPILER)
    #define ALIGNMENT 64
    #define ARCH "KNL"
    #define NUM_THREADS 256
    #include <omp.h>
    #include <mkl.h>
    #include <mpi.h>
    #include "MPI_init.h"

#elif defined(__MIPC) && defined(USE_MPI)
    #define ALIGNMENT 64
    #define ARCH "KNL"
    #define NUM_THREADS 256
    #include <omp.h>
    #include <mkl.h>
    #include <mpi.h>
    #include "MPI_init.h"


#elif defined(__clang__) || defined(__GNUC__)
    #define NUM_THREADS 16

#else
    #error "Unsupported compiler or configuration"
#endif

typedef double __attribute__((aligned(ALIGNMENT))) ldouble_a;

#if defined(__LINUX__) || defined(__linux__)

	#define Plateform "Linux"
#elif defined(__APPLE__) || defined(__MACH__)
	#define Plateform "Darwin (macOS)"
#endif


// double sqrt_asm(double n) {
//     double result;
//     asm("movsd %%xmm0, %1;"
//         "sqrtsd %%xmm0, %%xmm0;"
//         "movsd %0, %%xmm0;"
//         : "=m" (result)
//         : "m" (n)
//         : "%xmm0");
//     return result;
// }


void sincos(double x, double *sin, double *cos);
void christoffel(double g[4][4], double christoffel[4][4][4]);
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);
void Boyer_lindquist_coord(double *X, double *r, double *th);
void gcov(double *X, double gcov[][NDIM]);
void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);
void write_obj_file(const char *filename);
void write_hdf5(const char *filename);
void gcov(double *X, double gcov[][NDIM]);
void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4], double step_size, void (*store_point)(double[], double));


__m256d _mm256_exp_pd(__m256d x); 
__m256d _mm256_sin_pd(__m256d x); 
__m256d _mm256_cos_pd(__m256d x); 
void geodesic_AVX(__m256d x[4], __m256d v[4], float lambda_max, __m256d christoffel[4][4][4], __m256d step_size) ;
void store_geodesic_point_AVX(__m256d x[4], __m256d lambda);