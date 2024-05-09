#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>

#define MAX_POINTS 100000
#define c 299792458.0
#define G 6.67430e-11
#define M 1.9884e22
#define a 1.4
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024

#if defined(__INTEL_COMPILER) || defined(__ICC) || defined(__INTEL_LLVM_COMPILER)
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
    #define ALIGNMENT 32
    #define ARCH "x86_64"
    #define NUM_THREADS 16

#elif defined(USE_MPI)
    #define ALIGNMENT 32
    #define ARCH "x86_64"
    #define NUM_THREADS 16
    #include <omp.h>
    #include <mpi.h>
    #include "MPI_init.h"

#else
    #error "Unsupported compiler or configuration"
#endif

typedef double __attribute__((aligned(ALIGNMENT))) ldouble_a;

#if defined(__LINUX__) || defined(__linux__)
	#define sqrt_asm sqrt_asm_linux
	#define Plateform "Linux"
#elif defined(__APPLE__) || defined(__MACH__)
	#define sqrt_asm sqrt_asm_darwin
	#define Plateform "Darwin (macOS)"
#endif


inline double sqrt_asm(double n)
{
    double result;
    asm("fld %1;"
        "fsqrt;"
        "fstp %0;"
        : "=m" (result)
        : "m" (n));
    return result;
}

static inline double sqrt_asm_macos(double n)
{
	double result;
	asm("fld %1;"
		"fsqrt;"
		"fstp %0;"
		: "=m" (result)
		: "m" (n));
	return result;
}

void sincos(double x, double *sin, double *cos);
void christoffel(double g[4][4], double christoffel[4][4][4]);
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);

void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);

void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4], double step_size, void (*store_point)(double[], double));
