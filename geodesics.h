#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

#define MAX_POINTS 100000
#define c 299792458
#define G 6.67430e-11
#define M 1.9884e26
#define a 2.9
#define BLOCK_SIZE 1000

#if  defined(__INTEL_COMPILER) || defined(__ICC)
	#define ALIGNMENT 64
	#define ARCH "KNL"
#elif defined(__clang__) || defined(__GNUC__)
	#define ALIGNMENT 32
	#define ARCH "x86_64"
#else
	0
#endif
typedef double __attribute__((aligned(ALIGNMENT))) ldouble_a;
#define BUFFER_SIZE 1024
#if defined(__LINUX__) || defined(__linux__)
	#define sqrt_asm sqrt_asm_linux
	#define Plateform "Linux"
#elif defined(__APPLE__) || defined(__MACH__)
	#define sqrt_asm sqrt_asm_darwin
	#define Plateform "Darwin (macOS)"
#endif

static inline double sqrt_asm(double n)
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


void christoffel(double g[4][4], double christoffel[4][4][4]);
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);

void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);

void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4], double step_size, void (*store_point)(double[], double));
