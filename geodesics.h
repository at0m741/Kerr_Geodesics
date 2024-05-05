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
#define M 1.9884e22
#define a 2.9
#define ALIGNMENT 32
#define BLOCK_SIZE 1000
typedef long double __attribute__((aligned(32))) ldouble_a32;
typedef double __attribute__((aligned(32))) double_a32;

static inline double_a32 sqrt_asm(double n)
{
    double result;
    asm("fld %1;"
        "fsqrt;"
        "fstp %0;"
        : "=m" (result)
        : "m" (n));
    return result;
}

void christoffel(long double g[4][4], long double christoffel[4][4][4]);
void riemann(long double g[4][4], long double christoffel[4][4][4], long double riemann[4][4][4][4]);


void write_vtk_file(const char *filename);
void store_geodesic_point(long double x[4], long double lambda);

void geodesic(long double x[4], long double v[4], long double lambda_max, long double christoffel[4][4][4], long double step_size, void (*store_point)(long double[], long double));