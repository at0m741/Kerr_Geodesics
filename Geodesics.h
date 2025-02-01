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

#define MAX_POINTS 100000
#define c 299792458.0
#define G 6.67430e-11
#define M 5.0
#define a 0.2
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define TT 0
#define DT 0.0000005
#define max_dt 100840.0
#define ALIGNMENT 32
#define AVX2 1
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

#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)


#ifdef  AVX2
    #define VEC_TYPE __m256d
    #define VEC_MUL_PD(a, b) _mm256_mul_pd(a, b)
    #define VEC_SUB_PD(a, b) _mm256_sub_pd(a, b)
    #define VEC_ADD_PD(a, b) _mm256_add_pd(a, b)
    #define VEC_DIV_PD(a, b) _mm256_div_pd(a, b)
    #define VEC_SET1_PD(a)   _mm256_set1_pd(a)
    #define VEC_LOAD_PD(a)   _mm256_load_pd(a)
    #define VEC_STORE_PD(a, b) _mm256_store_pd(a, b)
    #define VEC_CVTSD_F64(a) _mm256_cvtsd_f64(a)
    #define CALCULATE_K_AVX2(k, src_v) \
    for (int mu = 0; mu < 4; mu++) { \
        k[mu] = src_v[mu]; \
        for (int alpha = 0; alpha < 4; alpha++) { \
            __m256d product1 = _mm256_mul_pd(christoffel[mu][alpha][0], \
                                _mm256_mul_pd(src_v[alpha], src_v[0])); \
            __m256d product2 = _mm256_mul_pd(christoffel[mu][alpha][1], \
                                _mm256_mul_pd(src_v[alpha], src_v[1])); \
            __m256d product3 = _mm256_mul_pd(christoffel[mu][alpha][2], \
                                _mm256_mul_pd(src_v[alpha], src_v[2])); \
            __m256d product4 = _mm256_mul_pd(christoffel[mu][alpha][3], \
                                _mm256_mul_pd(src_v[alpha], src_v[3])); \
            k[mu] = _mm256_sub_pd(k[mu], product1); \
            k[mu] = _mm256_sub_pd(k[mu], product2); \
            k[mu] = _mm256_sub_pd(k[mu], product3); \
            k[mu] = _mm256_sub_pd(k[mu], product4); \
        } \
    }

    #define UPDATE_POSITIONS_AVX2(x, v, k, step) \
        for (int mu = 0; mu < 4; mu++) { \
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step, k[mu])); \
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step, k[mu])); \
        }


#else
    #error "AVX2 or AVX512F support required"
#endif

void sincos(double x, double *sin, double *cos);
void calculate_christoffel(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]); 
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);
void Boyer_lindquist_coord(double *X, double *r, double *th);
void gcov(double *X, double gcov[][NDIM]);
void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);
void write_obj_file(const char *filename);
void write_hdf5(const char *filename);
void gcov(double *X, double gcov[][NDIM]);

void calculate_metric(double x[4], double g[4][4], double g_inv[4][4]);
void verify_metric(double g[4][4], double g_inv[4][4]);

__m256d _mm256_exp_pd(__m256d x); 
__m256d _mm256_sin_pd(__m256d x); 
__m256d _mm256_cos_pd(__m256d x); 
void geodesic_AVX(__m256d x[4], __m256d v[4], double lambda_max,\
				  __m256d christoffel[4][4][4], __m256d step_size);
void store_geodesic_point_AVX(__m256d x[4], double lambda);
