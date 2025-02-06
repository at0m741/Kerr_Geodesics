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
#define M 1.0
#define a 0.0
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
#define TOLERANCE 1e-4
#define DELTA 1e-4





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

void calculate_christoffel(double X[NDIM], double h, \
							double gamma[NDIM][NDIM][NDIM],
							double g[NDIM][NDIM],
							double g_inv[NDIM][NDIM], const char *metric); 
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);
void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);
void calculate_metric(double x[4], double g[4][4], double g_inv[4][4]);
void verify_metric(double g[4][4], double g_inv[4][4]);
void minkowski_metric(double g[NDIM][NDIM], double g_inv[NDIM][NDIM]); 
void geodesic_AVX(__m256d x[4], __m256d v[4], double lambda_max,\
				  __m256d christoffel[4][4][4], __m256d step_size);
void store_geodesic_point_AVX(__m256d x[4], double lambda);

void check_symmetry_christoffel(double gamma[NDIM][NDIM][NDIM]); 
double determinant3x3(double mat[3][3]);
double determinant4x4(double mat[4][4]);
void cofactor(double mat[NDIM][NDIM], double cofactorMat[NDIM][NDIM]);
void transpose(double mat[NDIM][NDIM], double transposed[NDIM][NDIM]);
int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]);
void check_inverse(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
void print_matrix(const char* name, double mat[NDIM][NDIM]);

void check_symmetry_christoffel(double gamma[NDIM][NDIM][NDIM]); 
void initialize_riemann_tensor(double R[NDIM][NDIM][NDIM][NDIM]);
void print_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM]);
void print_christoffel(double Gamma[NDIM][NDIM][NDIM]);
void print_christoffel_matrix(double gamma[NDIM][NDIM][NDIM]);


void calculate_Gamma_at_offset(double X[NDIM], int direction, 
						double offset, double delta,
						double gcov[NDIM][NDIM], 
						double gcon[NDIM][NDIM], 
						double Gamma_slice[NDIM][NDIM][NDIM], 
						const char *metric_type);
void calculate_riemann(double Gamma[NDIM][NDIM][NDIM], 
                       double Gamma_plus_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_minus_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM],
                       double Riemann[NDIM][NDIM][NDIM][NDIM], 
                       double h);
 
void check_riemann_symmetries(double Riemann[NDIM][NDIM][NDIM][NDIM], double tolerance);
void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM], double g_inv[NDIM][NDIM]);

