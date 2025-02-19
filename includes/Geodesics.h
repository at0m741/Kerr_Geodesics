#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>
#include <chrono>
#include <sys/time.h>


#define MAX_POINTS 100000
#define C 1.0
#define G 1.0 
#define M 1.0
#define BLOCK_SIZE 1024
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define TT 0
#define DT 0.0000005
#define max_dt 70000.0
#define ALIGNMENT 32
#define AVX2 1
#define ARCH "AVX2"
#define TOLERANCE 1e-4
#define DELTA 1e-4
#define NDIM3 3
#define DELTA3 1e-9

#include <Tensor.h>
#include <matrix.h>
#include <Metric.h>
#include <Connexion.h>
typedef struct {
    double x, y, z;
    double lambda;
} GeodesicPoint;

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


void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);
void geodesic_AVX(__m256d x[4], __m256d v[4], double lambda_max,\
				  __m256d christoffel[4][4][4], __m256d step_size);
void store_geodesic_point_AVX(__m256d x[4], double lambda);

double calculate_impact_parameter(double p_t, double p_phi, double g_tt, double g_tphi, double g_phiphi);
double calculate_emission_angle(double p_r, double p_phi, double g_rr, double g_phiphi);
double b_critique_kerr(double a, int sense);
int  compute_photon_properties(double g[4][4], double p[4]);

/* Problem specific functions */

int Riemann_tensor(const char *metric);
int Geodesics_prob();
int light_geodesics_prob(); 
int Metric_prob();
int grid_setup(); 
void generate_blackhole_image();
void generate_blackhole_shadow();
void calculate_christoffel_3D(
    double X[NDIM3],         
    double Gamma3[NDIM3][NDIM3][NDIM3]);
void calc_gamma_ij(const double X3D[3],
                   double gamma3[3][3],       
                   double gamma3_inv[3][3]);

void compute_extrinsic_curvature_stationary_3D(
    double X[3],       
    double alpha,
    double beta_cov[3],
    double Gamma3[3][3][3],
    double dbeta[3][3],
    double K[3][3]
);
void extract_3p1(
		double g[4][4],     
		double g_inv[4][4],
		double *alpha,    
		double beta_cov[3], 
		double beta_con[3],
		double gamma[3][3],
		double gamma_inv[3][3] );
void calculeBeta(double X[3], double beta_cov[3]);
void calculate_dbeta(double X[3], double dbeta[3][3]);
void compute_ricci_3d(
    const double X[3],     
    double Gamma3[3][3][3], 
    double R3[3][3]         
);
void print_ricci_tensor(double R3[3][3]);
double compute_K(double gamma_inv[3][3], double K[3][3]);
double compute_Kij_Kij(double gamma_inv[3][3], double K[3][3]);
double compute_hamiltonian_constraint(double gamma_inv[3][3], double K[3][3], double Ricci[3][3]);
void verify3x3(double matrix[3][3], double inv_matrix[3][3]);
void compute_partial_christoffel_3D(
    const double X[3],
    int m, 
    double dGamma[3][3][3], 
    double delta);
