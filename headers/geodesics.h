/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   geodesics.h                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 15:03:15 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/27 15:41:14 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <time.h>
#include <xmmintrin.h>
#include <sys/time.h>

#define MAX_POINTS 100000
#define c 299792458.0
#define G 6.67430e-11
#define M 1.0
#define a 0.935
#define BLOCK_SIZE 4
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define TT 0
#define DT 0.0000005
#define max_dt 4.0
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


#ifdef __AVX512F__

    #define VEC_TYPE __m512d
    #define VEC_MUL_PD(a, b) _mm512_mul_pd(a, b)
    
    #define VEC_SUB_PD(a, b) _mm512_sub_pd(a, b)
    #define VEC_ADD_PD(a, b) _mm512_add_pd(a, b)
    #define VEC_DIV_PD(a, b) _mm512_div_pd(a, b)
    #define VEC_SET1_PD(a)   _mm512_set1_pd(a)
    #define VEC_LOAD_PD(a)   _mm512_load_pd(a)
    #define VEC_STORE_PD(a, b) _mm512_store_pd(a, b)
    #define VEC_CVTSD_F64(a) _mm512_cvtsd_f64(a)
    #define CALCULATE_K_AVX2(k, src_v) \
        for (int mu = 0; mu < 4; mu++) { \
            k[mu] = src_v[mu]; \
            for (int alpha = 0; alpha < 4; alpha++) { \
                __m512d src_v_alpha_src_v0 = _mm512_mul_pd(src_v[alpha], src_v[0]); \
                __m512d src_v_alpha_src_v1 = _mm512_mul_pd(src_v[alpha], src_v[1]); \
                __m512d src_v_alpha_src_v2 = _mm512_mul_pd(src_v[alpha], src_v[2]); \
                __m512d src_v_alpha_src_v3 = _mm512_mul_pd(src_v[alpha], src_v[3]); \
                __m512d product1 = _mm512_mul_pd(christoffel[mu][alpha][0], src_v_alpha_src_v0); \
                __m512d product2 = _mm512_mul_pd(christoffel[mu][alpha][1], src_v_alpha_src_v1); \
                __m512d product3 = _mm512_mul_pd(christoffel[mu][alpha][2], src_v_alpha_src_v2); \
                __m512d product4 = _mm512_mul_pd(christoffel[mu][alpha][3], src_v_alpha_src_v3); \
                k[mu] = _mm512_sub_pd(k[mu], product1); \
                k[mu] = _mm512_sub_pd(k[mu], product2); \
                k[mu] = _mm512_sub_pd(k[mu], product3); \
                k[mu] = _mm512_sub_pd(k[mu], product4); \
            } \
        }


    #define UPDATE_POSITIONS_AVX2(x, v, k, step) \
        for (int mu = 0; mu < 4; mu++) { \
            temp_x[mu] = _mm512_add_pd(x[mu], _mm512_mul_pd(step, k[mu])); \
            temp_v[mu] = _mm512_add_pd(v[mu], _mm512_mul_pd(step, k[mu])); \
        }


#elif AVX2

    /* 
    *  Use of the AVX2 instruction set for the geodesic calculation
    *  The AVX2 instruction set is used to calculate the geodesic
    * @param x[4] - the position of the geodesic
    * @param v[4] - the velocity of the geodesic
    * @param lambda_max - the maximum value of the affine parameter
    * @param christoffel[4][4][4] - the Christoffel symbols
    * @param step_size - the step size for the geodesic calculation
    * the equation is : \frac{d^2 x^\mu}{d\lambda^2} + 
                        \Gamma^\mu_{\alpha \beta} 
                        \frac{d x^\alpha}{d\lambda}
                        \frac{d x^\beta}{d\lambda} = 0
    */


    #define ALIGNED_32 __attribute__((aligned(32)))
    #define LOOP_OVER_3_INDICES for (int mu = 0; mu < 4; mu++) \
                            for (int beta = 0; beta < 4; beta++) \
                            for (int nu = 0; nu < 4; nu++)

    #define LOOP_OVER_SIGMA for (int sigma = 0; sigma < 4; sigma++)
    /*
    *  AVX2 instruction set for the geodesic calculation
    *  The AVX2 instruction set is used to calculate the geodesic 
    */
    #define VEC_TYPE __m256d
    #define CHECK_NAN_PD(v) _mm256_cmp_pd(v, v, _CMP_UNORD_Q)
    #define HANDLE_NAN(v, replacement) _mm256_blendv_pd(v, replacement, CHECK_NAN_PD(v))   
    #define VEC_MUL_PD(a, b) _mm256_mul_pd(a, b)
    #define VEC_SUB_PD(a, b) _mm256_sub_pd(a, b)
    #define VEC_ADD_PD(a, b) _mm256_add_pd(a, b)
    #define VEC_DIV_PD(a, b) _mm256_div_pd(a, b)
    #define VEC_SET1_PD(a)   _mm256_set1_pd(a)
    #define VEC_LOAD_PD(a)   _mm256_load_pd(a)
    #define VEC_STORE_PD(a, b) _mm256_store_pd(a, b)
    #define VEC_CVTSD_F64(a) _mm256_cvtsd_f64(a)
    /*
    *  Defined loop for the calculation of the Christoffel symbols
    *  The loop is used to calculate the Geodesics equation
    *  using the Christoffel symbols
    */
    #define CALCULATE_K_AVX2(k, src_v) \
    for (int mu = 0; mu < 4; mu++) { \
        k[mu] = src_v[mu]; \
        for (int alpha_block = 0; alpha_block < 4; alpha_block += BLOCK_SIZE) { \
            for (int alpha = alpha_block; alpha < alpha_block + BLOCK_SIZE && alpha < 4; alpha++) { \
                _mm_prefetch((const char *)&christoffel[mu][alpha][0], _MM_HINT_NTA); \
                for (int beta = 0; beta < 4; beta++) { \
                    __m256d src_v_alpha_src_v = _mm256_mul_pd(src_v[alpha], src_v[beta]); \
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], src_v_alpha_src_v); \
                    k[mu] = _mm256_sub_pd(k[mu], product); \
                } \
            } \
        } \
    }

    /* 
    *  Update the positions of the geodesic using the AVX2 instruction set
    */
    #define UPDATE_POSITIONS_AVX2(x, v, k, step) \
        for (int mu = 0; mu < 4; mu++) { \
            _mm_prefetch((const char *)&temp_x[mu], _MM_HINT_NTA); \
            _mm_prefetch((const char *)&temp_v[mu], _MM_HINT_NTA); \
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step, k[mu])); \
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step, k[mu])); \
        }


#else
    #error "AVX2 or AVX512F support required"
#endif

void sincos(double x, double *sin, double *cos);
/** 
    * @brief riemann - Calculate the Riemann tensor (optional but can be called from main after the Christoffel symbols are calculated)
    * @param g - the metric tensor
    * @param christoffel - the Christoffel symbols
    * @param riemann - the Riemann tensor
*/
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4]);
void Boyer_lindquist_coord(double *X, double *r, double *th);
void gcov(double *X, double gcov[][NDIM]);
void write_vtk_file(const char *filename);
void store_geodesic_point(double x[4], double lambda);
void write_obj_file(const char *filename);
void write_hdf5(const char *filename);
void gcov(double *X, double gcov[][NDIM]);
void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4],\
              double step_size, void (*store_point)(double[], double));
/* 
*  AVX2 and AVX512F intrinsics for the exponential, sine and cosine functions
*/
__m256d _mm256_exp_pd(__m256d x); 
__m256d _mm256_sin_pd(__m256d x); 
__m256d _mm256_cos_pd(__m256d x); 
/*
*  AVX2 and AVX512F intrinsics for the geodesic calculation
* The AVX2 and AVX512F instruction sets are used to calculate the geodesic
* same for Christoffel symbols and step size
*/
void christoffel_AVX(__m256d g[4][4], __m256d christoffel[4][4][4]);
void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], double lambda_max, \
                  VEC_TYPE christoffel[4][4][4], VEC_TYPE step_size);
void store_geodesic_point_AVX(__m256d x[4], double lambda);