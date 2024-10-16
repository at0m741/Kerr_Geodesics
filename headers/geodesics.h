/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   geodesics.h                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 15:03:15 by ltouzali          #+#    #+#             */
/*   Updated: 2024/10/10 00:14:31 by babonnet         ###   ########.fr       */
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
#include <fcntl.h>

#undef a
#undef G
#undef M
#undef c
#undef BLOCK_SIZE
#undef BUFFER_SIZE
#undef SMALL
#undef NDIM
#undef TT
#undef DT
#undef max_dt

#define MAX_POINTS 100000
#define c 1.0 
#define G 1.0
#define M 1.0
#define a 0.9375
#define BLOCK_SIZE 4
#define BUFFER_SIZE 1024
#define SMALL 1.e-40
#define NDIM 4
#define TT 0
#define DT 0.0000002
#define max_dt 7.6

#define TOLERANCE 1e-6

#ifdef AVX512F
	#pragma message "AVX512"
    #define ARCH				"AVX512"
    #define ALIGNMENT			64
    #define VEC_TYPE			__m512d
    #define VEC_SET0_PD()		_mm512_setzero_pd()
    #define VEC_SET_PD(a)		_mm512_set1_pd(a)
    #define VEC_MUL_PD(a, b)	_mm512_mul_pd(a, b)
    #define VEC_SUB_PD(a, b)	_mm512_sub_pd(a, b)
    #define VEC_ADD_PD(a, b)	_mm512_add_pd(a, b)
    #define VEC_DIV_PD(a, b)	_mm512_div_pd(a, b)
    #define VEC_SET1_PD(a)		_mm512_set1_pd(a)
    #define VEC_LOAD_PD(a)		_mm512_load_pd(a)
	#define VEC_LOADU_PD(a)		_mm512_loadu_pd(a)
    #define VEC_STORE_PD(a, b)	_mm512_store_pd(a, b)
    #define VEC_CVTSD_F64(a)	_mm512_cvtsd_f64(a)
	#define VEC_EXTRACT_D(a)	_mm512_cvtsd_f64(a)
#endif

#ifdef AVX2
	#pragma message "AVX2"
    #define ARCH				"AVX2"
    #define ALIGNMENT			32
    #define VEC_TYPE			__m256d
    #define VEC_SET0_PD()		_mm256_setzero_pd()
    #define VEC_SET_PD(a)		_mm256_set1_pd(a)
    #define VEC_MUL_PD(a, b)	_mm256_mul_pd(a, b)
    #define VEC_SUB_PD(a, b)	_mm256_sub_pd(a, b)
    #define VEC_ADD_PD(a, b)	_mm256_add_pd(a, b)
    #define VEC_DIV_PD(a, b)	_mm256_div_pd(a, b)
    #define VEC_SET1_PD(a)		_mm256_set1_pd(a)
    #define VEC_LOAD_PD(a)		_mm256_load_pd(a)
	#define VEC_LOADU_PD(a)		_mm256_loadu_pd(a)
    #define VEC_STORE_PD(a, b)	_mm256_store_pd(a, b)
    #define VEC_CVTSD_F64(a)	_mm256_cvtsd_f64(a)
	#define VEC_EXTRACT_D(a)	_mm256_cvtsd_f64(a)
#endif

typedef double __attribute__((aligned(ALIGNMENT))) ldouble_a;

#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)

#define ALIGNED __attribute__((aligned(ALIGNMENT)))

#if defined(__LINUX__) || defined(__linux__)
	#define Plateform "Linux"
#elif defined(__APPLE__) || defined(__MACH__)
	#define Plateform "Darwin (macOS)"
#endif



#define CHECK_GEODESIC_STABILITY_AVX(x, v, g, energy_initial, L_initial, threshold) { \
    VEC_TYPE energy_vec = VEC_SET0_PD(); \
    VEC_TYPE L_vec = VEC_SET0_PD(); \
    VEC_TYPE norm_vec = VEC_SET0_PD(); \
    \
    for (int mu = 0; mu < 4; mu++) { \
        for (int nu = 0; nu < 4; nu++) { \
            VEC_TYPE g_mu_nu = g[mu][nu]; \
            VEC_TYPE v_mu = v[mu]; \
            VEC_TYPE v_nu = v[nu]; \
            VEC_TYPE energy_term = VEC_MUL_PD(g_mu_nu, VEC_MUL_PD(v_mu, v_nu)); \
            energy_vec = VEC_SUB_PD(energy_vec, energy_term); \
            norm_vec = VEC_ADD_PD(norm_vec, VEC_MUL_PD(g_mu_nu, VEC_MUL_PD(v_mu, v_nu))); \
        } \
    } \
    \
    VEC_TYPE g_phi_phi = g[3][3]; \
    VEC_TYPE v_phi = v[3]; \
    L_vec = VEC_MUL_PD(g_phi_phi, v_phi); \
    \
    double energy = VEC_EXTRACT_D(energy_vec); \
    double L = VEC_EXTRACT_D(L_vec); \
    double norm = VEC_EXTRACT_D(norm_vec); \
    \
  }


#define CALCULATE_K(k, src_v, christoffel) \
    for (int mu = 0; mu < 4; mu++) { \
        k[mu] = src_v[mu]; \
        for (int alpha = 0; alpha < 4; alpha++) { \
            VEC_TYPE src_v_alpha_src_v_0 = VEC_MUL_PD(src_v[alpha], src_v[0]); \
            VEC_TYPE product_0 = VEC_MUL_PD(christoffel[mu][alpha][0], src_v_alpha_src_v_0); \
            k[mu] = VEC_SUB_PD(k[mu], product_0); \
            VEC_TYPE src_v_alpha_src_v_1 = VEC_MUL_PD(src_v[alpha], src_v[1]); \
            VEC_TYPE product_1 = VEC_MUL_PD(christoffel[mu][alpha][1], src_v_alpha_src_v_1); \
            k[mu] = VEC_SUB_PD(k[mu], product_1); \
            VEC_TYPE src_v_alpha_src_v_2 = VEC_MUL_PD(src_v[alpha], src_v[2]); \
            VEC_TYPE product_2 = VEC_MUL_PD(christoffel[mu][alpha][2], src_v_alpha_src_v_2); \
            k[mu] = VEC_SUB_PD(k[mu], product_2); \
            VEC_TYPE src_v_alpha_src_v_3 = VEC_MUL_PD(src_v[alpha], src_v[3]); \
            VEC_TYPE product_3 = VEC_MUL_PD(christoffel[mu][alpha][3], src_v_alpha_src_v_3); \
            k[mu] = VEC_SUB_PD(k[mu], product_3); \
        } \
    }

#define UPDATE_POSITIONS(x, v, k, step, temp_x, temp_v) \
    for (int mu = 0; mu < 4; mu++) { \
        temp_x[mu] = VEC_ADD_PD(x[mu], VEC_MUL_PD(step, k[mu])); \
        temp_v[mu] = VEC_ADD_PD(v[mu], VEC_MUL_PD(step, k[mu])); \
    }
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
void gcon(double r, double th, double gcon[][NDIM]);
void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4],\
              double step_size, void (*store_point)(double[], double));
double calculate_angular_momentum(double v[4], double g[4][4]);
/* 
	* AVX2 and AVX512F intrinsics for the exponential, sine and cosine functions
*/

VEC_TYPE _mm256_exp_pd(VEC_TYPE x); 
VEC_TYPE _mm256_sin_pd(VEC_TYPE x); 
VEC_TYPE _mm256_cos_pd(VEC_TYPE x); 
void sincos(double x, double *sin, double *cos);

/*
	* AVX2 and AVX512F intrinsics for the geodesic calculation
	* The AVX2 and AVX512F instruction sets are used to calculate the geodesic
	* same for Christoffel symbols and step size
*/

void christoffel_symbols(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM]); 
void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], double lambda_max, \
                  VEC_TYPE christoffel[4][4][4], VEC_TYPE step_size, VEC_TYPE g[4][4]);
void store_geodesic_point_AVX(VEC_TYPE x[4], double lambda);
void invert_metric(double gcov[][NDIM], double gcon[][NDIM]);
int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]); 
void Boyer_lindquist_coord(double *X, double *r, double *th);
void verify_metric(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
