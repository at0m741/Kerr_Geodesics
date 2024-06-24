/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   christoffel.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 18:12:02 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/25 00:19:05 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
#define __AVX2__ 1
/*
*  Calculate the Christoffel symbols
*  The Christoffel symbols are calculated using the metric tensor
*  and the inverse metric tensor
*  All calculations are done in parallel using OpenMP and AVX2 instructions 
*/
#ifdef __AVX2__
    #include <immintrin.h>


    void christoffel_AVX(__m256d g[4][4], __m256d christoffel[4][4][4]) {
        printf("using AVX2 for Christoffel Symbols\n");
        __m256d (*g_aligned)[4] = aligned_alloc(ALIGNMENT, sizeof(__m256d[4][4]));
        __m256d (*christoffel_aligned)[4][4] = aligned_alloc(ALIGNMENT, sizeof(__m256d[4][4][4]));
        memcpy(g_aligned, g, sizeof(__m256d[4][4]));
        memcpy(christoffel_aligned, christoffel, sizeof(__m256d[4][4][4]));

        // Precompute frequently used products and broadcast them outside the innermost loop
        __m256d half = _mm256_set1_pd(0.5);

        for (int mu = 0; mu < 4; mu++) {
            for (int beta = 0; beta < 4; beta++) {
                for (int nu = 0; nu < 4; nu++) {
                    __m256d sum = _mm256_setzero_pd();

                    for (int sigma = 0; sigma < 4; sigma++) {
                        // Prefetch for next iteration of sigma if within bounds
                        if (sigma + 1 < 4) {
                            _mm_prefetch((const char *)&g_aligned[mu][sigma + 1], _MM_HINT_T0);
                            _mm_prefetch((const char *)&christoffel_aligned[mu][sigma + 1], _MM_HINT_T1);
                        }

                        __m256d g_mu_sigma = _mm256_broadcast_sd((double *)&g_aligned[mu][sigma]);
                        __m256d g_sigma_beta = _mm256_broadcast_sd((double *)&g_aligned[sigma][beta]);
                        __m256d g_beta_sigma = _mm256_broadcast_sd((double *)&g_aligned[beta][sigma]);
                        __m256d g_beta_nu = _mm256_broadcast_sd((double *)&g_aligned[beta][nu]);

                        __m256d term1 = _mm256_mul_pd(g_mu_sigma, g_sigma_beta);
                        __m256d term2 = _mm256_mul_pd(g_mu_sigma, g_beta_sigma);
                        __m256d term3 = _mm256_mul_pd(g_mu_sigma, g_beta_nu);

                        sum = _mm256_add_pd(sum, _mm256_mul_pd(half, \
                            _mm256_sub_pd(_mm256_add_pd(term1, term2), term3)));
                    }

                    christoffel_aligned[mu][beta][nu] = sum;
                }
            }
        }
        _mm256_zeroupper();
        _mm_empty();
        memcpy(christoffel, christoffel_aligned, sizeof(__m256d[4][4][4]));
        free(g_aligned);
        free(christoffel_aligned);
    }

#else
    #pragma omp declare simd
    void christoffel(double g[4][4], double christoffel[4][4][4])
    {
        #ifdef USE_MPI
            printf("MPI is defined for Christoffel\n");
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        #elif _OPENMP
            printf("OpenMP is supported and used\n");
        #else
            printf("OpenMP is not supported\n");
        #endif

        double (*g_aligned)[4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4]));
        double (*christoffel_aligned)[4][4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4][4]));
        memcpy(g_aligned, g, sizeof(double[4][4]));
        memcpy(christoffel_aligned, christoffel, sizeof(double[4][4][4]));

__m256d _mm256_exp_pd(__m256d x); 
// __m256d _mm256_sin_pd(__m256d x); 
// __m256d _mm256_cos_pd(__m256d x); 
        #pragma omp parallel for collapse(3)
        #pragma vector aligned
        for (int mu = 0; mu < 4; mu++) {
            for (int beta = 0; beta < 4; beta++) {
                for (int nu = 0; nu < 4; nu++) {
                    double sum = 0;
                    #pragma omp simd reduction(+:sum) aligned(g_aligned:ALIGNMENT)
                    #pragma vector aligned
                    for (int sigma = 0; sigma < 4; sigma++) {
                        sum += 0.5 * (g_aligned[mu][sigma] * (g_aligned[sigma][beta] + g_aligned[beta][sigma] - g_aligned[beta][nu]));                    
                    }
                    christoffel_aligned[mu][beta][nu] = sum;
                }
            }
        }

        memcpy(christoffel, christoffel_aligned, sizeof(double[4][4][4]));
        free(g_aligned);
        free(christoffel_aligned);
    }


    #pragma omp declare simd
    void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4])
    {
        #ifdef _OPENMP
            printf("OpenMP is supported and used\n");
        #else
            printf("OpenMP is not supported\n");
        #endif
        #pragma omp parallel for collapse(4)
        for (int i = 1; i < 4; i++) {
            for (int j = 1; j < 4; j++) {
                for (int k = 1; k < 4; k++) {
                    for (int l = 1; l < 4; l++) {
                        double sum = 0;
                        #pragma omp simd reduction(+:sum) aligned(christoffel,g:ALIGNMENT)
                        for (int m = 0; m < 4; m++) 
                        {
                            sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] * christoffel[m][j][k] - christoffel[k][j][m] * christoffel[m][i][k]);
                            sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - christoffel[k][j][m] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                            sum += (1 / (2 * g[0][0])) * (christoffel[m][i][k] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - christoffel[m][j][k] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                        }
                        riemann[i][j][k][l] = sum;
                        printf("Riemann[%d][%d][%d][%d] = %f\n", i, j, k, l, riemann[i][j][k][l]);
                    }
                }
            }
        }
    }

#endif