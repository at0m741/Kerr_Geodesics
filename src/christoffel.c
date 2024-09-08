/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   christoffel.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 18:12:02 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/08 23:19:19 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
#define DELTA 1.e-5
/*
*  Calculate the Christoffel symbols
*  The Christoffel symbols are calculated using the metric tensor
*  and the inverse metric tensor
*  All calculations are done in parallel using OpenMP and AVX2 instructions 
*/

void print_m256d(VEC_TYPE x) {
    double *ptr = (double *)&x;
    for (int i = 0; i < 4; i++)
        printf("%f ", ptr[i]);
    printf("\n");
}

#define DIFFERENCE_FINITE(g1, g2, delta) ((g2 - g1) / delta)

void christoffel_AVX(VEC_TYPE g[4][4], VEC_TYPE christoffel[4][4][4], VEC_TYPE gcon[4][4])
{
    VEC_TYPE (*g_aligned)[4];
    VEC_TYPE (*christoffel_aligned)[4][4];
    
    if (posix_memalign((void **)&g_aligned, ALIGNMENT, sizeof(VEC_TYPE[4][4])) != 0) 
	{
        perror("Failed to allocate memory for g_aligned");
        exit(EXIT_FAILURE);
    }

    if (posix_memalign((void **)&christoffel_aligned, ALIGNMENT, sizeof(VEC_TYPE[4][4][4])) != 0) 
	{
        perror("Failed to allocate memory for christoffel_aligned");
        free(g_aligned);
        exit(EXIT_FAILURE);
    }
    VEC_TYPE half = VEC_SET_PD(0.5);

    memcpy(g_aligned, g, sizeof(VEC_TYPE[4][4]));
    memcpy(christoffel_aligned, christoffel, sizeof(VEC_TYPE[4][4][4]));
	
	
	#pragma omp simd
	for (int mu = 0; mu < 4; mu++) 
	{
		for (int nu = 0; nu < 4; nu++) 
		{
			for (int beta = 0; beta < 4; beta++) 
			{
				VEC_TYPE sum = VEC_SET0_PD();

				for (int sigma = 0; sigma < 4; sigma++) 
				{
					VEC_TYPE dg_sigma_nu_dbeta = DIFFERENCE_FINITE(g_aligned[sigma][nu], g_aligned[sigma + 1][nu], beta);
					VEC_TYPE dg_sigma_beta_dnu = DIFFERENCE_FINITE(g_aligned[sigma][beta], g_aligned[sigma + 1][beta], nu);
					VEC_TYPE dg_nu_beta_dsigma = DIFFERENCE_FINITE(g_aligned[nu][beta], g_aligned[nu][beta + 1], sigma);

					VEC_TYPE term1 = VEC_MUL_PD(gcon[sigma][mu], VEC_ADD_PD(dg_sigma_nu_dbeta, dg_sigma_beta_dnu));
					VEC_TYPE term2 = VEC_MUL_PD(gcon[sigma][mu], dg_nu_beta_dsigma);

					sum = VEC_ADD_PD(sum, VEC_SUB_PD(term1, term2));
				}
				sum = VEC_MUL_PD(half, sum);
				christoffel_aligned[mu][nu][beta] = sum;
			}
		}
	}

#ifdef DEBUG
    #pragma omp simd
    for (int mu = 0; mu < 4; mu++) 
        for (int nu = 0; nu < 4; nu++) 
            for (int beta = 0; beta < 4; beta++) 
                printf("Christoffel[%d][%d][%d] = %f\n", mu, nu, beta, _mm256_cvtsd_f64(christoffel_aligned[mu][nu][beta]));
	for (int mu = 0; mu < 4; mu++) 
	{
		for (int nu = 0; nu < 4; nu++) 
		{
			for (int beta = 0; beta < 4; beta++) 
			{
				if (fabs(_mm256_cvtsd_f64(christoffel_aligned[mu][nu][beta]) - _mm256_cvtsd_f64(christoffel_aligned[mu][beta][nu])) > DELTA) 
				{
					printf("Christoffel[%d][%d][%d] != Christoffel[%d][%d][%d]\n", mu, nu, beta, mu, beta, nu);
				}
			}
		}
	}
#endif

#ifdef AVX512F 
	printf("AVX512");
#endif
#ifdef AVX2
    _mm256_zeroupper();
#endif
    memcpy(christoffel, christoffel_aligned, sizeof(VEC_TYPE[4][4][4]));
    free(g_aligned);
    free(christoffel_aligned);
}



void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4])
{

#ifdef _OPENMP
    printf("OpenMP is supported and used\n");
#else
    printf("OpenMP is not supported\n");
#endif
    
	#pragma omp parallel for collapse(4)
    for (int i = 1; i < 4; i++) 
	{
        for (int j = 1; j < 4; j++) 
		{
            for (int k = 1; k < 4; k++) 
			{
                for (int l = 1; l < 4; l++) 
				{
                    double sum = 0;
                    #pragma omp simd reduction(+:sum) aligned(christoffel,g:ALIGNMENT)
                    for (int m = 0; m < 4; m++) 
                    {
                        sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] *\
													  christoffel[m][j][k] - christoffel[k][j][m] * christoffel[m][i][k]);
                        sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - \
													  christoffel[k][j][m] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                        sum += (1 / (2 * g[0][0])) * (christoffel[m][i][k] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - \
													  christoffel[m][j][k] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                    }
					riemann[i][j][k][l] = sum;
                    printf("Riemann[%d][%d][%d][%d] = %f\n", i, j, k, l, riemann[i][j][k][l]);
                }
            }
        }
    }
}
    
