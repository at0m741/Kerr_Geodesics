/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   christoffel.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 18:12:02 by ltouzali          #+#    #+#             */
/*   Updated: 2024/12/12 02:48:23 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;

/*
*  Calculate the Christoffel symbols
*  The Christoffel symbols are calculated using the metric tensor
*  and the inverse metric tensor
*  All calculations are done in parallel using OpenMP and AVX2 instructions 
*/


void christoffel_symbols(double X[NDIM], double h, double gamma[NDIM][NDIM][NDIM]) {
    double _gcov[NDIM][NDIM];
    double _gcon[NDIM][NDIM];
    double dgcov_mu[NDIM][NDIM][NDIM];

    double gcov_forward[NDIM][NDIM], gcov_backward[NDIM][NDIM];

    gcov(X, _gcov);
	inverse_matrix(_gcov, _gcon);
	verify_metric(_gcov, _gcon);
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int sigma = 0; sigma < NDIM; sigma++) {
                dgcov_mu[mu][nu][sigma] = DERIV_CENTREE(gcov, X, h, gcov_forward, gcov_backward, mu, sigma);
                printf("dgcov_mu[%d][%d][%d] = %f\n", mu, nu, sigma, dgcov_mu[mu][nu][sigma]);
            }
        }
    }

    VEC_TYPE vec_gcon, vec_dgcov_mu_nu_sigma, vec_dgcov_nu_mu_sigma, vec_dgcov_sigma_mu_nu;
    VEC_TYPE vec_sum, vec_half, vec_gamma_lambda_mu_nu;

    vec_half = VEC_SET1_PD(0.5);

    for (int lambda = 0; lambda < NDIM; lambda++) {
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = 0; nu < NDIM; nu++) {
                vec_sum = VEC_SET0_PD();

                for (int sigma = 0; sigma < NDIM; sigma++) {
                    vec_gcon = VEC_SET1_PD(_gcon[lambda][sigma]);

                    vec_dgcov_mu_nu_sigma   = VEC_SET1_PD(dgcov_mu[mu][nu][sigma]);
                    vec_dgcov_nu_mu_sigma   = VEC_SET1_PD(dgcov_mu[nu][mu][sigma]);
                    vec_dgcov_sigma_mu_nu   = VEC_SET1_PD(dgcov_mu[sigma][mu][nu]);

                    VEC_TYPE vec_term = VEC_ADD_PD(vec_dgcov_mu_nu_sigma, vec_dgcov_nu_mu_sigma);
                    vec_term = VEC_SUB_PD(vec_term, vec_dgcov_sigma_mu_nu);
                    vec_term = VEC_MUL_PD(vec_gcon, vec_term);

                    vec_sum = VEC_ADD_PD(vec_sum, vec_term);
                }

                vec_gamma_lambda_mu_nu = VEC_MUL_PD(vec_sum, vec_half);
                gamma[lambda][mu][nu] = VEC_EXTRACT_D(vec_gamma_lambda_mu_nu);
            }
        }
    }

    printf("Christoffel symbols calculated\n");
    printf("result : \n");
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                printf("gamma[%d][%d][%d] = %f\n", i, j, k, gamma[i][j][k]);
            }
        }
    }
}
