/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   RK4_ode_solver.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 17:50:41 by ltouzali          #+#    #+#             */
/*   Updated: 2024/12/12 02:35:18 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"
#include <immintrin.h>

/**
	* @brief Calculate the k1, k2, k3, k4 terms for the geodesic equation
    * @param x the position vector
    * @param v the velocity vector
    * @param k the k term
    * @param step_size the step size

    *TODO: Adaptative time step and change for a rk45 method for better precision\
            then compare k4 and k5 to get the error and adjust the step size 
*/

void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], double lambda_max, VEC_TYPE christoffel[4][4][4], \
				  VEC_TYPE step_size, VEC_TYPE g[4][4]) 
{
    int step = 0;
    double lambda = 0.0;
    double energy_initial = 0.0, L_initial = 0.0;
    double threshold = 1e-6;
	double gcon[4][4];
    VEC_TYPE energy_vec = VEC_SET0_PD();
    VEC_TYPE L_vec = VEC_SET0_PD();

    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            VEC_TYPE g_mu_nu = g[mu][nu];
            VEC_TYPE v_mu = v[mu];
            VEC_TYPE v_nu = v[nu];
            VEC_TYPE term = VEC_MUL_PD(g_mu_nu, VEC_MUL_PD(v_mu, v_nu));
            energy_vec = VEC_ADD_PD(energy_vec, term);
        }
    }

    VEC_TYPE g_phi_phi = g[3][3]; 
    VEC_TYPE v_phi = v[3]; 
    L_vec = VEC_MUL_PD(g_phi_phi, v_phi);
    
    energy_initial = VEC_EXTRACT_D(energy_vec);
    L_initial = VEC_EXTRACT_D(L_vec);

    VEC_TYPE *k1, *k2, *k3, *k4;
    posix_memalign((void**)&k1, ALIGNMENT, 4 * sizeof(VEC_TYPE));
    posix_memalign((void**)&k2, ALIGNMENT, 4 * sizeof(VEC_TYPE));
    posix_memalign((void**)&k3, ALIGNMENT, 4 * sizeof(VEC_TYPE));
    posix_memalign((void**)&k4, ALIGNMENT, 4 * sizeof(VEC_TYPE));
	if (k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL)
	{
		printf("Memory allocation failed\n");
		exit(1);
	}
    __attribute__((aligned(ALIGNMENT))) VEC_TYPE temp_x[4], temp_v[4];

    while (lambda < lambda_max) {
        CALCULATE_K(k1, v, christoffel);
        UPDATE_POSITIONS(x, v, k1, step_size, temp_x, temp_v);

        CALCULATE_K(k2, temp_v, christoffel);
        UPDATE_POSITIONS(x, v, k2, step_size, temp_x, temp_v);

        CALCULATE_K(k3, temp_v, christoffel);
        UPDATE_POSITIONS(x, v, k3, step_size, temp_x, temp_v);

        CALCULATE_K(k4, temp_v, christoffel);
        UPDATE_POSITIONS(x, v, k4, step_size, temp_x, temp_v);

        for (int mu = 0; mu < 4; mu++) {
            VEC_TYPE k1_term = k1[mu];
            VEC_TYPE k2_term = VEC_MUL_PD(VEC_SET1_PD(2.0), k2[mu]);
            VEC_TYPE k3_term = VEC_MUL_PD(VEC_SET1_PD(2.0), k3[mu]);
            VEC_TYPE k4_term = k4[mu];
            VEC_TYPE sum_k = VEC_ADD_PD(VEC_ADD_PD(k1_term, k4_term), VEC_ADD_PD(k2_term, k3_term));
            VEC_TYPE step_sum = VEC_DIV_PD(VEC_MUL_PD(step_size, sum_k), VEC_SET1_PD(6.0));

            x[mu] = VEC_ADD_PD(x[mu], step_sum);
            v[mu] = VEC_ADD_PD(v[mu], step_sum);
        }

        CHECK_GEODESIC_STABILITY_AVX(x, v, g, energy_initial, L_initial, threshold);

        lambda += VEC_EXTRACT_D(step_size);
        store_geodesic_point_AVX(x, lambda);
        step++;
    }

    _mm256_zeroupper();
    free(k1);
    free(k2);
    free(k3);
    free(k4);
}
