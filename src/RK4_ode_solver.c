/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   RK4_ode_solver.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 17:50:41 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/01 20:04:23 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

/**
	* @brief Calculate the k1, k2, k3, k4 terms for the geodesic equation
    * @param x the position vector
    * @param v the velocity vector
    * @param k the k term
    * @param step_size the step size

    *TODO: Adaptative time step and change for a rk45 method for better precision\
            then compare k4 and k5 to get the error and adjust the step size 
*/

void geodesic_AVX(VEC_TYPE x[4], VEC_TYPE v[4], double lambda_max, VEC_TYPE christoffel[4][4][4], VEC_TYPE step_size) 
{
    int step = 0;
    double lambda = 0.0;
	
    __attribute__((aligned(ALIGNMENT))) VEC_TYPE k1[4], k2[4], k3[4], k4[4];
    __attribute__((aligned(ALIGNMENT))) VEC_TYPE temp_x[4], temp_v[4];
	#pragma omp parallel
	{	
        while (lambda < lambda_max) 
        {
            #pragma omp for simd aligned(temp_v, christoffel: ALIGNMENT) nowait
            CALCULATE_K(k1, v, christoffel)
            #pragma omp for simd aligned(temp_v, christoffel: ALIGNMENT) nowait
            UPDATE_POSITIONS(x, v, k1, step_size, temp_x, temp_v)

            #pragma omp for simd aligned(temp_v, christoffel: ALIGNMENT) nowait
            CALCULATE_K(k2, temp_v, christoffel)
            UPDATE_POSITIONS(x, v, k2, step_size, temp_x, temp_v)

            #pragma omp for simd aligned(temp_v, christoffel: ALIGNMENT) nowait
            CALCULATE_K(k3, temp_v, christoffel)
            UPDATE_POSITIONS(x, v, k3, step_size, temp_x, temp_v)

            #pragma omp for simd aligned(temp_v, christoffel: ALIGNMENT) nowait
            CALCULATE_K(k4, temp_v, christoffel)
            UPDATE_POSITIONS(x, v, k4, step_size, temp_x, temp_v)

            #pragma omp simd aligned(x, v, k1, k2, k3, k4: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) 
            {
                ALIGNED VEC_TYPE k1_term = k1[mu];
                ALIGNED VEC_TYPE k2_term = VEC_MUL_PD(VEC_SET1_PD(2.0), k2[mu]);
                ALIGNED VEC_TYPE k3_term = VEC_MUL_PD(VEC_SET1_PD(2.0), k3[mu]);
                ALIGNED VEC_TYPE k4_term = k4[mu];
                ALIGNED VEC_TYPE sum_k = VEC_ADD_PD(VEC_ADD_PD(k1_term, k4_term), VEC_ADD_PD(k2_term, k3_term));
                ALIGNED VEC_TYPE step_sum = VEC_DIV_PD(VEC_MUL_PD(step_size, sum_k), VEC_SET1_PD(6.0));

                x[mu] = VEC_ADD_PD(x[mu], step_sum);
                v[mu] = VEC_ADD_PD(v[mu], step_sum);
            }

            #pragma omp critical
            {
                lambda += VEC_EXTRACT_D(step_size);
                store_geodesic_point_AVX(x, lambda);
                step++;
            }
        }
	}

#ifdef AVX2
    _mm256_zeroupper();
#endif

    printf("Geodesics computed\n");
}
