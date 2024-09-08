/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   RK4_ode_solver.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 17:50:41 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/05 22:20:33 by at0m             ###   ########.fr       */
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
	
    VEC_TYPE *k1, *k2, *k3, *k4;
	if (posix_memalign((void**)&k1, ALIGNMENT, 4 * sizeof(VEC_TYPE)) != 0 ||
		posix_memalign((void**)&k2, ALIGNMENT, 4 * sizeof(VEC_TYPE)) != 0 ||
		posix_memalign((void**)&k3, ALIGNMENT, 4 * sizeof(VEC_TYPE)) != 0 ||
		posix_memalign((void**)&k4, ALIGNMENT, 4 * sizeof(VEC_TYPE)) != 0) 
	{
		perror("Memory allocation failed");
		exit(1);
	}    
	__attribute__((aligned(ALIGNMENT))) VEC_TYPE temp_x[4], temp_v[4];

	printf("Computing geodesics\n");
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
				/* if (step % 1000 == 0) */
				/* { */
				/* 	printf("progress: %f\n", lambda / lambda_max); */
				/* 	printf("step: %d\n", step); */
				/* 	printf("x[0]: %f\n", VEC_EXTRACT_D(x[0])); */
				/* 	printf("x[1]: %f\n", VEC_EXTRACT_D(x[1])); */
				/* }	 */
            }
        }
	}

#ifdef AVX2
    _mm256_zeroupper();
#endif
	free(k1);
	free(k2);
	free(k3);
	free(k4);

    printf("Geodesics computed\n");
}
