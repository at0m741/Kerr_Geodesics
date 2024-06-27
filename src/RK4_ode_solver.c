/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   RK4_ode_solver.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/23 17:50:41 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/27 15:42:36 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

#ifdef AVX2
/**
 * @brief Calculate the k1, k2, k3, k4 terms for the geodesic equation
    * @param x the position vector
    * @param v the velocity vector
    * @param k the k term
    * @param step_size the step size

    *TODO: Adaptative time step and change for a rk45 method for better precision\
            then compare k4 and k5 to get the error and adjust the step size 
*/

    double calculate_max_speed(__m256d christoffel[4][4][4], __m256d x[4], __m256d v[4]) {
        // Un exemple simple qui retourne toujours 1.0 (vitesse de la lumière en unités naturelles)
        return 1.0;
    }

    void adjust_step_size(__m256d christoffel[4][4][4], __m256d x[4], __m256d v[4], __m256d step_size, double delta_x, double safety_factor) {
        double max_speed = calculate_max_speed(christoffel, x, v);
        double new_step_size = safety_factor * delta_x / max_speed;
        if (new_step_size < _mm256_cvtsd_f64(step_size)) {
            step_size = _mm256_set1_pd(new_step_size);
        }
    }

    void geodesic_AVX(__m256d x[4], __m256d v[4], double lambda_max, __m256d christoffel[4][4][4], __m256d step_size) 
    {
        int step = 0;
        double lambda = 0.0;
        __attribute__((aligned(32))) __m256d k1[4], k2[4], k3[4], k4[4], k5[4],\
                                            temp_x[4], temp_v[4];
        printf("using AVX2 for geodesic\n");

        while (lambda < lambda_max) {
            adjust_step_size(christoffel, x, v, step_size, 0.01, 0.9);
            #pragma omp simd aligned(temp_v, christoffel: 32)
            CALCULATE_K_AVX2(k1, v)
            #pragma omp simd aligned(temp_v, christoffel: 32)
            UPDATE_POSITIONS_AVX2(x, v, k1, step_size)
            #pragma omp simd aligned(temp_v, christoffel: 32)
            CALCULATE_K_AVX2(k2, temp_v)
            UPDATE_POSITIONS_AVX2(x, v, k2, step_size)
            #pragma omp simd aligned(temp_v, christoffel: 32) 
            CALCULATE_K_AVX2(k3, temp_v)
            UPDATE_POSITIONS_AVX2(x, v, k3, step_size)
            #pragma omp simd aligned(temp_v, christoffel: 32)
            CALCULATE_K_AVX2(k4, temp_v)
            UPDATE_POSITIONS_AVX2(x, v, k4, step_size)
            /*
             *update the position and velocity vectors 
            */
            #pragma omp prefetch
            #pragma omp simd aligned(x, v, k1, k2, k3, k4: 32)
            for (int mu = 0; mu < 4; mu++) {
                _mm_prefetch((const char *)&x[mu], _MM_HINT_T0);
                _mm_prefetch((const char *)&v[mu], _MM_HINT_T0);
                ALIGNED_32 __m256d k1_term = k1[mu];
                ALIGNED_32 __m256d k2_term = _mm256_mul_pd(_mm256_set1_pd(2.0), k2[mu]);
                ALIGNED_32 __m256d k3_term = _mm256_mul_pd(_mm256_set1_pd(2.0), k3[mu]);
                ALIGNED_32 __m256d k4_term = k4[mu];
                ALIGNED_32 __m256d sum_k = _mm256_add_pd(_mm256_add_pd(k1_term, k4_term), \
                                _mm256_add_pd(k2_term, k3_term));
                ALIGNED_32 __m256d step_sum = _mm256_div_pd(_mm256_mul_pd(step_size, sum_k), \
                                _mm256_set1_pd(6.0));

                x[mu] = _mm256_add_pd(x[mu], step_sum);
                v[mu] = _mm256_add_pd(v[mu], step_sum);
            }

            lambda += _mm256_cvtsd_f64(step_size);
            store_geodesic_point_AVX(x, lambda);
        }
        printf("Geodesics computed\n");
    }
#elif __AVX512F__
void geodesic_AVX(__m512d x[4], __m512d v[4], double lambda_max, __m512d christoffel[4][4][4], __m512d step_size) 
{
    __attribute__((aligned(64))) __m512d k1[4], k2[4], k3[4], k4[4],\
                                         temp_x[4], temp_v[4];
    double lambda = 0.0;
    int step = 0;

    while (lambda < lambda_max) {
        #pragma omp simd aligned(temp_v, christoffel: 64)
        CALCULATE_K_AVX512(k1, v)
        UPDATE_POSITIONS_AVX512(x, v, k1, step_size)
        #pragma omp simd aligned(temp_v, christoffel: 64)
        CALCULATE_K_AVX512(k2, temp_v)
        UPDATE_POSITIONS_AVX512(x, v, k2, step_size)
        #pragma omp simd aligned(temp_v, christoffel: 64)
        CALCULATE_K_AVX512(k3, temp_v)
        UPDATE_POSITIONS_AVX512(x, v, k3, step_size)
        #pragma omp simd aligned(temp_v, christoffel: 64)
        CALCULATE_K_AVX512(k4, temp_v)
        #pragma omp simd aligned(x, v, k1, k2, k3, k4: 64)
        for (int mu = 0; mu < 4; mu++) {
            __m512d k1_term = k1[mu];
            __m512d k2_term = _mm512_mul_pd(_mm512_set1_pd(2.0), k2[mu]);
            __m512d k3_term = _mm512_mul_pd(_mm512_set1_pd(2.0), k3[mu]);
            __m512d k4_term = k4[mu];
            __m512d sum_k = _mm512_add_pd(_mm512_add_pd(k1_term, k4_term), _mm512_add_pd(k2_term, k3_term));
            __m512d step_sum = _mm512_div_pd(_mm512_mul_pd(step_size, sum_k), _mm512_set1_pd(6.0));

            x[mu] = _mm512_add_pd(x[mu], step_sum);
            v[mu] = _mm512_add_pd(v[mu], step_sum);
        }

        lambda += _mm512_cvtsd_f64(step_size);
        store_geodesic_point_AVX(x, lambda);
    }
    printf("using AVX512 for geodesic\n");
}
#endif
