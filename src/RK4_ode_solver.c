#include "../headers/geodesics.h"

void geodesic_AVX(__m256d x[4], __m256d v[4], double lambda_max, __m256d christoffel[4][4][4], __m256d step_size) 
{
    __attribute__((aligned(32))) __m256d k1[4], k2[4], k3[4], k4[4], temp_x[4], temp_v[4];
    double lambda = 0.0;
    int step = 0;

    while (lambda < lambda_max) {
        #pragma ivdep
        #pragma omp simd aligned(temp_v, christoffel: 32)
        CALCULATE_K_AVX2(k1, v)
        UPDATE_POSITIONS_AVX2(x, v, k1, step_size)
        #pragma ivdep
        #pragma omp simd aligned(temp_v, christoffel: 32)
        CALCULATE_K_AVX2(k2, temp_v)
        UPDATE_POSITIONS_AVX2(x, v, k2, step_size)
        #pragma ivdep
        #pragma omp simd aligned(temp_v, christoffel: 32)
        CALCULATE_K_AVX2(k3, temp_v)
        UPDATE_POSITIONS_AVX2(x, v, k3, step_size)
        #pragma ivdep
        #pragma omp simd aligned(temp_v, christoffel: 32)
        CALCULATE_K_AVX2(k4, temp_v)
        #pragma ivdep
        #pragma omp simd aligned(x, v, k1, k2, k3, k4: 32)
        for (int mu = 0; mu < 4; mu++) {
            __m256d k1_term = k1[mu];
            __m256d k2_term = _mm256_mul_pd(_mm256_set1_pd(2.0), k2[mu]);
            __m256d k3_term = _mm256_mul_pd(_mm256_set1_pd(2.0), k3[mu]);
            __m256d k4_term = k4[mu];
            __m256d sum_k = _mm256_add_pd(_mm256_add_pd(k1_term, k4_term), _mm256_add_pd(k2_term, k3_term));
            __m256d step_sum = _mm256_div_pd(_mm256_mul_pd(step_size, sum_k), _mm256_set1_pd(6.0));

            x[mu] = _mm256_add_pd(x[mu], step_sum);
            v[mu] = _mm256_add_pd(v[mu], step_sum);
        }

        lambda += _mm256_cvtsd_f64(step_size);
        store_geodesic_point_AVX(x, lambda);
    }
    printf("use AVX2 for geodesic\n");
}


