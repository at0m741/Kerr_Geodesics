
#include "../headers/geodesics.h"


#ifdef __AVX2__
#include <immintrin.h>
void geodesic_AVX(__m256d x[4], __m256d v[4], float lambda_max, __m256d christoffel[4][4][4], __m256d step_size) {
    __asm__ __volatile__("");  
    printf("use AVX2 for geodesic\n");
    __m256d k1[4], k2[4], k3[4], k4[4], temp_x[4], temp_v[4];
    float lambda = 0.0;
    int step = 0;

    while (lambda < lambda_max) {
        #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++) {
            k1[mu] = v[mu];
            for (int nu = 0; nu < 4; nu++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d product = _mm256_mul_pd(christoffel[mu][nu][beta], 
                                        _mm256_mul_pd(v[nu], v[beta]));
                    k1[mu] = _mm256_sub_pd(k1[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k1[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k1[mu]));
        }
        #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++) {
            k2[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], 
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k2[mu] = _mm256_sub_pd(k2[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k2[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k2[mu]));
        }
        #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++) {
            k3[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], 
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k3[mu] = _mm256_sub_pd(k3[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k3[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k3[mu]));
        }
        #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++) {
            k4[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], 
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k4[mu] = _mm256_sub_pd(k4[mu], product);
                }
            }
        }
        #pragma omp simd aligned(x, v, k1, k2, k3, k4: ALIGNMENT)
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
        __m256d lambda_cast = _mm256_set1_pd(lambda);
        store_geodesic_point_AVX(x, lambda_cast);
    }
    __asm__ __volatile__("");
}

#elif defined(__GNUC__) || defined(__clang__)

    extern double (*geodesic_points)[4];
    extern int num_points;

    void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4], double step_size, void (*store_point)(double[], double))
    {
        double k1[4], k2[4], k3[4], k4[4];
        double temp_x[4], temp_v[4];
        ldouble_a lambda = 0;
        int step = 0;
        #ifdef _OPENMP
            printf("OpenMP is supported and used\n");
        #else
            printf("OpenMP is not supported\n");
        #endif
        #pragma omp parallel shared(lambda, step, x, v) private(k1, k2, k3, k4, temp_x, temp_v)
        {
            for(; lambda < lambda_max;)
            {
                #pragma omp single
                {
                    #pragma omp simd aligned (v, christoffel: ALIGNMENT)
                    for (int mu = 0; mu < 4; mu++) {
                        k1[mu] = v[mu];
                        for (int alpha = 0; alpha < 4; alpha++) {
                            for (int beta = 0; beta < 4; beta++) {
                                k1[mu] -= christoffel[mu][alpha][beta] * v[alpha] * v[beta];
                            }
                        }
                        temp_x[mu] = x[mu] + 0.5 * step_size * k1[mu];
                        temp_v[mu] = v[mu] + 0.5 * step_size * k1[mu];
                    }
                    #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
                    for (int mu = 0; mu < 4; mu++) {
                        k2[mu] = temp_v[mu];
                        for (int alpha = 0; alpha < 4; alpha++) {
                            for (int beta = 0; beta < 4; beta++) {
                                k2[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                            }
                        }
                        temp_x[mu] = x[mu] + 0.5 * step_size * k2[mu];
                        temp_v[mu] = v[mu] + 0.5 * step_size * k2[mu];
                    }
                    #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
                    for (int mu = 0; mu < 4; mu++) {
                        k3[mu] = temp_v[mu];
                        for (int alpha = 0; alpha < 4; alpha++) {
                            for (int beta = 0; beta < 4; beta++) {
                                k3[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                            }
                        }
                        temp_x[mu] = x[mu] + step_size * k3[mu];
                        temp_v[mu] = v[mu] + step_size * k3[mu];
                    }
                    #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
                    for (int mu = 0; mu < 4; mu++) {
                        k4[mu] = temp_v[mu];
                        for (int alpha = 0; alpha < 4; alpha++) {
                            for (int beta = 0; beta < 4; beta++) {
                                k4[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                            }
                        }
                    }

                    #pragma omp simd aligned (x, v, k1, k2, k3, k4: ALIGNMENT)
                    for (int mu = 0; mu < 4; mu++) {
                        x[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
                        v[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
                    }
                }
                #pragma omp single
                {
                    store_geodesic_point(x, lambda);
                    if (step % 1000 == 0)
                        fprintf(stderr, "Lambda %f: x = %f, y = %f, z = %f, t = %f\n", lambda, x[1], x[2], x[3], x[0]);
                    step++;
                    lambda += step_size;
                }
            }
        }
    }
#else
    #error "Unsupported compiler or configuration"
#endif
