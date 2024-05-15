#include "../headers/geodesics.h"

extern double (*geodesic_points)[4];
extern int num_points;

void geodesic(double x[4], double v[4], double lambda_max, double christoffel[4][4][4], \
    double step_size, void (*store_point)(double[], double))
{
    double k1[4], k2[4], k3[4], k4[4], k5[4];
    double temp_x[4], temp_v[4];
    double lambda = 0;
    short step = 0;
    (void)store_point;
    #ifdef _OPENMP
        printf("OpenMP is supported and used\n");
    #else
        printf("OpenMP is not supported\n");
    #endif
    
    #pragma omp parallel shared(lambda, step, x, v) private(k1, k2, k3, k4, temp_x, temp_v)
    {
        for(; lambda < lambda_max;)
        {
            #pragma omp simd aligned(v, christoffel: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                k1[mu] = v[mu];
                for (int alpha = 0; alpha < 4; alpha++) {
                    for (int beta = 0; beta < 4; beta++) {
                        double temp = christoffel[mu][alpha][beta] * v[alpha] * v[beta];
                        k1[mu] -= temp;
                    }
                }
                temp_x[mu] = x[mu] + 0.5 * step_size * k1[mu];
                temp_v[mu] = v[mu] + 0.5 * step_size * k1[mu];
            }
            
            #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                k2[mu] = temp_v[mu];
                for (int alpha = 0; alpha < 4; alpha++) {
                    for (int beta = 0; beta < 4; beta++) {
                        double temp = christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                        k2[mu] -= temp;
                    }
                }
                temp_x[mu] = x[mu] + 0.5 * step_size * k2[mu];
                temp_v[mu] = v[mu] + 0.5 * step_size * k2[mu];
            }
            
            #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                k3[mu] = temp_v[mu];
                for (int alpha = 0; alpha < 4; alpha++) {
                    for (int beta = 0; beta < 4; beta++) {
                        double temp = christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                        k3[mu] -= temp;
                    }
                }
                temp_x[mu] = x[mu] + step_size * k3[mu];
                temp_v[mu] = v[mu] + step_size * k3[mu];
            }
            
            #pragma omp simd aligned(temp_v, christoffel: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                k4[mu] = temp_v[mu];
                for (int alpha = 0; alpha < 4; alpha++) {
                    for (int beta = 0; beta < 4; beta++) {
                        double temp = christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                        k4[mu] -= temp;
                    }
                }
            }
            
            #pragma omp simd aligned(x, v, k1, k2, k3, k4: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                x[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
                v[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
            }
            
            #pragma omp simd aligned(v, christoffel: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                k5[mu] = v[mu];
                for (int alpha = 0; alpha < 4; alpha++) {
                    for (int beta = 0; beta < 4; beta++) {
                        double temp = christoffel[mu][alpha][beta] * v[alpha] * v[beta];
                        k5[mu] -= temp;
                    }
                }
            }

            #pragma omp simd aligned(x, v, k1, k2, k3, k4, k5: ALIGNMENT)
            for (int mu = 0; mu < 4; mu++) {
                x[mu] += step_size * (k1[mu] + 4 * k2[mu] + k3[mu]) / 6.0;
                v[mu] += step_size * (k1[mu] + 4 * k2[mu] + k3[mu]) / 6.0;
            }

            #pragma omp single
            {
                store_geodesic_point(x, lambda);
                if (step % 1000 == 0)
                    printf("Cycle = %f, r = %f, θ = %f, φ = %f, t = %f\n", lambda, x[1], x[2], x[3], x[0]);
                step++;
                lambda += step_size;
            }
        }
    }
}
