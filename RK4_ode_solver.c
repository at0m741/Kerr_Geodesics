#include "geodesics.h"

extern long double (*geodesic_points)[4];
extern int num_points;

void geodesic(long double x[4], long double v[4], long double lambda_max, long double christoffel[4][4][4], long double step_size, void (*store_point)(long double[], long double))
{
    long double k1[4], k2[4], k3[4], k4[4];
    long double temp_x[4], temp_v[4];
    ldouble_a32 lambda = 0;
    int step = 0;

    #pragma omp parallel shared(lambda, step, x, v) private(k1, k2, k3, k4, temp_x, temp_v)
    {
        while (lambda < lambda_max)
        {
            #pragma omp single
            {
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

                for (int mu = 0; mu < 4; mu++) {
                    k4[mu] = temp_v[mu];
                    for (int alpha = 0; alpha < 4; alpha++) {
                        for (int beta = 0; beta < 4; beta++) {
                            k4[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                        }
                    }
                }

                #pragma omp parallel for simd
                for (int mu = 0; mu < 4; mu++) {
                    x[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
                    v[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
                }
            }
			if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3]) || isnan(v[0]) || isnan(v[1]) || isnan(v[2]) || isnan(v[3])
				|| isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3]) || isinf(v[0]) || isinf(v[1]) || isinf(v[2]) || isinf(v[3]))
			{
				printf("Error: NaN or Inf detected\n");
				exit(1);
			}
            #pragma omp single
            {
                store_geodesic_point(x, lambda);
                if (step % 1000 == 0)
                    printf("Lambda %Lf: x = %Lf, y = %Lf, z = %Lf\n", lambda, x[0], x[1], x[2]);
                step++;
                lambda += step_size;
            }
        }
    }
}