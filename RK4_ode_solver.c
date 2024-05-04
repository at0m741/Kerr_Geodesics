#include "geodesics.h"

void geodesic(long double x[4], long double v[4], long double lambda_max, long double christoffel[4][4][4], long double step_size, void (*store_point)(long double[], long double))
{
    long double k1[4], k2[4], k3[4], k4[4];
    long double temp_x[4], temp_v[4];
    ldouble_a32 lambda = 0;
	int step = 0;
    while (lambda < lambda_max) 
	{
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k1[mu] = v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k1[mu] -= christoffel[mu][alpha][beta] * v[alpha] * v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + 0.5 * step_size * k1[mu];
            temp_v[mu] = v[mu] + 0.5 * step_size * k1[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k2[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k2[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + 0.5 * step_size * k2[mu];
            temp_v[mu] = v[mu] + 0.5 * step_size * k2[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k3[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k3[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + step_size * k3[mu];
            temp_v[mu] = v[mu] + step_size * k3[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k4[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k4[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            x[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
            v[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
        }
		store_geodesic_point(x, lambda);
		if (step % 1000 == 0)
			printf("Lambda %Lf: x = %Lf, y = %Lf, z = %Lf\n", lambda, x[0], x[1], x[2]);
		step++;
        lambda += step_size;
	}
}