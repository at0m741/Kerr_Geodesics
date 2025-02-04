#include <immintrin.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Geodesics.h"
#include <chrono>
double (*geodesic_points)[5] = NULL;
int num_points = 0;
#define NDIM 4




void calculate_Gamma_at_offset(double X[NDIM], int direction, double offset, double delta,
    double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], double Gamma_slice[NDIM][NDIM][NDIM]) {
    double X_offset[NDIM];
    memcpy(X_offset, X, sizeof(double) * NDIM);
    X_offset[direction] += offset;
    
    double tempGamma[NDIM][NDIM][NDIM];
    calculate_metric(X_offset, gcov, gcon);
    calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon);
    
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                Gamma_slice[rho][mu][nu] = tempGamma[rho][mu][nu];
            }
        }
    }
}

void calculate_ricci(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM]) {
    memset(Ricci, 0, sizeof(double) * NDIM * NDIM);
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int lambda = 0; lambda < NDIM; lambda++) {
                Ricci[mu][nu] += Riemann[lambda][mu][lambda][nu];
            }
        }
    }
	printf("Ricci tensor calculation completed\n");
	for (int mu = 0; mu < NDIM; mu++) {
		for (int nu = 0; nu < NDIM; nu++) {
			printf("Ricci[%d][%d] = %f\n", mu, nu, Ricci[mu][nu]);
		}
	}
}


int main(int argc, char **argv)
{
	double r0 = 20.0;
    double X[NDIM] = {0.0, r0, M_PI/2.0, 0.2};
    
	if (argc < 2) {
		printf("Usage: <options>\n");
		printf("       -G  -  geodesic calculation\n");
		printf("       -R  -  Riemann tensor calculation\n");
		printf("       -M  -  Metric tensor calculation (a = 0 -> Schwarzschild or a > 0 -> Kerr)\n");
		return 0;
	}


   if (strcmp(argv[1], "-R") == 0) {
        double gcovR[NDIM][NDIM], gconR[NDIM][NDIM];
        double christoffelR[NDIM][NDIM][NDIM];
        double Gamma_plus_h[NDIM][NDIM][NDIM][NDIM], Gamma_minus_h[NDIM][NDIM][NDIM][NDIM];
        double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM], Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM];
        double Riemann[NDIM][NDIM][NDIM][NDIM];
        double gcov_half[NDIM][NDIM], gcon_half[NDIM][NDIM];
		double Ricci[NDIM][NDIM];
        calculate_metric(X, gcovR, gconR);
        calculate_christoffel(X, DELTA, christoffelR, gcovR, gconR);

        for (int d = 0; d < NDIM; d++) {
            calculate_Gamma_at_offset(X, d, DELTA, DELTA, gcovR, gconR, Gamma_plus_h[d]);
            calculate_Gamma_at_offset(X, d, -DELTA, DELTA, gcovR, gconR, Gamma_minus_h[d]);
        }

        for (int d = 0; d < NDIM; d++) {
            calculate_Gamma_at_offset(X, d, DELTA / 2.0, DELTA, gcov_half, gcon_half, Gamma_plus_half_h[d]);
            calculate_Gamma_at_offset(X, d, -DELTA / 2.0, DELTA, gcov_half, gcon_half, Gamma_minus_half_h[d]);
        }

        calculate_riemann(christoffelR, Gamma_plus_h, Gamma_minus_h, Gamma_plus_half_h, Gamma_minus_half_h, Riemann, DELTA);
        printf("Riemann tensor calculation completed\n");
        print_riemann(Riemann);
        check_riemann_symmetries(Riemann, 1e-6);
		contract_riemann(Riemann, Ricci, gconR);
		return 0;
    }else if (strcmp(argv[1], "-G") == 0) {
		double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
		double g_tt = gcov[0][0];
		double g_tphi = gcov[0][3];
		double g_phiphi = gcov[3][3];
		double omega = 1.0 / (pow(r0, 1.5) + a);
		double denom = -(g_tt + 2.0 * g_tphi * omega + g_phiphi * omega * omega);
		double vt = 1.0 / sqrt(denom);
		double v[NDIM] = {vt, 0.0, 0.0, omega * vt};
		double christoffel[NDIM][NDIM][NDIM];
		double dt = 0.00910;
		__m256d X_avx[NDIM], v_avx[NDIM];
		for (int i = 0; i < NDIM; i++) {
			X_avx[i] = _mm256_set1_pd(X[i]);
			v_avx[i] = _mm256_set1_pd(v[i]);
		}
		__m256d christoffel_avx[NDIM][NDIM][NDIM];
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				for (int k = 0; k < NDIM; k++) {
					christoffel_avx[i][j][k] = _mm256_set1_pd(christoffel[i][j][k]);
				}
			}
		}


	    calculate_metric(X, gcov, gcon);
		calculate_christoffel(X, DELTA, christoffel, gcov, gcon);
		auto start = std::chrono::high_resolution_clock::now();
		geodesic_AVX(X_avx, v_avx, max_dt + 4, ( __m256d (*)[NDIM][NDIM] )christoffel_avx, _mm256_set1_pd(dt));
		auto end = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;
		printf("Elapsed time: %f\n", elapsed_seconds.count());
		write_vtk_file("geodesic.vtk");
		if (geodesic_points != NULL) {
			free(geodesic_points);
		}
		return 0;
	} else if (strcmp(argv[1], "-M") == 0) {
		double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
		calculate_metric(X, gcov, gcon);
	} else {
		printf("Invalid option\n");
	}
}
