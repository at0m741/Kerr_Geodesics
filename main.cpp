#include <immintrin.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Geodesics.h"
#include <chrono>
double (*geodesic_points)[5] = NULL;
int num_points = 0;
#define DELTA 1e-6
#define NDIM 4

int main(int argc, char **argv)
{
    double r0 = 20.0;
    double X[NDIM] = {0.4, r0, M_PI/4.0, 0.2};
    double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
    calculate_metric(X, gcov, gcon);
    double g_tt = gcov[0][0];
    double g_tphi = gcov[0][3];
    double g_phiphi = gcov[3][3];
    double Omega = 1.0 / (pow(r0, 1.5) + a);
    double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
    double vt = 1.0 / sqrt(denom);
    double v[NDIM] = {vt, 0.0, 0.0, Omega * vt};
    double dt = 0.00910;
    double christoffel[NDIM][NDIM][NDIM];
    calculate_christoffel(X, DELTA, christoffel, gcov, gcon);
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
}
