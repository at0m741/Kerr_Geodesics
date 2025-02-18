#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

int Riemann_tensor(const char *metric) {
	double r0 = 20.0;
	double X[NDIM] = {0.0, r0, M_PI/4.0, 0.0};
	double gcovR[NDIM][NDIM], gconR[NDIM][NDIM];
	double christoffelR[NDIM][NDIM][NDIM];
	double Gamma_plus_h[NDIM][NDIM][NDIM][NDIM], Gamma_minus_h[NDIM][NDIM][NDIM][NDIM];
	double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM], Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM];
	double Riemann[NDIM][NDIM][NDIM][NDIM];
	double gcov_half[NDIM][NDIM], gcon_half[NDIM][NDIM];
	double Ricci[NDIM][NDIM];
	if (strcmp(metric, "ds") == 0) {
		printf("KDS metric calculation\n");
		calculate_metric_kds(X, gcovR, gconR);
	} else {
		calculate_metric(X, gcovR, gconR);
	}
	calculate_christoffel(X, DELTA, christoffelR, gcovR, gconR, "kerr");

	for (int d = 0; d < NDIM; d++) {
		calculate_Gamma_at_offset(X, d, DELTA, DELTA, gcovR, gconR, Gamma_plus_h[d], "kerr");
		calculate_Gamma_at_offset(X, d, -DELTA, DELTA, gcovR, gconR, Gamma_minus_h[d], "kerr");
	}

	for (int d = 0; d < NDIM; d++) {
		calculate_Gamma_at_offset(X, d, DELTA / 2.0, DELTA, gcov_half, \
									gcon_half, Gamma_plus_half_h[d], "kerr");
		calculate_Gamma_at_offset(X, d, -DELTA / 2.0, DELTA, gcov_half, \
									gcon_half, Gamma_minus_half_h[d], "kerr");
	}

	calculate_riemann(christoffelR, Gamma_plus_h, Gamma_minus_h, \
					  Gamma_plus_half_h, Gamma_minus_half_h, Riemann, DELTA);
	printf("Riemann tensor calculation completed\n");
	print_riemann(Riemann);
	contract_riemann(Riemann, Ricci, gconR);
	return 0;
}

int Geodesics_prob() {
	double r0 = 20.0;
	double X[NDIM] = {0.4, r0, M_PI/3.8, 0.2};
	double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
	calculate_metric(X, gcov, gcon);
	double g_tt = gcov[0][0];
	double g_tphi = gcov[0][3];
	double g_phiphi = gcov[3][3];
	double Omega = 1.0 / (pow(r0, 1.5) + a);
	double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
	double vt = 1.0 / sqrt(denom);
	double v[NDIM] = {vt, 0.0, 0.0, Omega * vt};
	double norm = g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
	printf("Norme quadrivecteur = %e (doit être proche de 0)\n", norm);
	double dt = 0.00910;
	double christoffel[NDIM][NDIM][NDIM];
	calculate_christoffel(X, DELTA, christoffel, gcov, gcon, "kerr");
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
	compute_photon_properties(gcov, v);
	auto start = std::chrono::high_resolution_clock::now();
	geodesic_AVX(X_avx, v_avx, max_dt + 4, ( __m256d (*)[NDIM][NDIM] )christoffel_avx, _mm256_set1_pd(dt));
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	printf("Elapsed time: %f\n", elapsed_seconds.count());
	write_vtk_file("output/geodesic.vtk");
	if (geodesic_points != NULL) {
		free(geodesic_points);
	}
	return 0;
}

int light_geodesics_prob() {
	double r0 = 20.0;
	double X[NDIM] = {0.4, r0, M_PI/3.0, 0.2};
	double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
	calculate_metric(X, gcov, gcon);
	double g_tt = gcov[0][0];
	double g_tphi = gcov[0][3];
	double g_phiphi = gcov[3][3];
	double Omega = 1.0 / (pow(r0, 1.5) + a);
	double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
	double vt = 1.0 / sqrt(fabs(denom));
	double v[NDIM] = {vt, 0.0, 0.0, 3.5 * Omega * vt};
	double norm = g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
	printf("Norme quadrivecteur = %e (doit être proche de 0)\n", norm);
	double dt = 0.00910;
	double christoffel[NDIM][NDIM][NDIM];
	calculate_christoffel(X, DELTA, christoffel, gcov, gcon, "kerr");
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
	generate_blackhole_image();
	if (geodesic_points != NULL) {
		free(geodesic_points);
	}
	return 0;
}

int Metric_prob() {
	double r0 = 20.0;
	double X[NDIM] = {0.4, r0, M_PI/4.0, 0.2};
	double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
	printf("=====================================================\n");
	calculate_metric(X, gcov, gcon);

	printf("=====================================================\n");
	double gcov_KN[NDIM][NDIM], gcon_KN[NDIM][NDIM];
	calculate_metric_kerr_newman(X, gcov_KN, gcon_KN);

	printf("=====================================================\n");
	double gcov_kds[NDIM][NDIM], gcon_kds[NDIM][NDIM];
	calculate_metric_kds(X, gcov_kds, gcon_kds);
	return 0;
}

