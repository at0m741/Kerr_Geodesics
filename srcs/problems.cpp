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

#define WIDTH 64
#define HEIGHT 64
unsigned char image[HEIGHT][WIDTH];



void init_photon_global(
    double alpha, 
    double beta,
    double e_t[4], 
    double e_x[4], 
    double e_y[4], 
    double e_z[4],
    double gcov[4][4],
    double p[4]
) {
    double px_loc = -9.0;
    double py_loc = alpha;
    double pz_loc = beta;
    double pt_loc = sqrt(px_loc*px_loc + py_loc*py_loc + pz_loc*pz_loc);

    double p_local[4];
    p_local[0] = pt_loc;
    p_local[1] = px_loc;
    p_local[2] = py_loc;
    p_local[3] = pz_loc;

    for (int mu = 0; mu < 4; mu++) {
        p[mu] = p_local[0]*e_t[mu]
              + p_local[1]*e_x[mu]
              + p_local[2]*e_y[mu]
              + p_local[3]*e_z[mu];
		printf("p[%d] = %f\n", mu, p[mu]);
    }

    double norm = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            norm += gcov[mu][nu] * p[mu] * p[nu];
			printf("gcov[%d][%d] = %f\n", mu, nu, gcov[mu][nu]);
        }
    }

    if (fabs(norm) > 1e-12) {
        double scale = 1.0 / sqrt(fabs(norm));
        for (int mu = 0; mu < 4; mu++) {
            p[mu] *= scale;
			printf("p[%d] = %f\n", mu, p[mu]);
        }
    }
}


void build_minkowski_tetrad(double e_t[4], double e_x[4], double e_y[4], double e_z[4]) {
    e_t[0] = 1.0; e_t[1] = 0.0; e_t[2] = 0.0; e_t[3] = 0.0;
    e_x[0] = 0.0; e_x[1] = 1.0; e_x[2] = 0.0; e_x[3] = 0.0;
    e_y[0] = 0.0; e_y[1] = 0.0; e_y[2] = 1.0; e_y[3] = 0.0;
    e_z[0] = 0.0; e_z[1] = 0.0; e_z[2] = 0.0; e_z[3] = 1.0;
	for (int i = 0; i < 4; i++) {
		printf("e_t[%d] = %f\n", i, e_t[i]);
		printf("e_x[%d] = %f\n", i, e_x[i]);
		printf("e_y[%d] = %f\n", i, e_y[i]);
		printf("e_z[%d] = %f\n", i, e_z[i]);
	}
}

int solve_geodesic_AVX(double X[NDIM], __m256d p[NDIM]) {
    double lambda_max = 100.0;  
    __m256d step_size = _mm256_set1_pd(0.00910);
    double r_horizon = 1.0 + sqrt(1.0 - a * a); 

    double christoffel[NDIM][NDIM][NDIM];
    double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

    calculate_metric(X, gcov, gcon);
    calculate_christoffel(X, DELTA, christoffel, gcov, gcon, "kerr");

    __m256d X_avx[NDIM], p_avx[NDIM], christoffel_avx[NDIM][NDIM][NDIM];
    
    for (int i = 0; i < NDIM; i++) {
        X_avx[i] = _mm256_set1_pd(X[i]);
        p_avx[i] = p[i];
        for (int j = 0; j < NDIM; j++)
            for (int k = 0; k < NDIM; k++)
                christoffel_avx[i][j][k] = _mm256_set1_pd(christoffel[i][j][k]);
    }

    num_points = 0;  

    geodesic_AVX(X_avx, p_avx, lambda_max, christoffel_avx, step_size);
    for (int i = 0; i < num_points; i++) {
        double x = geodesic_points[i][0];
		double y = geodesic_points[i][1];
		double z = geodesic_points[i][2];

		double r = sqrt(x*x + y*y + z*z);
		/* printf("r = %f\n", r); */
		if (r < r_horizon) {
			printf("Photon absorbed at r = %f\n", r);
			return 1;
		}    
	}

	printf("Photon escaped, r = %f\n", sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3]));
    return 0; 
}

void save_image(const char *filename, unsigned char image[HEIGHT][WIDTH]) {
    FILE *f = fopen(filename, "wb");
    fprintf(f, "P5\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(image, 1, WIDTH * HEIGHT, f);
    fclose(f);
}

void generate_blackhole_shadow() {
    double r_obs = 10.0;
    double theta_obs = M_PI / 2.0;
    double phi_obs   = 0.0;

    double X[4] = {0.0, r_obs, theta_obs, phi_obs};

    double gcov[4][4], gcon[4][4];
    calculate_metric(X, gcov, gcon);

    double e_t[4], e_x[4], e_y[4], e_z[4];
    build_minkowski_tetrad(e_t, e_x, e_y, e_z);

    double fov = 1.0;
    for(int i = 0; i < HEIGHT; i++) {
        for(int j = 0; j < WIDTH; j++) {
            double alpha = 0.0;
            double beta  = 0.0; 

            double p[4];
            init_photon_global(alpha, beta, e_t, e_x, e_y, e_z, gcov, p);

            __m256d p_avx[4], X_avx[4];
            for(int mu=0; mu<4; mu++) {
                p_avx[mu] = _mm256_set1_pd(p[mu]);
                X_avx[mu] = _mm256_set1_pd(X[mu]); 
            }

            int hit_horizon = solve_geodesic_AVX(X, p_avx);
            image[i][j] = (hit_horizon ? 0 : 255);
        }
    }

    FILE *f = fopen("blackhole_shadow.ppm","wb");
    fprintf(f,"P5\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(image, 1, WIDTH*HEIGHT, f);
    fclose(f);
    printf("Image saved: blackhole_shadow.ppm\n");
}
