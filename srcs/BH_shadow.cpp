#include <Geodesics.h>

#define WIDTH 64	
#define HEIGHT 64
unsigned char image[HEIGHT][WIDTH];
extern int num_points;
extern double a;
extern double (*geodesic_points)[5];

double calculate_impact_parameter(double p_t, double p_phi, double g_tt, double g_tphi, double g_phiphi) {
    return - (g_tphi * p_t + g_phiphi * p_phi) / (g_tt * p_t + g_tphi * p_phi);
}

double calculate_emission_angle(double p_r, double p_phi, double g_rr, double g_phiphi) {
    return atan2(sqrt(g_phiphi) * p_phi, sqrt(g_rr) * p_r) * 180.0 / M_PI; 
}

double b_critique_kerr(double a, int sense) {
    double r_ph = 2.0 * M;
    for (int i = 0; i < 10; i++) {  
        double f = 2 * r_ph - 3 * M + 4 * a * sqrt(M / r_ph) * sense;
        double df = 2 - 2 * a * sqrt(M / (r_ph * r_ph * r_ph)) * sense;
        r_ph -= f / df;
    }
    return (r_ph * r_ph + a * a) / (a + sense * sqrt(r_ph * r_ph * r_ph / M));
}

int compute_photon_properties(double g[4][4], double p[4]) {
    double b = calculate_impact_parameter(p[0], p[3], g[0][0], g[0][3], g[3][3]);
    double alpha = calculate_emission_angle(p[1], p[3], g[1][1], g[3][3]);

    printf("Impact parameter b = %f\n", b);
    printf("Emission angle alpha = %f\n", alpha);
	if (b < b_critique_kerr(a, 1) && b > b_critique_kerr(a, -1)) {
		printf("Photon is captured by the black hole\n");
	} else {
		printf("Photon is not captured by the black hole\n");
	}
	return 0;
}

/*  */
/* int init_photon_global( */
/*     double alpha,  */
/*     double beta, */
/*     double e_t[4],  */
/*     double e_x[4],  */
/*     double e_y[4],  */
/*     double e_z[4], */
/*     double gcov[4][4], */
/*     double p[4], */
/*     double X[4]) { */
/*      */
/*     double p_local[4]; */
/*     p_local[0] = -1.0;  */
/*     p_local[1] = 0.0; */
/*  */
/*     p_local[2] = alpha * sqrt(fabs(gcov[2][2])) / fabs(1.0 + fabs(gcov[0][0]) + fabs(gcov[3][3])); */
/*     p_local[3] = -beta * sqrt(fabs(gcov[3][3])) / fabs(1.0 + fabs(gcov[0][0]) + fabs(gcov[3][3]));   */
/*  */
/*     printf("Debug: p_theta = %f, p_phi = %f\n", p_local[2], p_local[3]); */
/* 	printf("Debug: p_t = %f, p_r = %f, p_theta = %f, p_phi = %f\n", p_local[0], p_local[1], p_local[2], p_local[3]); */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         p[mu] = p_local[0] * e_t[mu] */
/*               + p_local[1] * e_x[mu] */
/*               + p_local[2] * e_y[mu] */
/*               + p_local[3] * e_z[mu]; */
/*     } */
/*  */
/*     double E  = fabs(-gcov[0][0] * p[0] - gcov[0][3] * p[3]); */
/*     double Lz = gcov[0][3] * p[0] + gcov[3][3] * p[3]; */
/*  */
/*     p_local[2] = alpha * sqrt(fabs(gcov[2][2])) / fabs(E);   */
/*     p_local[3] = beta / fabs(E);  */
/*  */
/*     if (fabs(E) > 2.0) { */
/*         E = copysign(2.0, E); */
/*     } */
/*     if (fabs(Lz) > 10.0) { */
/*         Lz = copysign(10.0, Lz); */
/*     } */
/*  */
/*     double p_theta = p[2];  */
/*     double theta = X[2]; */
/*     double cosTheta2 = cos(theta) * cos(theta); */
/*     if (cosTheta2 < 1e-10) cosTheta2 = 1e-10;   */
/*  */
/*     double Q = fabs(gcov[2][2] * p_theta * p_theta) + */
/*         (a * a * (E * E) - (Lz * Lz)) * cosTheta2; */
/*  */
/*     if (Q < 0) { */
/* 		image[0][0] = 255; */
/*         Q = fabs(gcov[2][2]) * p_theta * p_theta +  */
/*             (a * a * (E * E) - (Lz * Lz)) * cosTheta2; */
/*     } */
/*  */
/*     printf("Q = %f\n", Q); */
/*  */
/*     double sum = 0.0; */
/*     for (int i = 1; i < 4; i++) {  */
/*         for (int j = 1; j < 4; j++) { */
/*             sum += gcov[i][j] * p[i] * p[j]; */
/*         } */
/*     } */
/*  */
/*     printf("sum = %f\n", sum); */
/*  */
/*     p[0] = (-gcov[0][3] * p[3] - sqrt(fabs(-gcov[0][0] * sum))) / gcov[0][0]; */
/*  */
/*     printf("p[0] = %f\n", p[0]); */
/*     printf("Corrected Q = %f\n", Q); */
/* 	return 0; */
/* } */
/*  */
/*  */
/*  */
/* void build_kerr_tetrad(double X[4], double gcov[4][4], */
/*                        double e_t[4], double e_r[4], */
/*                        double e_theta[4], double e_phi[4])  */
/* { */
/*     double r     = X[1]; */
/*     double theta = X[2]; */
/*  */
/*     double Delta = r*r - 2*M*r + a*a; */
/*     double Sigma = r*r + a*a*cos(theta)*cos(theta); */
/*     double A     = (r*r + a*a)*(r*r + a*a) - a*a * Delta * sin(theta)*sin(theta); */
/*  */
/*     e_t[0] = sqrt(A / (Sigma * Delta)); */
/*     e_t[1] = 0.0; */
/*     e_t[2] = 0.0; */
/*     e_t[3] = (2 * a * M * r) / sqrt(Sigma * A * Delta); */
/*  */
/*     e_r[0] = 0.0; */
/*     e_r[1] = sqrt(Delta / Sigma); */
/*     e_r[2] = 0.0; */
/*     e_r[3] = 0.0; */
/*  */
/*     e_theta[0] = 0.0; */
/*     e_theta[1] = 0.0; */
/*     e_theta[2] = 1.0 / sqrt(Sigma); */
/*     e_theta[3] = 0.0; */
/*  */
/*     e_phi[0] = 0.0; */
/*     e_phi[1] = 0.0; */
/*     e_phi[2] = 0.0; */
/*     if (fabs(sin(theta)) > 1e-10) { */
/*         e_phi[3] = sqrt(Sigma / A) / sin(theta); */
/*     } else { */
/*         e_phi[3] = 0.0;  */
/*     } */
/*  */
/*     double dot_et_et = 0, dot_er_er = 0, dot_et_er = 0, dot_er_et = 0; */
/*     for (int mu = 0; mu < 4; mu++) { */
/*         for (int nu = 0; nu < 4; nu++) { */
/*             dot_et_et += gcov[mu][nu] * e_t[mu] * e_t[nu]; */
/*             dot_er_er += gcov[mu][nu] * e_r[mu] * e_r[nu]; */
/*             dot_et_er += gcov[mu][nu] * e_t[mu] * e_r[nu]; */
/*             dot_er_et += gcov[mu][nu] * e_r[mu] * e_t[nu]; */
/*         } */
/*     } */
/*     printf("Test normalisation e_t.e_t = %f (devrait être -1)\n", dot_et_et); */
/*     printf("Test normalisation e_r.e_r = %f (devrait être 1)\n", dot_er_er); */
/*     printf("Test normalisation e_t.e_r = %f (devrait être 0)\n", dot_et_er); */
/*     printf("Test normalisation e_r.e_t = %f (devrait être 0)\n", dot_er_et); */
/* } */
/*  */
/* int solve_geodesic_AVX(double X[NDIM], __m256d p[NDIM], double *Q_out) { */
/* 	Connexion connexion; */
/* 	Metric metric_obj; */
/*     double lambda_max = 6.0;   */
/*     __m256d step_size = _mm256_set1_pd(0.00019); */
/*     double r_horizon = 1.0 + sqrt(1.0 - a * a);  */
/*     double christoffel[NDIM][NDIM][NDIM]; */
/*     double gcov[NDIM][NDIM], gcon[NDIM][NDIM]; */
/*  */
/*     metric_obj.calculate_metric(X, gcov, gcon); */
/*     connexion.calculate_christoffel(X, DELTA, christoffel, gcov, gcon, "kerr"); */
/*  */
/*     __m256d X_avx[NDIM], p_avx[NDIM], christoffel_avx[NDIM][NDIM][NDIM]; */
/*      */
/*     for (int i = 0; i < NDIM; i++) { */
/*         X_avx[i] = _mm256_set1_pd(X[i]); */
/*         p_avx[i] = p[i]; */
/*         for (int j = 0; j < NDIM; j++) */
/*             for (int k = 0; k < NDIM; k++) */
/*                 christoffel_avx[i][j][k] = _mm256_set1_pd(christoffel[i][j][k]); */
/*     } */
/*  */
/*     num_points = 0;   */
/*  */
/*     geodesic_AVX(X_avx, p_avx, lambda_max, christoffel_avx, step_size); */
/*      */
/*     double E = -gcov[0][0] * p[0][0] - gcov[0][3] * p[3][0]; */
/*     double Lz = gcov[0][3] * p[0][0] + gcov[3][3] * p[3][0]; */
/*  */
/*     double r_ps = 2 + 2 * cos( (2.0 / 3.0) * acos(-a) ); */
/*     double Lz_crit = - (r_ps * r_ps + a * a - 2 * r_ps) / a; */
/*     double Q = ((Lz - a * E) * (Lz - a * E)) / (a * a) - r_ps * r_ps; */
/* 	*Q_out = Q;	 */
/*     for (int i = 0; i < num_points; i++) { */
/*         double x = geodesic_points[i][0]; */
/*         double y = geodesic_points[i][1]; */
/*         double z = geodesic_points[i][2]; */
/*  */
/*         double r = sqrt(x*x + y*y + z*z); */
/*  */
/*  */
/* 		if (r < r_horizon || Q < -1e-10  ) { */
/* 			printf("Photon absorbed at r = %f, E = %f, Lz = %f, Q = %f, r_horizon = %f, r_ps = %f\n", r, E, Lz, Q, r_horizon, r_ps); */
/* 			return 1; */
/* 		} */
/*  */
/*     } */
/*  */
/*     printf("Photon escaped, r = %f, E = %f, Lz = %f, Q = %f\n", sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3]), E, Lz, Q);  */
/*     return 0;  */
/* } */
/*  */
/*  */
/* void generate_blackhole_shadow() { */
/* 	Metric metric_obj; */
/*     double r_obs = 5 * M; */
/*     double theta_obs = M_PI / 2.0;  */
/*     double phi_obs   = 0.0; */
/*  */
/*     double X[4] = {1.0, r_obs, theta_obs, phi_obs}; */
/*  */
/*     double gcov[4][4], gcon[4][4]; */
/*     metric_obj.calculate_metric(X, gcov, gcon); */
/*  */
/*      */
/* 	double e_t_loc[4], e_r_loc[4], e_theta_loc[4], e_phi_loc[4]; */
/* 	build_kerr_tetrad(X, gcov, e_t_loc, e_r_loc, e_theta_loc, e_phi_loc); */
/*  */
/* 	double e_t[4], e_x[4], e_y[4], e_z[4]; */
/* 	for (int mu = 0; mu < 4; mu++) { */
/* 		e_t[mu] = e_t_loc[mu]; */
/* 		e_x[mu] = -e_r_loc[mu]; */
/* 		e_y[mu] = e_phi_loc[mu]; */
/* 		e_z[mu] = e_theta_loc[mu]; */
/* 	} */
/*  */
/*  */
/* 	double fov = 4.4; */
/* 	double aspect_ratio = (double)WIDTH / (double)HEIGHT; */
/* 	for (int i = 0; i < HEIGHT; i++) { */
/* 		for (int j = 0; j < WIDTH; j++) { */
/* 			double x_shift = 0.0; */
/* 			double alpha = fov * (2.0 * ((double)j / (WIDTH - 1)) - 1.0); */
/* 			double beta  = -fov * (2.0 * ((double)i / (HEIGHT - 1)) - 1.0); */
/* 			printf("alpha = %f, beta = %f\n", alpha, beta); */
/*  */
/* 			double p[4]; */
/* 			init_photon_global(alpha, beta, e_t, e_x, e_y, e_z, gcov, p, X); */
/*  */
/* 			__m256d p_avx[4], X_avx[4]; */
/* 			for (int mu = 0; mu < 4; mu++) { */
/* 				p_avx[mu] = _mm256_set1_pd(p[mu]); */
/* 				X_avx[mu] = _mm256_set1_pd(X[mu]);  */
/* 			} */
/*  */
/* 			double r_horizon = 1.0 + sqrt(1.0 - a * a); */
/* 			double Q;  */
/* 			int total_hit = solve_geodesic_AVX(X, p_avx, &Q);  */
/* 			double r = sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3]); */
/* 			double intensity = exp(-fabs(Q) / (2 * a * a)) * (1.0 - sqrt(r_horizon / r)); */
/* 			if (total_hit >= 1) { */
/* 				printf("r = %f, Q = %f, intensity = %f\n", r, Q, intensity); */
/* 				image[i][j] = (unsigned char)(255.0 * intensity); */
/* 			} else { */
/* 				image[i][j] = 0;  */
/* 				printf("Intensity = %f\n", intensity); */
/*  */
/* 			} */
/* 		} */
/* 	} */
/*  */
/*     FILE *f = fopen("blackhole_shadow.ppm","wb"); */
/*     fprintf(f,"P5\n%d %d\n255\n", WIDTH, HEIGHT); */
/*     fwrite(image, 1, WIDTH*HEIGHT, f); */
/*     fclose(f); */
/*     printf("Image saved: blackhole_shadow.ppm\n"); */
/* } */
