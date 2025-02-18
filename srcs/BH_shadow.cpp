#include <Geodesics.h>

#define WIDTH 128	
#define HEIGHT 128
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

void compute_photon_properties(double g[4][4], double p[4]) {
    double b = calculate_impact_parameter(p[0], p[3], g[0][0], g[0][3], g[3][3]);
    double alpha = calculate_emission_angle(p[1], p[3], g[1][1], g[3][3]);

    printf("Impact parameter b = %f\n", b);
    printf("Emission angle alpha = %f\n", alpha);
	if (b < b_critique_kerr(a, 1) && b > b_critique_kerr(a, -1)) {
		printf("Photon is captured by the black hole\n");
	} else {
		printf("Photon is not captured by the black hole\n");
	}
}


void init_photon_global(
    double alpha, 
    double beta,
    double e_t[4], 
    double e_x[4], 
    double e_y[4], 
    double e_z[4],
    double gcov[4][4],
    double p[4]) {
    double p_local[4];
    p_local[0] = 1.0; 
    p_local[1] = 0.0;
    p_local[2] = alpha;
    p_local[3] = beta;

    for (int mu = 0; mu < 4; mu++) {
        p[mu] = p_local[0] * e_t[mu]
              + p_local[1] * e_x[mu]
              + p_local[2] * e_y[mu]
              + p_local[3] * e_z[mu];
    }




    double E = -gcov[0][0] * p[0] - gcov[0][3] * p[3];
    double Lz = gcov[0][3] * p[0] + gcov[3][3] * p[3];

    double r_ps = 2 + 2 * cos((2.0 / 3.0) * acos(-a));
    double Lz_crit = - (r_ps * r_ps + a * a - 2 * r_ps) / a;
    double Q = ((Lz - a * E) * (Lz - a * E)) / (a * a) - r_ps * r_ps;

	double sum = 0.0;
	for (int i = 1; i < 4; i++) { 
		for (int j = 1; j < 4; j++) {
			sum += gcov[i][j] * p[i] * p[j];
		}
	}
	printf("sum = %f\n", sum);
	p[0] = ( - gcov[0][3]*p[3] - sqrt( fabs(-gcov[0][0]*sum) ) ) / gcov[0][0];
	printf("p[0] = %f\n", p[0]);

}


void build_kerr_tetrad(double X[4], double gcov[4][4],
                       double e_t[4], double e_r[4],
                       double e_theta[4], double e_phi[4]) 
{
    double g_tt         = gcov[0][0];
    double g_tphi       = gcov[0][3];
    double g_phiphi     = gcov[3][3];
    double g_rr         = gcov[1][1];
    double g_thetatheta = gcov[2][2];

    double N_t     = sqrt(-g_tt + (g_tphi*g_tphi)/g_phiphi);
    double N_phi   = sqrt(g_phiphi);
    double N_r     = sqrt(g_rr);
    double N_theta = sqrt(g_thetatheta);

    e_t[0] = 1.0 / N_t;
    e_t[1] = 0.0;
    e_t[2] = 0.0;
    e_t[3] = g_tphi / (g_phiphi * N_t);

    e_r[0] = 0.0;
    e_r[1] = 1.0 / N_r;
    e_r[2] = 0.0;
    e_r[3] = 0.0;

    e_theta[0] = 0.0;
    e_theta[1] = 0.0;
    e_theta[2] = 1.0 / N_theta;
    e_theta[3] = 0.0;

    e_phi[0] = 0.0;
    e_phi[1] = 0.0;
    e_phi[2] = 0.0;
    e_phi[3] = 1.0 / N_phi;
}

int solve_geodesic_AVX(double X[NDIM], __m256d p[NDIM]) {
    double lambda_max = 20.0;  
    __m256d step_size = _mm256_set1_pd(0.00009);
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
    
    double E = -gcov[0][0] * p[0][0] - gcov[0][3] * p[3][0];
    double Lz = gcov[0][3] * p[0][0] + gcov[3][3] * p[3][0];

    double r_ps = 2 + 2 * cos( (2.0 / 3.0) * acos(-a) );
    double Lz_crit = - (r_ps * r_ps + a * a - 2 * r_ps) / a;
    double Q = ((Lz - a * E) * (Lz - a * E)) / (a * a) - r_ps * r_ps;

    for (int i = 0; i < num_points; i++) {
        double x = geodesic_points[i][0];
        double y = geodesic_points[i][1];
        double z = geodesic_points[i][2];

        double r = sqrt(x*x + y*y + z*z);


		if (r < r_horizon || Q < -1e-10  ) {
			printf("Photon absorbed at r = %f, E = %f, Lz = %f, Q = %f, r_horizon = %f, r_ps = %f\n", r, E, Lz, Q, r_horizon, r_ps);
			return 1;
		}
  
    }

    printf("Photon escaped, r = %f, E = %f, Lz = %f, Q = %f\n", sqrt(X[1]*X[1] + X[2]*X[2] + X[3]*X[3]), E, Lz, Q); 
    return 0; 
}


void generate_blackhole_shadow() {
    double r_obs = 5.0;
    double theta_obs = M_PI / 2.0;
    double phi_obs   = 0.0;

    double X[4] = {0.0, r_obs, theta_obs, phi_obs};

    double gcov[4][4], gcon[4][4];
    calculate_metric(X, gcov, gcon);

    
	double e_t_loc[4], e_r_loc[4], e_theta_loc[4], e_phi_loc[4];
	build_kerr_tetrad(X, gcov, e_t_loc, e_r_loc, e_theta_loc, e_phi_loc);

	double e_t[4], e_x[4], e_y[4], e_z[4];
	for (int mu = 0; mu < 4; mu++) {
		e_t[mu] = e_t_loc[mu];
		e_x[mu] = e_r_loc[mu];
		e_y[mu] = e_phi_loc[mu];
		e_z[mu] = e_theta_loc[mu];
	}

    double fov = 1.0;
    
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			int total_hit = 0;
			for (int sx = 0; sx < 2; sx++) {
				for (int sy = 0; sy < 2; sy++) {
					double alpha = fov * (2.0 * ((j + (sx / 2.0)) / WIDTH) - 1.0);
					double beta  = fov * (2.0 * ((i + (sy / 2.0)) / HEIGHT) - 1.0);

					double p[4];
					init_photon_global(alpha, beta, e_t, e_x, e_y, e_z, gcov, p);

					__m256d p_avx[4], X_avx[4];
					for(int mu = 0; mu < 4; mu++) {
						p_avx[mu] = _mm256_set1_pd(p[mu]);
						X_avx[mu] = _mm256_set1_pd(X[mu]); 
					}

					total_hit += solve_geodesic_AVX(X, p_avx);
				}
			}
			image[i][j] = (total_hit >= 2) ? 0 : 255;
		}
	}

    FILE *f = fopen("blackhole_shadow.ppm","wb");
    fprintf(f,"P5\n%d %d\n255\n", WIDTH, HEIGHT);
    fwrite(image, 1, WIDTH*HEIGHT, f);
    fclose(f);
    printf("Image saved: blackhole_shadow.ppm\n");
}
