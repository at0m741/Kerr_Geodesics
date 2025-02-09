#include <Geodesics.h>



void calculate_Gamma_at_offset(double X[NDIM], int direction, 
						double offset, double delta,
						double gcov[NDIM][NDIM], 
						double gcon[NDIM][NDIM], 
						double Gamma_slice[NDIM][NDIM][NDIM], 
						const char *metric_type) {
    double X_offset[NDIM];
    memcpy(X_offset, X, sizeof(double) * NDIM);
    X_offset[direction] += offset;
    
    double tempGamma[NDIM][NDIM][NDIM];
	if (strcmp(metric_type, "minkowski") == 0) {
		calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, "minkowski");
	} else if (strcmp(metric_type, "kerr") == 0 || strcmp(metric_type, "schwarzschild") == 0) {
		calculate_metric(X_offset, gcov, gcon);
	}
    calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, metric_type);
    
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                Gamma_slice[rho][mu][nu] = tempGamma[rho][mu][nu];
            }
        }
    }
}

void calculate_christoffel(double X[NDIM], double h, \
							double gamma[NDIM][NDIM][NDIM],
							double g[NDIM][NDIM],
							double g_inv[NDIM][NDIM], const char *metric) {
    double tmp[NDIM][NDIM][NDIM];
    double Xh[NDIM], Xl[NDIM]; 
    double gh[NDIM][NDIM], gl[NDIM][NDIM]; 

    memset(gamma, 0, sizeof(double) * NDIM * NDIM * NDIM);

    for (int mu = 0; mu < NDIM; mu++) {
        for (int kap = 0; kap < NDIM; kap++) {
            Xh[kap] = X[kap];
            Xl[kap] = X[kap];
        }

        Xh[mu] += DELTA;
        Xl[mu] -= DELTA;
		if (strcmp(metric, "minkowski") == 0)
			minkowski_metric(gh, g_inv);
		else if (strcmp(metric, "schwarzschild") == 0 || strcmp(metric, "kerr") == 0) {
			calculate_metric(Xh, gh, g_inv);
			calculate_metric(Xl, gl, g_inv);
			verify_metric(gh, g_inv);
			verify_metric(gl, g_inv);
		}
		else {
			printf("Invalid metric type\n");
			return;
		}
        for (int lam = 0; lam < NDIM; lam++) {
            for (int nu = 0; nu < NDIM; nu++) {
                gamma[lam][nu][mu] = (gh[lam][nu] - gl[lam][nu]) / (2 * DELTA);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                tmp[lam][nu][mu] = 0.5 * (gamma[nu][lam][mu] + 
                                          gamma[mu][lam][nu] - 
                                          gamma[mu][nu][lam]);
            }
        }
    }

    for (int lam = 0; lam < NDIM; lam++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int mu = 0; mu < NDIM; mu++) {
                gamma[lam][nu][mu] = 0.0;
                for (int kap = 0; kap < NDIM; kap++) {
                    gamma[lam][nu][mu] += g_inv[lam][kap] * tmp[kap][nu][mu];
                }
            }
        }
    }    
	printf("Christoffel symbols calculated\n");
	print_christoffel_matrix(gamma); 
	check_symmetry_christoffel(gamma);
}


