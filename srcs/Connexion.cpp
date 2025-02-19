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
		if (strcmp(metric, "schwarzschild") == 0 || 
			strcmp(metric, "kerr") == 0 ||
			strcmp(metric, "kerr-newman") == 0 ||
			strcmp(metric, "ds") == 0) {
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



static double dgamma[NDIM3][NDIM3][NDIM3];

void calc_gamma_ij(const double X3D[3],
                   double gamma3[3][3],       
                   double gamma3_inv[3][3])  
{
    double X4D[4] = {0.0, X3D[0], X3D[1], X3D[2]};

    double g[4][4], g_inv[4][4];
    memset(g, 0, sizeof(g));
    memset(g_inv, 0, sizeof(g_inv));

    calculate_metric(X4D, g, g_inv);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            gamma3[i][j] = g[i+1][j+1];
        }
    }

    int status = inverse_3x3(gamma3, gamma3_inv);
    if (!status) {
        printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
    }
}


void calculate_christoffel_3D(
    double X[NDIM3],         
    double Gamma3[NDIM3][NDIM3][NDIM3] 
) {
    double gamma_ij[NDIM3][NDIM3], gamma_inv[NDIM3][NDIM3];
    calc_gamma_ij(X, gamma_ij, gamma_inv);

    for (int m = 0; m < NDIM3; m++) {
        double Xp[NDIM3], Xm[NDIM3];
        memcpy(Xp, X, sizeof(Xp));
        memcpy(Xm, X, sizeof(Xm));

        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;

        double gamma_p[NDIM3][NDIM3], gamma_m[NDIM3][NDIM3];
        double gamma_p_inv[NDIM3][NDIM3], gamma_m_inv[NDIM3][NDIM3];
        calc_gamma_ij(Xp, gamma_p, gamma_p_inv);
        calc_gamma_ij(Xm, gamma_m, gamma_m_inv);

        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                dgamma[m][i][j] = (gamma_p[i][j] - gamma_m[i][j]) / (2.0 * DELTA3);
            }
        }
    }

    for (int k = 0; k < NDIM3; k++) {
        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                double sum = 0.0;
                for (int l = 0; l < NDIM3; l++) {
                    sum += gamma_inv[k][l] * (
                        dgamma[i][l][j] + dgamma[j][l][i] - dgamma[l][i][j]
                    );
                }
                Gamma3[k][i][j] = 0.5 * sum;
            }
        }
    }
}

