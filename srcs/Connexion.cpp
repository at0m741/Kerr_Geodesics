#include <Geodesics.h>


void Connexion::calculate_christoffel(const VectorNDIM& X, double h,
                                      Christoffel3D& gamma,
                                      std::array<std::array<double, NDIM>, NDIM>& g,
                                      std::array<std::array<double, NDIM>, NDIM>& g_inv, 
                                      const char* metric) {
    Christoffel3D tmp{};
    VectorNDIM Xh = X, Xl = X;
    MatrixNDIM gh{}, gl{};
    Metric metric_obj;

    gamma.fill({}); 

    for (int mu = 0; mu < NDIM; mu++) {
        Xh = X;
        Xl = X;
        Xh[mu] += DELTA;
        Xl[mu] -= DELTA;

        if (strcmp(metric, "schwarzschild") == 0 || 
            strcmp(metric, "kerr") == 0 ||
            strcmp(metric, "kerr-newman") == 0 ||
            strcmp(metric, "ds") == 0) {
            
			metric_obj.calculate_metric(Xh, gh, g_inv); 
            metric_obj.calculate_metric(Xl, gl, g_inv);
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

    std::cout << "Christoffel symbols calculated\n";
    print_christoffel_matrix(gamma); 
    check_symmetry_christoffel(gamma);
}



static double dgamma[NDIM3][NDIM3][NDIM3];

void calc_gamma_ij(const Vector3& X3D, Matrix3x3& gamma3, Matrix3x3& gamma3_inv) {
    Vector4 X4D = { 0.0, X3D[0], X3D[1], X3D[2] };
    Matrix4x4 g{};    
    Matrix4x4 g_inv{}; 
    Metric metric;
    Matrix matrix_obj;
    metric.calculate_metric(X4D, g, g_inv);

    Matrix3x3 gamma3_arr{};
    for (int i = 0; i < NDIM3; i++) {
        for (int j = 0; j < NDIM3; j++) {
            gamma3_arr[i][j] = g[i+1][j+1];
        }
    }

    Matrix3x3 gamma3_inv_arr{};
    if (!matrix_obj.inverse_3x3(gamma3_arr, gamma3_inv_arr)) {
        printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
    }

    gamma3 = gamma3_arr;
    gamma3_inv = gamma3_inv_arr;
}


void Grid::calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3, 
                                     const Matrix3x3& gamma, Matrix3x3 gamma_inv) {
    Matrix3x3 gamma_p{}, gamma_m{};
    Tensor3D dgamma{};
    Gamma3.fill({});  

    for (int m = 0; m < NDIM3; m++) {
        Vector3 Xp = X;
        Vector3 Xm = X;
        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;

        calc_gamma_ij(Xp, gamma_p, gamma_inv);
        calc_gamma_ij(Xm, gamma_m, gamma_inv);
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
                    sum += gamma_inv[k][l] * ( dgamma[i][l][j] + dgamma[j][l][i] - dgamma[l][i][j] );
                }
                Gamma3[k][i][j] = 0.5 * sum;
            }
        }
    }

    const double tol = 1e-12;
    for (int k = 0; k < NDIM3; k++) {
        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                if (fabs(Gamma3[k][i][j] - Gamma3[k][j][i]) > tol) {
                    printf("Erreur: Gamma[%d][%d][%d] != Gamma[%d][%d][%d]\n", 
                           k, i, j, k, j, i);
                    printf("Gamma[%d][%d][%d] = %f, Gamma[%d][%d][%d] = %f\n",
                           k, i, j, Gamma3[k][i][j], k, j, i, Gamma3[k][j][i]);
                }
            }
        }
    }
}



void Grid::calculate_christoffel_3D_grid(
    std::vector<std::vector<Grid::Cell2D>>& grid,
    int Nx, int Ny,
    double dr, double dtheta,
    double r_min,
    double theta_min
) {
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            const Matrix3x3 &gamma = grid[i][j].gamma;
            const Matrix3x3 &gamma_inv = grid[i][j].gamma_inv;
            double dgamma[3][3][3] = {{{0.0}}};

            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    double gp = grid[i+1][j].gamma[a][b];
                    double gm = grid[i-1][j].gamma[a][b];
                    dgamma[0][a][b] = (gp - gm) / (2.0 * dr);
                }
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    double gp = grid[i][j+1].gamma[a][b];
                    double gm = grid[i][j-1].gamma[a][b];
                    dgamma[1][a][b] = (gp - gm) / (2.0 * dtheta);
                }
            }
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    dgamma[2][a][b] = 0.0;
                }
            }

            Tensor3D localGamma3 = {{{{0.0}}}};
            for (int k = 0; k < 3; k++) {
                for (int i_idx = 0; i_idx < 3; i_idx++) {
                    for (int j_idx = 0; j_idx < 3; j_idx++) {
                        double sum = 0.0;
                        for (int l = 0; l < 3; l++) {
                            double d_i_gamma_jl = (i_idx < 2) ? dgamma[i_idx][j_idx][l] : 0.0;
                            double d_j_gamma_il = (j_idx < 2) ? dgamma[j_idx][i_idx][l] : 0.0;
                            double d_l_gamma_ij = (l < 2)      ? dgamma[l][i_idx][j_idx] : 0.0;
                            sum += gamma_inv[k][l] * (d_i_gamma_jl + d_j_gamma_il - d_l_gamma_ij);
                        }
                        localGamma3[k][i_idx][j_idx] = 0.5 * sum;
                    }
                }
            }
            grid[i][j].Gamma3 = localGamma3;
        }
    }

#ifdef DEBUG
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						if (grid[i][j].Gamma3[k][l][m] != 0.0)
							printf("Gamma3[%d][%d][%d][%d][%d] = %e\n", i, j, k, l, m, grid[i][j].Gamma3[k][l][m]);

						if (fabs(grid[i][j].Gamma3[k][l][m] - grid[i][j].Gamma3[k][m][l]) > 1e-12) {
							printf("Erreur: Gamma3[%d][%d][%d][%d][%d] != Gamma3[%d][%d][%d][%d][%d]\n", 
									i, j, k, l, m, i, j, k, m, l);
							printf("Gamma3[%d][%d][%d][%d][%d] = %f, Gamma3[%d][%d][%d][%d][%d] = %f\n",
									i, j, k, l, m, grid[i][j].Gamma3[k][l][m],
									i, j, k, m, l, grid[i][j].Gamma3[k][m][l]);
						}
					}
				}
			}
		}
	}
#endif
}
