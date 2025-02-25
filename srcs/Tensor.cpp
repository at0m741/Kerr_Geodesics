#include <Geodesics.h>

/* 
 * Richardson Extrapolation for the derivative of christoffel symbols
 * It is used to calculate the gamma derivative in the Riemann tensor
 */

double Tensor::richardson_derivative(
    const Tensor3D& Gamma_plus_h, 
    const Tensor3D& Gamma_minus_h,
    const Tensor3D& Gamma_plus_half_h,
    const Tensor3D& Gamma_minus_half_h,
    int rho, int mu, int nu, double h) 
{
    double diff_h = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / (2 * h);
    double diff_half_h = (Gamma_plus_half_h[rho][mu][nu] - Gamma_minus_half_h[rho][mu][nu]) / h;
    return (4 * diff_half_h - diff_h) / 3;
}

void Tensor::calculate_Gamma_at_offset(const std::array<double, NDIM>& X, int direction, 
                                         double offset, double delta,
                                         Tensor::MatrixNDIM& gcov, 
                                         Tensor::MatrixNDIM& gcon, 
                                         Tensor::Christoffel3D& Gamma_slice, 
                                         const char* metric_type) {
    std::array<double, NDIM> X_offset = X;
    X_offset[direction] += offset;
    Tensor::Christoffel3D tempGamma{};
    Connexion connexion;
    Metric metric;
    
    if (strcmp(metric_type, "minkowski") == 0) {
        connexion.calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, "minkowski");
    } else if (strcmp(metric_type, "kerr") == 0 || strcmp(metric_type, "schwarzschild") == 0) {
        metric.calculate_metric(X_offset, gcov, gcon);
    }
    connexion.calculate_christoffel(X_offset, delta, tempGamma, gcov, gcon, metric_type);
    
    Gamma_slice = tempGamma;
}



/* 
 * Calculate the Riemann tensor using the Christoffel symbols
 * The Riemann tensor is calculated using the formula:
 * R^rho_sigma_mu_nu = dGamma^rho_mu_nu/dx^sigma - dGamma^rho_nu/sigma
 * + Gamma^rho_mu_lambda * Gamma^lambda_nu_sigma - Gamma^rho_nu_lambda * Gamma^lambda_mu_sigma
 */

void Tensor::calculate_riemann(const Christoffel3D& Gamma, 
				const Christoffel4D& Gamma_plus_h, 
				const Christoffel4D& Gamma_minus_h, 
				const Christoffel4D& Gamma_plus_half_h, 
				const Christoffel4D& Gamma_minus_half_h,
				Riemann4D& Riemann, 
				double h) {
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    double dGamma_mu = richardson_derivative(
                        Gamma_plus_h[mu], Gamma_minus_h[mu], 
                        Gamma_plus_half_h[mu], Gamma_minus_half_h[mu], 
                        rho, nu, sigma, h);
                    
                    double dGamma_nu = richardson_derivative(
                        Gamma_plus_h[nu], Gamma_minus_h[nu],
                        Gamma_plus_half_h[nu], Gamma_minus_half_h[nu],
                        rho, mu, sigma, h);

                    double Gamma_terms = 0.0;
                    for (int lambda = 0; lambda < NDIM; lambda++) {
                        Gamma_terms += Gamma[rho][mu][lambda] * Gamma[lambda][nu][sigma]
                                     - Gamma[rho][nu][lambda] * Gamma[lambda][mu][sigma];
                    }

                    Riemann[rho][sigma][mu][nu] = dGamma_mu - dGamma_nu + Gamma_terms;
                }
            }
        }
    }
	double Kretschmann_scalar = 0.0;
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    Kretschmann_scalar += Riemann[rho][sigma][mu][nu] * \
										  Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
    if (Kretschmann_scalar > 1e10) {
        printf("Kretschmann Scalar: INF (Singularity detected)\n");
    } else {
        printf("Kretschmann Scalar: %12.6f\n", Kretschmann_scalar);
    }
}

/* 
 * Contract the Riemann tensor to calculate the Ricci tensor
 * The Ricci tensor is calculated using the formula:
 * R_mu_nu = g^rho_sigma * R^sigma_rho_mu_nu
 */

void Tensor::contract_riemann(const Riemann4D& Riemann, MatrixNDIM& Ricci, const MatrixNDIM& g_inv) {
    for (auto &row : Ricci) {
        row.fill(0.0);
    }
    
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    Ricci[mu][nu] += g_inv[rho][sigma] * Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
    
    printf("\nRicci tensor:\n");
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            printf("%12.6f\t", Ricci[mu][nu]);
        }
        printf("\n");
    }
    
    double Ricci_scalar = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            Ricci_scalar += g_inv[mu][nu] * Ricci[mu][nu];
        }
    }
    printf("Ricci Scalar: %12.6f\n", Ricci_scalar);
}




double Grid::richardson_derivative_ricci(
    const Tensor3D &Gamma_plus_h, 
    const Tensor3D &Gamma_minus_h,
    const Tensor3D &Gamma_plus_half_h,
    const Tensor3D &Gamma_minus_half_h,
    int mu, int nu, int sigma, 
    double h) 
{
    double diff_h = (Gamma_plus_h[mu][nu][sigma] - Gamma_minus_h[mu][nu][sigma]) / (2.0 * h);
    double diff_half_h = (Gamma_plus_half_h[mu][nu][sigma] - Gamma_minus_half_h[mu][nu][sigma]) / h;
    return (4.0 * diff_half_h - diff_h) / 3.0;
}



void Grid::calculate_riemann_3d(
    const Christoffel3D& Gamma, 
    const std::array<Christoffel3D, 3>& Gamma_plus_h,
    const std::array<Christoffel3D, 3>& Gamma_minus_h,
    const std::array<Christoffel3D, 3>& Gamma_plus_half_h,
    const std::array<Christoffel3D, 3>& Gamma_minus_half_h,
    Riemann3D& Riemann,
    double h,
    double scale  
) {
    double effective_h = h * scale;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double dGamma_k = richardson_derivative_ricci(
                        Gamma_plus_h[k], Gamma_minus_h[k],
                        Gamma_plus_half_h[k], Gamma_minus_half_h[k],
                        i, j, l, 
                        effective_h
                    );

                    double dGamma_l = richardson_derivative_ricci(
                        Gamma_plus_h[l], Gamma_minus_h[l],
                        Gamma_plus_half_h[l], Gamma_minus_half_h[l],
                        i, j, k, 
                        effective_h
                    );

                    double Gamma_terms = 0.0;
                    for (int m = 0; m < 3; m++) {
                        Gamma_terms += Gamma[i][k][m] * Gamma[m][l][j]
                                      - Gamma[i][l][m] * Gamma[m][k][j];
                    }
                    
                    Riemann[i][j][k][l] = dGamma_k - dGamma_l + Gamma_terms;
                }
            }
        }
    }
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				printf("Riemann[%d][%d][%d][%d] = %e\n", i, j, k, k, Riemann[i][j][k][k]);
			}
		}
	}
}

void Grid::calculate_riemann_4d_from_3d(
    const Riemann3D &Riemann3,
    const Matrix3x3 &K,     
    Riemann3D &Riemann4    
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    Riemann4[i][j][k][l] = Riemann3[i][j][k][l]
                                          + K[i][k] * K[j][l]
                                          - K[i][l] * K[j][k];

                    if (i == j || k == l) {
                        Riemann4[i][j][k][l] = 0.0;  
                    }
                }
            }
        }
    }
    verify_riemann_symmetries(Riemann4);  
}


void Grid::calculate_ricci_3d_from_riemann(const Riemann3D& Riemann, Matrix3x3& Ricci) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += Riemann[k][i][k][j];
            }
            Ricci[i][j] = sum;
        }
    }
    printf("\nRicci tensor:\n");
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            printf("%12.6f\t", Ricci[mu][nu]);
        }
        printf("\n");
    }
    

}


/* void Grid::compute_ricci_3d_grid( */
/*     std::vector<std::vector<Cell2D>>& grid, */
/*     int Nr, int Ntheta, */
/*     double dr, double dtheta, */
/*     double r_min, double theta_min, */
/*     double delta   */
/* ) { */
/*     double h_riemann = delta;   */
/*  */
/*     for (int i = 0; i < Nr; i++) { */
/*         for (int j = 0; j < Ntheta; j++) { */
/*             double r_i = r_min + i * dr; */
/*             double th_j = theta_min + j * dtheta; */
/*             Vector3 X3D = { r_i, th_j, 0.0 }; */
/*  */
/*             Christoffel3D Gamma = grid[i][j].Gamma3; */
/*             Christoffel3D Gamma_plus_h, Gamma_minus_h; */
/*             Christoffel3D Gamma_plus_half_h, Gamma_minus_half_h; */
/*  */
/*             { */
/*                 Vector3 Xp = X3D, Xm = X3D, Xp_half = X3D, Xm_half = X3D; */
/*                 Xp[0] += h_riemann; */
/*                 Xm[0] -= h_riemann; */
/*                 Xp_half[0] += h_riemann / 2.0; */
/*                 Xm_half[0] -= h_riemann / 2.0; */
/*  */
/*                 calculate_christoffel_3D(Xp, Gamma_plus_h, grid[i][j].gamma, grid[i][j].gamma_inv); */
/*                 calculate_christoffel_3D(Xm, Gamma_minus_h, grid[i][j].gamma, grid[i][j].gamma_inv); */
/*                 calculate_christoffel_3D(Xp_half, Gamma_plus_half_h, grid[i][j].gamma, grid[i][j].gamma_inv); */
/*                 calculate_christoffel_3D(Xm_half, Gamma_minus_half_h, grid[i][j].gamma, grid[i][j].gamma_inv); */
/*             } */
/*  */
/*             Riemann3D Riemann; */
/*             calculate_riemann_3d(Gamma, Gamma_plus_h, Gamma_minus_h, */
/*                                  Gamma_plus_half_h, Gamma_minus_half_h, */
/*                                  Riemann, h_riemann, DELTA); */
/*  */
/*             Matrix3x3 Ricci; */
/*             calculate_ricci_3d_from_riemann(Riemann, Ricci); */
/*  */
/*             grid[i][j].Ricci = Ricci; */
/*             printf("Ricci tensor computed at r=%e, theta=%e\n", r_i, th_j); */
/*         } */
/*     } */
/* } */



void Grid::print_ricci_tensor(const Matrix3x3& R3) {
    printf("\nRicci tensor:\n");
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            printf("%12.6f\t", R3[i][j]);
        }
        printf("\n");
    }
}



void compute_partial_christoffel_3D(
    const Vector3& X,   
    int m,             
    Tensor3D& dGamma,  
    double delta
) {
    Vector3 Xp = X;
    Vector3 Xm = X;
    Grid grid_obj; 
    Xp[m] += delta;
    Xm[m] -= delta;
    
    Tensor3D Gamma_p{}; 
    Tensor3D Gamma_m{};
    Matrix3x3 gamma{};   
    Matrix3x3 gamma_inv{};

    grid_obj.calculate_christoffel_3D(Xp, Gamma_p, gamma, gamma_inv);
    grid_obj.calculate_christoffel_3D(Xm, Gamma_m, gamma, gamma_inv);

    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dGamma[k][i][j] = (Gamma_p[k][i][j] - Gamma_m[k][i][j]) / (2.0 * delta);
            }
        }
    }
}



void Grid::compute_ricci_3d(
    const Vector3& X,       
    const Tensor3D& Gamma3, 
    Matrix3x3& R3   
) {
    for (auto &row : R3) {
        row.fill(0.0);
    }

    static Tensor4D partialGamma{}; 

    double delta = 1e-5;

    for (int m = 0; m < 3; m++) {
        std::array<std::array<std::array<double, 3>, 3>, 3> dG{};
        compute_partial_christoffel_3D(X, m, dG, delta);
        
        for (int k = 0; k < 3; k++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    partialGamma[m][k][i][j] = dG[k][i][j];
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            for (int k = 0; k < 3; k++) {
                term1 += partialGamma[k][k][i][j];
                term2 += partialGamma[j][k][i][k];
            }
            for (int k = 0; k < 3; k++) {
                for (int m = 0; m < 3; m++) {
                    term3 += Gamma3[k][i][j] * Gamma3[m][k][m]; 
                }
            }
            for (int m = 0; m < 3; m++) {
                for (int k = 0; k < 3; k++) {
                    term4 += Gamma3[m][i][k] * Gamma3[k][j][m]; 
                }
            }
            R3[i][j] = term1 - term2 + term3 - term4;
        }
    }

    print_ricci_tensor(R3);
}

