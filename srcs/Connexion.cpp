#include <Geodesics.h>


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
    Vector4 X4D = {0.0, X3D[0], X3D[1], X3D[2]};
    
    Matrix4x4 g{};    
    Matrix4x4 g_inv{}; 

    Metric metric;
    Matrix matrix_obj;

    metric.calculate_metric(X4D, g, g_inv);

    Matrix3x3 gamma3_arr{};
    Matrix3x3 gamma3_inv_arr{};
    for (int i = 0; i < NDIM3; i++) {
        for (int j = 0; j < NDIM3; j++) {
            gamma3_arr[i][j] = g[i+1][j+1];
        }
    }

    if (!matrix_obj.inverse_3x3(gamma3_arr, gamma3_inv_arr)) {
        printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
    }

    gamma3 = gamma3_arr;
    gamma3_inv = gamma3_inv_arr;
}


void Grid::calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3) {
    Matrix3x3 gamma_ij{}, gamma_inv{};
    calc_gamma_ij(X, gamma_ij, gamma_inv);
	 
    for (int m = 0; m < NDIM3; m++) {
        Vector3 Xp = X;
        Vector3 Xm = X;
        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;
        
        Matrix3x3 gamma_p{}, gamma_m{}, dummy_inv;
        calc_gamma_ij(Xp, gamma_p, dummy_inv);
        calc_gamma_ij(Xm, gamma_m, dummy_inv);
        
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
}


void Grid::calculate_christoffel_3D_grid(
    std::vector<std::vector<Cell2D>>& grid,
    int Nx, int Ny,
    double dx, double dy
){
    for(int i=1; i<Nx-1; i++){
        for(int j=1; j<Ny-1; j++){
            const auto& gamma_ij = grid[i][j].gamma;
            const auto& gamma_inv = grid[i][j].gamma_inv;

            double dgamma[3][3][3] = {{{0.0}}};

            for(int a=0; a<3; a++){
                for(int b=0; b<3; b++){
                    double g_p = grid[i+1][j].gamma[a][b];
                    double g_m = grid[i-1][j].gamma[a][b];
                    dgamma[0][a][b] = (g_p - g_m)/(2.0 * dx);
                }
            }

            for(int a=0; a<3; a++){
                for(int b=0; b<3; b++){
                    double g_p = grid[i][j+1].gamma[a][b];
                    double g_m = grid[i][j-1].gamma[a][b];
                    dgamma[1][a][b] = (g_p - g_m)/(2.0 * dy);
                }
            }

			Tensor3D localGamma3 = {{{{0.0}}}};

            for(int k=0; k<3; k++){
                for(int iidx=0; iidx<3; iidx++){
                    for(int jidx=0; jidx<3; jidx++){
                        double sum = 0.0;
                        for(int l=0; l<3; l++){
                            double val_d1 = dgamma[iidx][l][jidx]; // ∂i gamma_{l j}
                            double val_d2 = dgamma[jidx][l][iidx]; // ∂j gamma_{l i}
                            double val_d3 = dgamma[l][iidx][jidx]; // ∂l gamma_{i j}
                            sum += gamma_inv[k][l] * (val_d1 + val_d2 - val_d3);
                        }
                        localGamma3[k][iidx][jidx] = 0.5 * sum;
						if (localGamma3[k][iidx][jidx] != 0.0) 
							printf("localGamma3[%d][%d][%d] = %e\n", k, iidx, jidx, localGamma3[k][iidx][jidx]);
                    }
                }
            }

            grid[i][j].Gamma3 = localGamma3;
        }
    }
}

