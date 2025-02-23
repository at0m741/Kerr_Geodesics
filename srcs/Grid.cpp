#include <Geodesics.h>


void Grid::extract_3p1(const Matrix4x4& g,
                 const Matrix4x4& ,
                 double* alpha,
                 Vector3& beta_cov,
                 Vector3& beta_con,
                 Matrix3x3& gamma,
                 Matrix3x3& gamma_inv) {

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            gamma[i][j] = g[i+1][j+1];
        }
    }

    double det_gamma = gamma[0][0] * (gamma[1][1] * gamma[2][2] - gamma[1][2] * gamma[2][1]) -
                       gamma[0][1] * (gamma[1][0] * gamma[2][2] - gamma[1][2] * gamma[2][0]) +
                       gamma[0][2] * (gamma[1][0] * gamma[2][1] - gamma[1][1] * gamma[2][0]);


	Matrix matrix_obj; 
	if (matrix_obj.inverse_3x3(gamma, gamma_inv) == 0) {
		printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
	} else {
		printf("gamma_{ij} est inversible\n");
	}


    for (int i = 0; i < 3; i++) {
        beta_cov[i] = g[0][i+1];
    }

    for (int i = 0; i < 3; i++) {
        double sum = 0.0;
        for (int j = 0; j < 3; j++) {
            sum += gamma_inv[i][j] * beta_cov[j];
        }
        beta_con[i] = sum;
    }

    double betabeta = 0.0;
    for (int i = 0; i < 3; i++) {
        betabeta += beta_con[i] * beta_cov[i];
    }

    *alpha = sqrt(betabeta - g[0][0]);
    matrix_obj.print_matrix_3x3("gamma", gamma);
    matrix_obj.print_matrix_3x3("gamma_inv", gamma_inv);
    printf("beta_i = (%e, %e, %e)\n", beta_cov[0], beta_cov[1], beta_cov[2]);
	printf("beta^i = (%e, %e, %e)\n", beta_con[0], beta_con[1], beta_con[2]);
	printf("alpha = %e\n", *alpha);
}


void Grid::calculate_dbeta(const Vector3& X, Matrix3x3& dbeta) {
    for (auto& row : dbeta) {
        row.fill(0.0);
    }


    for (int m = 0; m < DIM3; m++) {
        Vector3 Xp = X;
        Vector3 Xm = X;
        
        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;
        
        Vector3 beta_p, beta_m;
        calculeBeta(Xp, beta_p); 
        calculeBeta(Xm, beta_m);
        
        for (int j = 0; j < DIM3; j++) {
            dbeta[m][j] = (beta_p[j] - beta_m[j]) / (2.0 * DELTA3);
            printf("dbeta[%d][%d] = %e\n", m, j, dbeta[m][j]);
        }
    }
}


void Grid::calculeBeta(const Vector3& X, Vector3& beta_cov) {
    std::array<double, 4> X4D = {0.0, X[0], X[1], X[2]};
    std::array<std::array<double, 4>, 4> gcov{};
    std::array<std::array<double, 4>, 4> gcon{};
    
    Metric metric_obj;
    metric_obj.calculate_metric(X4D, gcov, gcon);
    
    for (int i = 0; i < DIM3; i++){
        beta_cov[i] = gcov[0][i+1];
        printf("beta_cov[%d] = %e\n", i, beta_cov[i]);
    }
}

double Grid::compute_K(const Matrix3x3& gamma_inv, const Matrix3x3& K) {
    double K_trace = 0.0;
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            K_trace += gamma_inv[i][j] * K[i][j];
        }
    }
    return K_trace;
}

double Grid::compute_Kij_Kij(const Matrix3x3& gamma_inv, const Matrix3x3& K) {
    double K_sq = 0.0;
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            for (int k = 0; k < DIM3; k++) {
                for (int l = 0; l < DIM3; l++) {
                    K_sq += gamma_inv[i][k] * gamma_inv[j][l] * K[i][j] * K[k][l];
                }
            }
        }
    }
    return K_sq;
}


Matrix3x3 Grid::compute_beta_gradient(int i, int j) {
    Matrix3x3 beta_grad;
    for (auto& row : beta_grad) {
        row.fill(0.0);
    }

    beta_grad[0][0] = (grid[i+1][j].beta_con[0] - grid[i-1][j].beta_con[0]) / (2.0 * dr);
    beta_grad[1][1] = (grid[i][j+1].beta_con[1] - grid[i][j-1].beta_con[1]) / (2.0 * dtheta);

    return beta_grad;
}

Matrix3x3 Grid::compute_second_derivative_alpha(int i, int j) {
    Matrix3x3 d2alpha;
    for (auto& row : d2alpha) {
        row.fill(0.0);
    }

    d2alpha[0][0] = (grid[i+1][j].alpha - 2.0 * grid[i][j].alpha + grid[i-1][j].alpha) / (dr * dr);
    d2alpha[1][1] = (grid[i][j+1].alpha - 2.0 * grid[i][j].alpha + grid[i][j-1].alpha) / (dtheta * dtheta);

    return d2alpha;
}



double Grid::compute_KijKij_component(const Matrix3x3& gamma_inv, const Matrix3x3& K, int a, int b) {
    double sum = 0.0;
    for (int k = 0; k < DIM3; k++) {
        sum += gamma_inv[a][k] * K[k][b];
		printf("sum = %e\n", sum);
    }
    return sum;
}

void Grid::evolve_Kij(double dt) {
    for (int i = 1; i < Nr - 1; i++) {
        for (int j = 1; j < Ntheta - 1; j++) {
            Cell2D& cell = grid[i][j];

            double alpha = cell.alpha;
            Vector3 beta_con = cell.beta_con;
            Matrix3x3& K = cell.K;
            Matrix3x3& Ricci = cell.Ricci;

            double K_trace = compute_K(cell.gamma_inv, K);
            Matrix3x3 d2alpha = compute_second_derivative_alpha(i, j);
            Matrix3x3 beta_grad = compute_beta_gradient(i, j);
		
            Matrix3x3 K_new;
            for (int a = 0; a < DIM3; a++) {
                for (int b = 0; b < DIM3; b++) {
                    double term1 = -d2alpha[a][b];
					double term2 = alpha * (Ricci[a][b] + K_trace * K[a][b] - 2.0 * compute_KijKij_component(cell.gamma_inv, K, a, b));
                    double term3 = beta_con[0] * (grid[i+1][j].K[a][b] - grid[i-1][j].K[a][b]) / (2.0 * dr) +
                                   beta_con[1] * (grid[i][j+1].K[a][b] - grid[i][j-1].K[a][b]) / (2.0 * dtheta);
                    double term4 = K[a][0] * beta_grad[0][b] + K[a][1] * beta_grad[1][b] + K[a][2] * beta_grad[2][b];

                    K_new[a][b] = K[a][b] + dt * (term1 + term2 + term3 + term4);
					printf("K_new[%d][%d] = %e\n", a, b, K_new[a][b]);
                }
            }

            cell.K = K_new;
        }
    }
}

double Grid::compute_hamiltonian_constraint(const Matrix3x3& gamma_inv, const Matrix3x3& K, const Matrix3x3& Ricci) {
    double R = 0.0;
	Grid grid_obj;
    double K_trace = grid_obj.compute_K(gamma_inv, K);
    double K_sq = grid_obj.compute_Kij_Kij(gamma_inv, K);
    
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            R += gamma_inv[i][j] * Ricci[i][j];
        }
    }
	printf("Ricci = %e\n", R); 
	printf("K_trace = %e\n", K_trace);
	printf("K_sq = %e\n", K_sq);
	printf("R + K_trace * K_trace - K_sq = %e\n", R + K_trace * K_trace - K_sq);
    return R + K_trace * K_trace - K_sq;
}

void Grid::compute_extrinsic_curvature_stationary_3D(
    const Vector3& X,     
    double alpha,
    const Vector3& beta_cov,
    std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>& Gamma3,
    const Matrix3x3& dbeta,
    Matrix3x3& K)
{
    Matrix3x3 nabla{};
    for (auto& row : nabla) {
        row.fill(0.0);
    }
    
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            double partial = dbeta[i][j];
            double chris = 0.0;
            for (int k = 0; k < DIM3; k++) {
                chris += Gamma3[k][i][j] * beta_cov[k];
            }
            nabla[i][j] = partial - chris;
        }
    }
    
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            double val = nabla[i][j] + nabla[j][i];
            K[i][j] = val / (2.0 * alpha);
        }
    }
}
