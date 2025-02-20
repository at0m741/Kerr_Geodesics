
#include "matrix.h"
#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;
double Lambda = 1e-4;
double Q = 0.9;



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

    Matrix matrix_obj; 
    if (matrix_obj.inverse_3x3(gamma, gamma_inv) == 0) {
        printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
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
}


void calculate_dbeta(double X[3], double dbeta[3][3]) {
    memset(dbeta, 0, sizeof(double)*3*3);

    for(int m=0; m<3; m++){
        double Xp[3], Xm[3];
        memcpy(Xp, X, sizeof(Xp));
        memcpy(Xm, X, sizeof(Xm));

        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;

        double beta_p[3], beta_m[3];
        calculeBeta(Xp, beta_p); 
        calculeBeta(Xm, beta_m); 

        for(int j=0; j<3; j++){
            dbeta[m][j] = (beta_p[j] - beta_m[j]) / (2.0*DELTA3);
            printf("dbeta[%d][%d] = %e\n", m, j, dbeta[m][j]);
        }
    }
}


void calculeBeta(double X[3], double beta_cov[3]) {
	Metric metric_obj;
    std::array<double, 4> X4D = {0.0, X[0], X[1], X[2]};
    std::array<std::array<double, 4>, 4> gcov{}, gcon{};
    for(int i=0; i<3; i++){
        beta_cov[i] = gcov[0][i+1];
		printf("beta_cov[%d] = %e\n", i, beta_cov[i]);
    }
}

double compute_K(double gamma_inv[3][3], double K[3][3]) {
    double K_trace = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            K_trace += gamma_inv[i][j] * K[i][j];
        }
    }
    return K_trace;
}

double compute_Kij_Kij(double gamma_inv[3][3], double K[3][3]) {
    double K_sq = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    K_sq += gamma_inv[i][k] * gamma_inv[j][l] * K[i][j] * K[k][l];
                }
            }
        }
    }
    return K_sq;
}


double compute_hamiltonian_constraint(double gamma_inv[3][3], double K[3][3], double Ricci[3][3]) {
    double R = 0.0;
    double K_trace = compute_K(gamma_inv, K);
    double K_sq = compute_Kij_Kij(gamma_inv, K);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R += gamma_inv[i][j] * Ricci[i][j];
        }
    }

    return R + K_trace * K_trace - K_sq;
}

void compute_extrinsic_curvature_stationary_3D(
    double X[3],       
    double alpha,
    double beta_cov[3],
    double Gamma3[3][3][3],
    double dbeta[3][3],
    double K[3][3]
){

    double nabla[3][3]; 
    memset(nabla, 0, sizeof(nabla));

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            double partial = dbeta[i][j];
            double chris = 0.0;
            for(int k=0; k<3; k++){
                chris += Gamma3[k][i][j] * beta_cov[k];
            }
            nabla[i][j] = partial - chris;
        }
    }

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            double val = nabla[i][j] + nabla[j][i];
            K[i][j] = val / (2.0 * alpha);
        }
    }
}


void Metric::calculate_metric(const std::array<double, NDIM>& x, 
                              std::array<std::array<double, NDIM>, NDIM>& g,
                              std::array<std::array<double, NDIM>, NDIM>& g_inv) {
    Matrix matrix_obj;
    Metric metric_obj;
    double r = x[1];
    double theta = x[2];
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta2 = cos_theta * cos_theta;
    
    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2.0 * M * r + a * a;
    
    g[0][0] = -(1.0 - (2.0 * M * r) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r * a * a * sin_theta2) / Sigma);
    g[0][3] = - (2.0 * M * r * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);
}


void Metric::calculate_metric_kds(const std::array<double, NDIM>& x, 
                                  std::array<std::array<double, NDIM>, NDIM>& g,
                                  std::array<std::array<double, NDIM>, NDIM>& g_inv) {
    Matrix matrix_obj;
    Metric metric_obj;
    double r = x[1];
    double theta = x[2];
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double Sigma = r * r + a * a * cos_theta * cos_theta;
    double Delta_r = (1.0 - (Lambda * r * r) / 3.0) * (r * r + a * a) - 2.0 * M * r;
    double Delta_theta = 1.0 + (Lambda * a * a / 3.0) * cos_theta * cos_theta;
    double Xi = 1.0 - (Lambda * a * a) / 3.0;

    for (auto& row : g) {
        row.fill(0.0);
    }
    for (auto& row : g_inv) {
        row.fill(0.0);
    }

    g[0][0] = - (Delta_r / (Sigma * Xi * Xi));
    g[1][1] = Sigma / Delta_r;
    g[2][2] = Sigma / Delta_theta;
    g[3][3] = (sin_theta * sin_theta / (Sigma * Xi * Xi)) * (r * r + a * a) * (r * r + a * a);
    g[0][3] = - (2.0 * M * r * a * sin_theta * sin_theta) / (Sigma * Xi * Xi);
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);

    if (a == 0.0 && Lambda == 0.0)
        printf("Schwarzschild metric calculated\n");
    else if (a != 0.0 && Lambda == 0.0)
        printf("Kerr metric calculated\n");
    else 
        printf("Kerr-de Sitter metric calculated\n");

    matrix_obj.print_matrix("g", g);
}

void Metric::calculate_metric_kerr_newman(const std::array<double, NDIM>& x, 
				std::array<std::array<double, NDIM>, NDIM>& g,
				std::array<std::array<double, NDIM>, NDIM>& g_inv) {
	Matrix matrix_obj;
	Metric metric_obj;
	double r = x[1];
    double theta = x[2];
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta2 = cos_theta * cos_theta;

    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2.0 * M * r + a * a + Q * Q;

	for (auto& row : g) {
		row.fill(0.0);
	} 
	for (auto& row : g_inv) {
		row.fill(0.0);
	}

    g[0][0] = -(1.0 - (2.0 * M * r - Q * Q) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r - Q * Q) * a * a * sin_theta2 / Sigma);
    g[0][3] = -((2.0 * M * r - Q * Q) * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);

    if (a == 0.0 && Q == 0.0) 
        printf("Schwarzschild metric calculated\n");
    else if (Q == 0.0) 
        printf("Kerr metric calculated\n");
	else
        printf("Kerr-Newman metric calculated\n");

    matrix_obj.print_matrix("g", g);
    matrix_obj.print_matrix("g_inv", g_inv);
}

void Metric::verify_metric(const std::array<std::array<double, NDIM>, NDIM>& g,
				const std::array<std::array<double, NDIM>, NDIM>& g_inv){
	
	Matrix matrix_obj;
    int i, j, k;
    double identity[NDIM][NDIM] = {0};
    double delta; 

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            identity[i][j] = 0.0;
            for (k = 0; k < NDIM; k++) {
                identity[i][j] += gcon[i][k] * gcov[k][j];
            }
        }
    }

    for (i = 0; i < NDIM; i++) {
        for (j = 0; j < NDIM; j++) {
            if (i == j) {
                delta = 1.0;
            }
            else {
                delta = 0.0;
            }

            if (fabs(identity[i][j] - delta) > TOLERANCE) {
                printf("Erreur: identity[%d][%d] = %e with %e\n", i, j, identity[i][j], delta);
            }
        }
    }
	matrix_obj.check_inverse(gcov, gcon);
}

void verify3x3(const Matrix3x3& matrix, const Matrix3x3& inv_matrix) {
    Matrix matrix_obj;
    Matrix3x3 identity{};
    double delta = 0.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            identity[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                identity[i][j] += inv_matrix[i][k] * matrix[k][j];
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            delta = (i == j) ? 1.0 : 0.0;
            if (std::fabs(identity[i][j] - delta) > TOLERANCE) {
                printf("Erreur: identity[%d][%d] = %e with %e\n", i, j, identity[i][j], delta);
            }
        }
    }
    
    matrix_obj.print_matrix_3x3("identity", identity);
}
