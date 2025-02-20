
#include "matrix.h"
#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;
double Lambda = 1e-4;
double Q = 0.9;

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

