#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

void minkowski_metric(double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
	memset(g, 0, sizeof(double) * NDIM * NDIM);
	memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

	g[0][0] = -1.0;	
	g[1][1] = 1.0;
	g[2][2] = 1.0;
	g[3][3] = 1.0;

	inverse_matrix(g, g_inv);
	printf("Minkowski metric calculated\n");

}

void calculate_metric(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[1];
    double theta = x[2];
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta2 = cos_theta * cos_theta;
    
    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2.0 * M * r + a * a;
    
    memset(g, 0, sizeof(double) * NDIM * NDIM);
    memset(g_inv, 0, sizeof(double) * NDIM * NDIM);
    
    g[0][0] = -(1.0 - (2.0 * M * r) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r * a * a * sin_theta2) / Sigma);
    g[0][3] = - (2.0 * M * r * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];
    inverse_matrix(g, g_inv);
}



void verify_metric(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM])
{
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
	printf("\n");
	check_inverse(gcov, gcon);
	printf("\n");
	print_matrix("gcov", gcov);
	print_matrix("gcon", gcon);
	printf("\n");
	printf("Identity matrix:\n");
	print_matrix("identity", identity);

}


