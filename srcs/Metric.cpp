#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;
double Lambda = 1e-4;
double Q = 0.9;


void print_matrix(char *name, double matrix[NDIM][NDIM]) {
	int i, j;
	printf("%s = \n", name);
	for (i = 0; i < NDIM; i++) {
		for (j = 0; j < NDIM; j++) {
			printf("%e ", matrix[i][j]);
		}
		printf("\n");
	}
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
	verify_metric(g, g_inv);	
	if (a == 0.0){
		printf("Schwarzschild metric calculated\n");
	}
	else {
		printf("Kerr metric calculated\n");
	}
	/* print_matrix("g", g); */
	/* print_matrix("g_inv", g_inv); */
}

void calculate_metric_kds(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[1];
    double theta = x[2];
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    double Sigma = r * r + a * a * cos_theta * cos_theta;
    double Delta_r = (1.0 - (Lambda * r * r) / 3.0) * (r * r + a * a) - 2.0 * M * r;
    double Delta_theta = 1.0 + (Lambda * a * a / 3.0) * cos_theta * cos_theta;
    double Xi = 1.0 - (Lambda * a * a) / 3.0;

    memset(g, 0, sizeof(double) * NDIM * NDIM);
    memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

    g[0][0] = - (Delta_r / (Sigma * Xi * Xi));
    g[1][1] = Sigma / Delta_r;
    g[2][2] = Sigma / Delta_theta;
    g[3][3] = (sin_theta * sin_theta / (Sigma * Xi * Xi)) * (r * r + a * a) * (r * r + a * a);
    g[0][3] = - (2.0 * M * r * a * sin_theta * sin_theta) / (Sigma * Xi * Xi);
    g[3][0] = g[0][3];

    inverse_matrix(g, g_inv);
	verify_metric(g, g_inv);
	if (a == 0.0 && Lambda == 0.0){
		printf("Schwarzschild metric calculated\n");
	}
	else if (a != 0.0 && Lambda == 0.0){
		printf("Kerr metric calculated\n");
	}
	else {
		printf("Kerr-de Sitter metric calculated\n");
	}
	print_matrix("g", g);
}

void calculate_metric_kerr_newman(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]) {
    double r = x[1];
    double theta = x[2];

    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta2 = cos_theta * cos_theta;

    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2.0 * M * r + a * a + Q * Q;

    memset(g, 0, sizeof(double) * NDIM * NDIM);
    memset(g_inv, 0, sizeof(double) * NDIM * NDIM);

    g[0][0] = -(1.0 - (2.0 * M * r - Q * Q) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r - Q * Q) * a * a * sin_theta2 / Sigma);
    g[0][3] = -((2.0 * M * r - Q * Q) * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    inverse_matrix(g, g_inv);
	check_inverse(g, g_inv);

    if (a == 0.0 && Q == 0.0) {
        printf("Schwarzschild metric calculated\n");
    } else if (Q == 0.0) {
        printf("Kerr metric calculated\n");
    } else {
        printf("Kerr-Newman metric calculated\n");
    }

    print_matrix("g", g);
    print_matrix("g_inv", g_inv);
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


