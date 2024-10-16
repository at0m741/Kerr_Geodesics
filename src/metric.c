/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metric.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/22 14:44:57 by ltouzali          #+#    #+#             */
/*   Updated: 2024/10/10 00:09:41 by babonnet         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"
extern double		(*geodesic_points)[5];
extern int			num_points;
int					i, j, k;

double determinant3x3(double mat[3][3]) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
         - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
         + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

double determinant4x4(double mat[NDIM][NDIM]) {
    double minor[3][3];
    double det = 0.0;
    for (int i = 0; i < NDIM; i++) {
        int subi = 0; 
        for (int j = 1; j < NDIM; j++) {
            int subj = 0;
            for (int k = 0; k < NDIM; k++) {
                if (k == i)
                    continue;
                minor[subi][subj] = mat[j][k];
                subj++;
            }
            subi++;
        }
        det += pow(-1, i) * mat[0][i] * determinant3x3(minor);
    }
    return det;
}

void cofactor(double mat[NDIM][NDIM], double cofactorMat[NDIM][NDIM]) {
    double minor[3][3];
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            int subi = 0;
            for (int x = 0; x < NDIM; x++) {
                if (x == i) continue;
                int subj = 0;
                for (int y = 0; y < NDIM; y++) {
                    if (y == j) continue;
                    minor[subi][subj] = mat[x][y];
                    subj++;
                }
                subi++;
            }
            cofactorMat[i][j] = pow(-1, i + j) * determinant3x3(minor);
        }
    }
}

void transpose(double mat[NDIM][NDIM], double transposed[NDIM][NDIM]) {
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            transposed[j][i] = mat[i][j];
        }
    }
}

int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]) {
    double det = determinant4x4(mat);
    double cofactorMat[NDIM][NDIM];
    cofactor(mat, cofactorMat);

    double adjugate[NDIM][NDIM];
    transpose(cofactorMat, adjugate);

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            inverse[i][j] = adjugate[i][j] / det;
        }
    }

    return 1;
}


void check_inverse(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]) {
    double identity[NDIM][NDIM] = {0};
    for (int i = 0; i < NDIM; i++) {
        identity[i][i] = 1.0;
    }

    double product[NDIM][NDIM] = {0};

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                product[i][j] += gcov[i][k] * gcon[k][j];
            }
        }
    }

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            if (fabs(product[i][j] - identity[i][j]) > TOLERANCE) {
                printf("Matrix inversion check failed at element (%d, %d): %f\n", i, j, product[i][j]);
            }
        }
    }
}

void print_matrix(const char* name, double mat[NDIM][NDIM]) {
    printf("%s =\n", name);
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("%.12f ", mat[i][j]);
        }
        printf("\n");
    }
}




inline void gcov(double *X, double gcov[][NDIM])
{
    DLOOP gcov[j][k] = 0.;

    double r, th;
    double sth, cth, s2, rho2;
	r = X[1];
	th = M_PI * X[2];
    


    double sin_theta = sin(th);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta = cos(th);
    double cos_theta2 = cos_theta * cos_theta;

    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2 * M * r + a * a;

    memset(gcov, 0, sizeof(double) * NDIM * NDIM);

    gcov[0][0] = -(1.0 - (2.0 * M * r) / Sigma);   
    gcov[1][1] = Sigma / Delta;   
    gcov[2][2] = Sigma;              
    gcov[3][3] = (r * r + a * a + (2.0 * M * r * a * a * sin_theta2) / Sigma) * sin_theta2; 
    gcov[0][3] = gcov[3][0] = -((2.0 * M * r * a * sin_theta2) / Sigma);

}


inline double calculate_angular_momentum(double v[4], double g[4][4]) {
    double g_phi_phi = g[3][3];
    double g_t_phi = g[0][3];

    double L_z = g_phi_phi * v[3] + g_t_phi * v[0];

    return L_z;
}



inline void gcon(double r, double th, double gcon[][NDIM])
{
    DLOOP gcon[j][k] = 0.;  

    double sth = sin(th);
    double cth = cos(th);
    double a2 = a * a;
    double r2 = r * r;
    double r3 = r2 * r;
    double DD = 1. - 2. / r + a2 / r2;  
    double mu = 1. + a2 * cth * cth / r2;  

    gcon[TT][TT] = -1. - 2. * (1. + a2 / r2) / (r * DD * mu);
    gcon[TT][3]  = -2. * a / (r3 * DD * mu);
    gcon[3][TT]  = gcon[TT][3];
    gcon[1][1]   = DD / mu;
    gcon[2][2]   = 1. / (r2 * mu);
    gcon[3][3]   = (1. - 2. / (r * mu)) / (r2 * sth * sth * DD);
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

