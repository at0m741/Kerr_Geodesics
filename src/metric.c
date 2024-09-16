/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metric.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/22 14:44:57 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/09 21:05:26 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
extern double		(*geodesic_points)[5];
extern int			num_points;
int					i, j, k;

void print_matrix(const char* name, double mat[NDIM][NDIM]) {
    printf("%s =\n", name);
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("%.12f ", mat[i][j]);
        }
        printf("\n");
    }
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
            if (fabs(product[i][j] - identity[i][j]) > 1e-10) {
                printf("Matrix inversion check failed at element (%d, %d): %f\n", i, j, product[i][j]);
            }
        }
    }
}

void invert_using_gsl(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]) {
    gsl_matrix *gcov_matrix = gsl_matrix_alloc(NDIM, NDIM);
    gsl_matrix *inverse_matrix = gsl_matrix_alloc(NDIM, NDIM);
    gsl_permutation *perm = gsl_permutation_alloc(NDIM);
    int signum;

    for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
            gsl_matrix_set(gcov_matrix, i, j, gcov[i][j]);
        }
    }

    int status = gsl_linalg_LU_decomp(gcov_matrix, perm, &signum);
    if (status) {
        printf("LU decomposition failed: %d\n", status);
        gsl_matrix_free(gcov_matrix);
        gsl_matrix_free(inverse_matrix);
        gsl_permutation_free(perm);
        return;
    }

    status = gsl_linalg_LU_invert(gcov_matrix, perm, inverse_matrix);
    if (status) {
        printf("Matrix inversion failed: %d\n", status);
        gsl_matrix_free(gcov_matrix);
        gsl_matrix_free(inverse_matrix);
        gsl_permutation_free(perm);
        return;
    }

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            gcon[i][j] = gsl_matrix_get(inverse_matrix, i, j);
        }
    }

    gsl_matrix_free(gcov_matrix);
    gsl_matrix_free(inverse_matrix);
    gsl_permutation_free(perm);
	printf("\n");
	check_inverse(gcov, gcon);
		printf("\n");
		print_matrix("gcov", gcov);
		print_matrix("gcon", gcon);
}



inline void gcov(double *X, double gcov[][NDIM])
{
    DLOOP gcov[j][k] = 0.;

    double r, th;
    double sth, cth, s2, rho2;

    Boyer_lindquist_coord(X, &r, &th);
    
    cth = cos(th);
    sth = sin(th);
    s2  = sth * sth;
    
    rho2 = r * r + a * a * cth * cth;

    gcov[TT][TT] = -(1. - 2. * r / rho2);
    gcov[TT][1]  = 2. * r / rho2;
    gcov[TT][3]  = -2. * a * r * s2 / rho2;
    gcov[1][TT]  = gcov[TT][1];
    gcov[1][1]   = 1. + 2. * r / rho2;
    gcov[1][3]   = -a * s2 * (1. + 2. * r / rho2);
    gcov[2][2]   = rho2;
    gcov[3][TT]  = gcov[TT][3];
    gcov[3][1]   = gcov[1][3];
    gcov[3][3]   = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
}



inline void Boyer_lindquist_coord(double *X, double *r, double *th)
{
    double R0 = 1.0, hslope = 0.3;
    
    *r = exp(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    
    if (fabs(*th) < SMALL)
    {
        *th = (*th >= 0) ? SMALL : -SMALL;
    }
    if (fabs(M_PI - *th) < SMALL)
    {
        *th = (*th >= M_PI) ? M_PI + SMALL : M_PI - SMALL;
    }
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
                printf("Erreur: identity[%d][%d] = %e au lieu de %e\n", i, j, identity[i][j], delta);
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

    printf("Vérification terminée.\n");
}

