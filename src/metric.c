/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metric.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/22 14:44:57 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/30 15:26:51 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
int i, j , k;
#pragma omp declare simd
void gcov(double *X, double gcov[][NDIM])
{
    DLOOP gcov[j][k] = 0.;
    ALIGNED_32 double dXpdX[NDIM];
    ALIGNED_32 double r, th;
	ALIGNED_32 double R0 = 1.0;
	ALIGNED_32 double hslope = 0.3;
    ALIGNED_32 double sth, cth, s2, rho2;
    Boyer_lindquist_coord(X, &r, &th);
    cth  = cosf(th);
    sth  = sinf(th);
    s2   = sth * sth;
    rho2 = r * r + a * a * cth * cth;

    gcov[TT][TT] = (-1. + 2. * r / rho2);
    gcov[TT][1]  = (2. * r / rho2);
    gcov[TT][3]  = (-2. * a * r * s2 / rho2);
    gcov[1][TT]  = gcov[TT][1];
    gcov[1][1]   = (1. + 2. * r / rho2);
    gcov[1][3]   = (-a * s2 * (1. + 2. * r / rho2));
    gcov[2][2]   = rho2;
    gcov[3][TT]  = gcov[TT][3];
    gcov[3][1]   = gcov[1][3];
    gcov[3][3]   = s2 * (rho2 + a * a * s2 *\
                 (1. + 2. * r / rho2));
                 
    dXpdX[0] = 1. ;
	dXpdX[1] = r - R0;
	dXpdX[2] = M_PI + (1. - hslope) * M_PI * cosf(2. * M_PI * X[2]);
	dXpdX[3] = 1.;
    #pragma omp simd
	for(j=0;j<NDIM;j++) 
	    for(k=0;k<NDIM;k++)
        {
    		gcov[j][k] *= dXpdX[j]*dXpdX[k];
            printf("gcov[%d][%d]: %f\n", j, k, gcov[j][k]);
        } 
}


void Boyer_lindquist_coord(double *X, double *r, double *th)
{
	ALIGNED_32 double  R0 = 1.0, Rin = 2.0, Rout = 10.0, hslope = 0.3;
    *r = expf(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sinf(2. * M_PI * X[2]);
    if (fabs(*th) < SMALL)
    {
        if ((*th) >= 0)
            *th = SMALL;
        if ((*th) < 0)
            *th = -SMALL;
    }
    if (fabs(M_PI - (*th)) < SMALL)
    {
        if ((*th) >= M_PI)
            *th = M_PI + SMALL;
        if ((*th) < M_PI)
            *th = M_PI - SMALL;
    }
    return;
}
void gcon(double r, double th, double gcon[][NDIM])
{
	int j, k;
	ALIGNED_32 double sth, cth, a2, r2, r3, DD, mu;
	DLOOP gcon[j][k] = 0.;

	sth = sinf(th);
	cth = cosf(th);

	a2 = a * a;
	r2 = r * r;
	r3 = r2 * r;
	DD = 1. - 2. / r + a2 / r2;
	mu = 1. + a2 * cth * cth / r2;

	gcon[TT][TT] = -1. - 2. * (1. + a2 / r2) / (r * DD * mu);
	gcon[TT][3] = -2. * a / (r3 * DD * mu);
	gcon[3][TT] = gcon[TT][3];
	gcon[1][1] = DD / mu;
	gcon[2][2] = 1. / (r2 * mu);
	gcon[3][3] = (1. - 2. / (r * mu)) / (r2 * sth * sth * DD);

    for (j = 0; j < NDIM; j++)
        for (k = 0; k < NDIM; k++)
            printf("gcon[%d][%d]: %f\n", j, k, gcon[j][k]);
    printf("\n");
}