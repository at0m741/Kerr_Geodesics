/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   metric.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/22 14:44:57 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/22 14:46:09 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
int j, k;
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)

void gcov(double *X, double gcov[][NDIM])
{
    DLOOP gcov[j][k] = 0.;
    double dXpdX[NDIM];

    double r, th;
	double R0 = 1.0;
	double hslope = 0.3;
    Boyer_lindquist_coord(X, &r, &th);
    double sth, cth, s2, rho2;
    cth = cos(th);
    sth = sin(th);
    s2 = sth * sth;
    rho2 = r * r + a * a * cth * cth;

    gcov[TT][TT] = (-1. + 2. * r / rho2);
    gcov[TT][1] = (2. * r / rho2);
    gcov[TT][3] = (-2. * a * r * s2 / rho2);
    gcov[1][TT] = gcov[TT][1];
    gcov[1][1] = (1. + 2. * r / rho2);
    gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2));
    gcov[2][2] = rho2;
    gcov[3][TT] = gcov[TT][3];
    gcov[3][1] = gcov[1][3];
    gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));

}


void Boyer_lindquist_coord(double *X, double *r, double *th)
{
	double  R0 = 1.0, Rin = 2.0, Rout = 10.0, hslope = 0.3;
    *r = exp(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
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