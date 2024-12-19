/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 15:03:26 by ltouzali          #+#    #+#             */
/*   Updated: 2024/12/12 02:43:19 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

extern double	(*geodesic_points)[5];
extern int		num_points;

void print_arch()
{
	printf("Architecture: %s\n", ARCH);
	if (ALIGNMENT == 64)
	{
		printf("AVX-512 is supported\n");
		printf("AVX2 is supported\n");
		printf("SSE is supported\n");
		printf("SSE2 is supported\n");
	}
	else if (ALIGNMENT == 32)
	{
		printf("AVX2 is supported\n");
		printf("SSE2 is supported\n");
		printf("SSE is supported\n");
	}
	else
	{
		printf("SSE is supported\n");
	}
	printf("Alignment: %d bits\n", ALIGNMENT);
}


void initialize_light_geodesic(double v[4], double g[4][4]) {
    v[0] = 1.0;
    v[1] = 0.5;
    v[3] = 0.1;

    double term = g[0][0] * v[0] * v[0]
                + 2 * g[0][1] * v[0] * v[1]
                + g[1][1] * v[1] * v[1]
                + 2 * g[0][3] * v[0] * v[3]
                + 2 * g[1][3] * v[1] * v[3]
                + g[3][3] * v[3] * v[3];

    v[2] = sqrt(fabs(-term / g[2][2]));
	v[2] = 0.0;
    printf("Initialized velocity for light: v[0] = %f, v[1] = %f, v[2] = %f, v[3] = %f\n",
           v[0], v[1], v[2], v[3]);
}

int main(int argc, char **argv)
{
	VEC_TYPE	dt = VEC_SET_PD(DT);
	VEC_TYPE	x[4], v[4], g[4][4], christoffel_avx[4][4][4], g_contr[4][4];
	double		x_vals[4] = {0.4, M_PI / 2.0, M_PI, 1.0};; // {r, theta, phi, t}
	double		v_vals[4] = {1.0, 0.4, 0.0, 1.0};
	double		g_vals[NDIM][NDIM] = {0};
	double		g_con[NDIM][NDIM] = {0};
	double		h = 1e-5;  
    double		gcov_output_forward[NDIM][NDIM], gcov_output_backward[NDIM][NDIM];
	double		dgcov_mu[NDIM][NDIM][NDIM];

	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			g[i][j] = VEC_SET_PD(g_vals[i][j]);
	
	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			g_contr[i][j] = VEC_SET_PD(g_con[i][j]);

	christoffel_symbols(x_vals, h, dgcov_mu);
	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
	{
		x[i] = VEC_SET_PD(x_vals[i]);
		v[i] = VEC_SET_PD(v_vals[i]);
	}


	#pragma omp for simd	
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			for (int k = 0; k < NDIM; k++)
				christoffel_avx[i][j][k] = VEC_SET_PD(dgcov_mu[i][j][k]);
	gcov(x_vals, g_vals);

	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			g[i][j] = VEC_SET_PD(g_vals[i][j]);

	printf("Compute geodesics equations using Runge Kutta..\n");
	geodesic_AVX(x, v, max_dt, christoffel_avx, dt, g); 
	printf("writing to file..\n");
	write_vtk_file("geodesic.vtk");

	if (geodesic_points != NULL)
		free(geodesic_points);

	return 0;
}
