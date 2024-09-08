/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 15:03:26 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/08 23:17:20 by at0m             ###   ########.fr       */
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


int main(int argc, char **argv)
{
	VEC_TYPE	dt = VEC_SET_PD(DT);
	VEC_TYPE	x[4], v[4], g[4][4], christoffel_avx[4][4][4], g_contr[4][4];
	double		x_vals[4] = {20.0, M_PI * 2.3, 1.1, 0.0};
	double		v_vals[4] = {10.2, 1.0, -1.0, 27.0};
	double		g_vals[NDIM][NDIM] = {0};
	double		g_con[NDIM][NDIM] = {0};
	
	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
	{
		x[i] = VEC_SET_PD(x_vals[i]);
		v[i] = VEC_SET_PD(v_vals[i]);
	}

	gcov(x_vals, g_vals);
	gcon(x_vals[1], x_vals[2], g_con);
	printf("Compute Christoffel symbols\n");
	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			g[i][j] = VEC_SET_PD(g_vals[i][j]);
	
	#pragma omp for simd
	for (int i = 0; i < NDIM; i++) 
		for (int j = 0; j < NDIM; j++)
			g_contr[i][j] = VEC_SET_PD(g_con[i][j]);

	christoffel_AVX(g, christoffel_avx, g_contr);
	printf("Compute geodesics equations using Runge Kutta..\n");
	geodesic_AVX(x, v, max_dt, christoffel_avx, dt);
	printf("writing to file..\n");
	write_vtk_file("geodesic.vtk");

	if (geodesic_points != NULL)
		free(geodesic_points);

	return 0;
}
