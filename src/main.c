/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 15:03:26 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/25 18:32:34 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"


extern double (*geodesic_points)[5];
extern int num_points;

void print_arch()
{
	printf("Architecture: %s\n", ARCH);
	printf("OS: %s\n", Plateform);
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



void Write_Geodesics_to_binary_file(const char *filename, double (*geodesic_points)[5], int num_points)
{
    FILE *file = fopen(filename, "wb");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }

    fwrite(geodesic_points, sizeof(double), num_points * 5, file);
    fclose(file);
}


int main(int argc, char **argv)
{
	#ifdef USE_MPI
		printf("MPI is defined\n");
		MPI_Init(&argc, &argv);
		omp_get_max_threads();
		printf("Number of threads: %d\n", omp_get_max_threads());
	    double start_time = MPI_Wtime();
		double x[4] = {0.0, M_PI / 2, 60 * M_PI, 0.0};
		double v[4] = {-70.0, 1.01, 1.0, 27.0};
		double christoffel_sym[4][4][4] = {0};
		int rank = 0;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		double g[4][4] = {0};
		gcov(x, g);
		if (rank == 0)
		{
			christoffel(g, christoffel_sym);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		geodesic(x, v, max_dt, christoffel_sym, DT, store_geodesic_point);
		if (rank == 0)
		{	
			write_vtk_file("geodesic.vtk");
			Write_Geodesics_to_binary_file("geodesic.bin", geodesic_points, num_points);
		}
		double end_time = MPI_Wtime();
		double elapsed_time = end_time - start_time;

		int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		if (world_rank == 0) {
			printf("Execution time: %f seconds\n", elapsed_time);
		}
	#elif AVX2 
		__m256d dt = _mm256_set1_pd(DT);
		__m256d x[4], v[4], g[4][4], christoffel_avx[4][4][4];
		double x_vals[4] = {1600.0, M_PI / 2, M_PI, 20.0};
		double v_vals[4] = {80.2, 10.0, 12.0, 27.0};
		double g_vals[NDIM][NDIM] = {0};

		for (int i = 0; i < NDIM; i++) {
			x[i] = _mm256_set1_pd(x_vals[i]);
			v[i] = _mm256_set1_pd(v_vals[i]);
		}

		gcov(x_vals, g_vals);
		for (int i = 0; i < NDIM; i++) {
			for (int j = 0; j < NDIM; j++) {
				g[i][j] = _mm256_set1_pd(g_vals[i][j]);
			}
		}
		
		christoffel_AVX(g, christoffel_avx);
		geodesic_AVX(x, v, max_dt, christoffel_avx, dt);
		write_vtk_file("geodesic.vtk");

	#endif
	print_arch();
	printf("num threads: %d\n", omp_get_max_threads());
	if (geodesic_points != NULL)
	{
		free(geodesic_points);
	}
	return 0;
}
