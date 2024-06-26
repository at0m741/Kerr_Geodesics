/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   parser.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/17 15:51:31 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/28 18:11:53 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */



#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
#include <fcntl.h>
int capacity = 0;
#define __AVX2__ 1

/* 
*  Write the geodesic points to a VTK file
*  The VTK file will contain the geodesic points and the lambda values
*  The lambda values are the affine parameter values
*  The VTK file can be visualized using Paraview
*/
#ifdef __AVX2__
void write_vtk_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Geodesic Points\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    fprintf(file, "POINTS %d double\n", num_points);
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%f %f %f\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);
    }

    fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));
    for (int i = 0; i < num_points - 1; ++i)
    {
        fprintf(file, "2 %d %d\n", i, i + 1);
    }

    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS lambda double\n");
    fprintf(file, "LOOKUP_TABLE default\n");
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%f\n", geodesic_points[i][3]);
    }
    printf("Number of points: %d\n", num_points);
    printf("VTK file %s has been written\n", filename);

    fclose(file);
}


inline void store_geodesic_point_AVX(__m256d x[4], double lambda) {
    if (num_points >= capacity) {
        capacity = (capacity == 0) ? 1000 : capacity * 2;
        double (*new_geodesic_points)[5];
        if (posix_memalign((void **)&new_geodesic_points, ALIGNMENT, capacity * sizeof(*geodesic_points)) != 0) {
            fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
            exit(EXIT_FAILURE);
        }

        if (num_points > 0) {
            memcpy(new_geodesic_points, geodesic_points, num_points * sizeof(*geodesic_points));
            free(geodesic_points);
        }

        geodesic_points = new_geodesic_points;
    }



    #pragma omp for simd aligned(geodesic_points: ALIGNMENT)
    for (int i = 0; i < 4; i+=4) 
    {
        __m256d r = x[1];
        __m256d theta = x[2];
        __m256d phi = x[3];
        _mm_prefetch((const char *)&r, _MM_HINT_NTA);
        _mm_prefetch((const char *)&theta, _MM_HINT_NTA);
        _mm_prefetch((const char *)&phi, _MM_HINT_NTA);
        
        __m256d sin_theta = _mm256_sin_pd(theta);
        __m256d cos_theta = _mm256_cos_pd(theta);
        __m256d sin_phi = _mm256_sin_pd(phi);
        __m256d cos_phi = _mm256_cos_pd(phi);
        double r_vals[4], theta_vals[4], phi_vals[4];
        _mm256_storeu_pd(r_vals, r);
        _mm256_storeu_pd(theta_vals, theta);
        _mm256_storeu_pd(phi_vals, phi);
        
        geodesic_points[num_points][0] = r_vals[i] * sin_theta[i] * cos_phi[i];
        geodesic_points[num_points][1] = r_vals[i] * sin_theta[i] * sin_phi[i];
        geodesic_points[num_points][2] = r_vals[i] * cos_theta[i];
        geodesic_points[num_points][3] = lambda;
        num_points++;
    }
}

#else


#pragma omp declare simd
void write_vtk_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Geodesic Points\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    fprintf(file, "POINTS %d double\n", num_points);
    #pragma omp for simd aligned(geodesic_points: ALIGNMENT)
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%f %f %f\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);
    }

    fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));
	#pragma omp simd aligned(geodesic_points: ALIGNMENT)
    for (int i = 0; i < num_points - 1; ++i)
    {
        fprintf(file, "2 %d %d\n", i, i + 1);
    }

    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS lambda double\n");
    fprintf(file, "LOOKUP_TABLE default\n");
	#pragma omp simd aligned(geodesic_points: ALIGNMENT)
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%f\n", geodesic_points[i][3]);
    }
	printf("Number of points: %d\n", num_points);
	printf("VTK file %s has been written\n", filename);

    fclose(file);
}


#pragma omp declare simd
void store_geodesic_point(double x[4], double lambda)
{
    if (num_points >= capacity)
    {
        capacity = (capacity == 0) ? 1000 : capacity * 2;
        double (*new_geodesic_points)[5] = aligned_alloc(ALIGNMENT, capacity * sizeof(*geodesic_points));
        if (new_geodesic_points == NULL)
        {
            fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
            exit(1);
        }
        
        if (num_points > 0)
        {
            memcpy(new_geodesic_points, geodesic_points, num_points * sizeof(*geodesic_points));
            free(geodesic_points);
        }
        
        geodesic_points = new_geodesic_points;
    }

    ldouble_a r = x[1];
    ldouble_a theta = x[2];
    ldouble_a phi = x[3];
    ldouble_a sin_theta, cos_theta, sin_phi, cos_phi;
	ldouble_a dt = x[4];

    sincos(theta, &sin_theta, &cos_theta);
    sincos(phi, &sin_phi, &cos_phi);

    geodesic_points[num_points][0] = r * sin_theta * cos_phi;
    geodesic_points[num_points][1] = r * sin_theta * sin_phi;
    geodesic_points[num_points][2] = r * cos_theta;
    geodesic_points[num_points][3] = lambda;
	geodesic_points[num_points][4] = dt;


    num_points++;
}

#endif

