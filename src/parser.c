/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   parser.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/17 15:51:31 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/01 02:39:30 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */



#include "../headers/geodesics.h"

extern double	(*geodesic_points)[5];
extern int		num_points;
int				capacity = 0;

/* 
	* Write the geodesic points to a VTK file
	* The VTK file will contain the geodesic points and the lambda values
	* The lambda values are the affine parameter values
	* The VTK file can be visualized using Paraview
*/

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
        fprintf(file, "%f %f %f\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);

    fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));

    for (int i = 0; i < num_points - 1; ++i)
        fprintf(file, "2 %d %d\n", i, i + 1);

    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS lambda double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int i = 0; i < num_points; ++i)
		fprintf(file, "%f\n", geodesic_points[i][3]);

    printf("Number of points: %d\n", num_points);
    printf("VTK file %s has been written\n", filename);

    fclose(file);
}


void store_geodesic_point_AVX(VEC_TYPE x[4], double lambda) {
    if (num_points >= capacity) 
	{
        capacity = (capacity == 0) ? 1000 : capacity * 2;
        double (*new_geodesic_points)[5];
        if (posix_memalign((void **)&new_geodesic_points, ALIGNMENT, capacity * sizeof(*geodesic_points)) != 0) 
		{
            fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
            exit(EXIT_FAILURE);
        }

        if (num_points > 0) 
		{
            memcpy(new_geodesic_points, geodesic_points, num_points * sizeof(*geodesic_points));
            free(geodesic_points);
        }

        geodesic_points = new_geodesic_points;
    }



    #pragma omp simd aligned(geodesic_points: ALIGNMENT)
    for (int i = 0; i < 4; i+=4) 
    {
        VEC_TYPE r = x[1];
        VEC_TYPE theta = x[2];
        VEC_TYPE phi = x[3];
        _mm_prefetch((const char *)&r, _MM_HINT_NTA);
        _mm_prefetch((const char *)&theta, _MM_HINT_NTA);
        _mm_prefetch((const char *)&phi, _MM_HINT_NTA);
        
        VEC_TYPE sin_theta = _mm256_sin_pd(theta);
        VEC_TYPE cos_theta = _mm256_cos_pd(theta);
        VEC_TYPE sin_phi = _mm256_sin_pd(phi);
        VEC_TYPE cos_phi = _mm256_cos_pd(phi);
        double r_vals[4], theta_vals[4], phi_vals[4];
        VEC_STORE_PD(r_vals, r);
        VEC_STORE_PD(theta_vals, theta);
        VEC_STORE_PD(phi_vals, phi);
        
        geodesic_points[num_points][0] = r_vals[i] * sin_theta[i] * cos_phi[i];
        geodesic_points[num_points][1] = r_vals[i] * sin_theta[i] * sin_phi[i];
        geodesic_points[num_points][2] = r_vals[i] * cos_theta[i];
        geodesic_points[num_points][3] = lambda;
        num_points++;
    }
}

