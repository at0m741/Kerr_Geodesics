/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   parser.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/17 15:51:31 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/09 17:44:07 by at0m             ###   ########.fr       */
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

    if (num_points > 1) {
        fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));
        for (int i = 0; i < num_points - 1; ++i)
            fprintf(file, "2 %d %d\n", i, i + 1);
    }

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

    __attribute__((aligned(32))) double r_vals[4], theta_vals[4], phi_vals[4];
    for (int i = 0; i < 4; ++i) 
    {
        VEC_TYPE r = x[1];
        VEC_TYPE theta = x[2];
        VEC_TYPE phi = x[3];

        VEC_STORE_PD(r_vals, r);
        VEC_STORE_PD(theta_vals, theta);
        VEC_STORE_PD(phi_vals, phi);
        
        double sin_theta_vals[4], cos_theta_vals[4], sin_phi_vals[4], cos_phi_vals[4];
        for (int j = 0; j < 4; ++j) {
            sin_theta_vals[j] = sin(theta_vals[j]);
            cos_theta_vals[j] = cos(theta_vals[j]);
            sin_phi_vals[j] = sin(phi_vals[j]);
            cos_phi_vals[j] = cos(phi_vals[j]);
        }

        for (int j = 0; j < 4; ++j) {
            geodesic_points[num_points][0] = r_vals[j] * sin_theta_vals[j] * cos_phi_vals[j];
            geodesic_points[num_points][1] = r_vals[j] * sin_theta_vals[j] * sin_phi_vals[j];
            geodesic_points[num_points][2] = r_vals[j] * cos_theta_vals[j];
            geodesic_points[num_points][3] = lambda;
        }
    }
    num_points++;
}
