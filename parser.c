#include "geodesics.h"


void write_vtk_file(const char *filename)
{
	long double (*geodesic_points)[4] = NULL;
	int num_points = 0;
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
        fprintf(file, "%Lf %Lf %Lf\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);
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
        fprintf(file, "%Lf\n", geodesic_points[i][3]);
    }
    fclose(file);
	free(geodesic_points);
}


void store_geodesic_point(long double x[4], long double lambda)
{
	long double (*geodesic_points)[4] = NULL;
	int num_points = 0;
    long double r = x[1];
    long double theta = x[2];
    long double phi = x[3];
    long double sin_theta = sin(theta);
    long double cos_theta = cos(theta);
    long double sin_phi = sin(phi);
    long double cos_phi = cos(phi);

    geodesic_points = realloc(geodesic_points, (num_points + 1000) * sizeof(*geodesic_points));
    if (geodesic_points == NULL)
    {
        fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
        exit(1);
    }

    geodesic_points[num_points][0] = r * sin_theta * cos_phi;
    geodesic_points[num_points][1] = r * sin_theta * sin_phi;
    geodesic_points[num_points][2] = r * cos_theta;
	geodesic_points[num_points][3] = lambda; 

    num_points++;
}