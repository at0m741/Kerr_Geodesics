#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;
int capacity = 0;


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
    // store the points of the geodesic
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
    //add RS 

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

void write_obj_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }
    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "v %f %f %f\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);
    }

    // set the number of points in each line
    for (int i = 1; i < num_points; ++i)
    {
        fprintf(file, "f %d %d %d\n", i, i + 1, i + 2);
    }
    fprintf(file, "\n");
    printf("Number of points: %d\n", num_points);
    printf("OBJ file %s has been written\n", filename);

    fclose(file);
}

// void write_hdf5(const char *fileneame)
// {
//     hid_t file_id, dataset_id, dataspace_id;
//     hsize_t dims[2] = {num_points, 5};
//     herr_t status;
//     file_id = H5Fcreate(fileneame, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//     dataspace_id = H5Screate_simple(2, dims, NULL);
//     dataset_id = H5Dcreate(file_id, "geodesic_points", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, geodesic_points);
//     H5Dclose(dataset_id);
//     H5Sclose(dataspace_id);
//     H5Fclose(file_id);
//     printf("HDF5 file %s has been written\n", fileneame);
// }

void write_binary(const char *filename)
{
    FILE *file = fopen(filename, "wb");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }
    fwrite(&num_points, sizeof(int), 1, file);
    fwrite(geodesic_points, sizeof(double), num_points * 5, file);
    printf("Binary file %s has been written\n", filename);
    fclose(file);
}