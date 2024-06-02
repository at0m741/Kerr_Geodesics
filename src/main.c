#include "../headers/geodesics.h"
int j, k;
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)

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
	#else
		printf("MPI is not defined\n");
	#endif

	print_arch();
	printf("num threads: %d\n", omp_get_max_threads());
	if (geodesic_points != NULL)
	{
		free(geodesic_points);
	}
	return 0;
}
