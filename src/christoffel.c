#include "../headers/geodesics.h"

extern double (*geodesic_points)[5];
extern int num_points;

#pragma omp declare simd
void christoffel(double g[4][4], double christoffel[4][4][4])
{
    #ifdef USE_MPI
        printf("MPI is defined for Christoffel\n");
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#elif _OPENMP
		printf("OpenMP is supported and used\n");
	#else
		printf("OpenMP is not supported\n");
	#endif

    double (*g_aligned)[4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4]));
    double (*christoffel_aligned)[4][4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4][4]));
    memcpy(g_aligned, g, sizeof(double[4][4]));
    memcpy(christoffel_aligned, christoffel, sizeof(double[4][4][4]));


    #pragma omp parallel for collapse(3)
    #pragma vector aligned
    for (int mu = 0; mu < 4; mu++) {
        for (int beta = 0; beta < 4; beta++) {
            for (int nu = 0; nu < 4; nu++) {
                double sum = 0;
                #pragma omp simd reduction(+:sum) aligned(g_aligned:ALIGNMENT)
                #pragma vector aligned
                for (int sigma = 0; sigma < 4; sigma++) {
                    sum += 0.5 * (g_aligned[mu][sigma] * (g_aligned[sigma][beta] + g_aligned[beta][sigma] - g_aligned[beta][nu]));
                }
                christoffel_aligned[mu][beta][nu] = sum;
            }
        }
    }

    memcpy(christoffel, christoffel_aligned, sizeof(double[4][4][4]));
    free(g_aligned);
    free(christoffel_aligned);
}


#pragma omp declare simd
void riemann(double g[4][4], double christoffel[4][4][4], double riemann[4][4][4][4])
{
	#ifdef _OPENMP
		printf("OpenMP is supported and used\n");
	#else
		printf("OpenMP is not supported\n");
	#endif
    #pragma omp parallel for collapse(4)
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            for (int k = 1; k < 4; k++) {
                for (int l = 1; l < 4; l++) {
                    double sum = 0;
                    #pragma omp simd reduction(+:sum) aligned(christoffel,g:ALIGNMENT)
                    for (int m = 0; m < 4; m++) 
					{
                        sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] * christoffel[m][j][k] - christoffel[k][j][m] * christoffel[m][i][k]);
                        sum += (1 / (2 * g[0][0])) * (christoffel[k][i][m] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - christoffel[k][j][m] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                        sum += (1 / (2 * g[0][0])) * (christoffel[m][i][k] * (g[m][j] * g[k][k] - g[m][k] * g[j][k]) - christoffel[m][j][k] * (g[m][i] * g[k][k] - g[m][k] * g[i][k]));
                    }
                    riemann[i][j][k][l] = sum;
                    printf("Riemann[%d][%d][%d][%d] = %f\n", i, j, k, l, riemann[i][j][k][l]);
                }
            }
        }
    }
}
