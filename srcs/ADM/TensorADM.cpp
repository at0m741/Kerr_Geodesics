#include <Geodesics.h>

static void print_christoffel(const char* name, const double Gamma[3][3][3])
{
    printf("\n%s:\n", name);
    for(int k = 0; k < 3; k++)
    {
        printf("  Gamma^%d_{ij}:\n", k);
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                printf("%12.6f ", Gamma[k][i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

static void print_partialGamma(const char* name, const double partialGamma[3][3][3][3])
{
    printf("\n%s:\n", name);
    for(int dim = 0; dim < 3; dim++)
    {
        printf("  Dimension = %d:\n", dim);
        for(int k = 0; k < 3; k++)
        {
            printf("    k = %d:\n", k);
            for(int a = 0; a < 3; a++)
            {
                for(int b = 0; b < 3; b++)
                {
                    printf("%12.6f ", partialGamma[dim][k][a][b]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}

static void print_matrix_2D(const char* name, const double mat[3][3])
{
    printf("\n%s (3x3):\n", name);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            printf("%12.6f ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

static void print_3Darray(const char* name, const double arr[3][3][3])
{
    printf("\n%s (3x3x3):\n", name);
    for(int i = 0; i < 3; i++)
    {
        printf("  i=%d:\n", i);
        for(int j = 0; j < 3; j++)
        {
            for(int k = 0; k < 3; k++)
            {
                printf("%12.6f ", arr[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
}



void Grid::compute_ricci_3D(int i, int j, int k, double Ricci[3][3])
{
    double Gamma[3][3][3];
    compute_christoffel_3D(i, j, k, Gamma);
    double partialGamma[3][3][3][3];
    memset(partialGamma, 0, sizeof(partialGamma));

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (i >= 2 && i <= NX - 3) {
            compute_christoffel_3D(i-2, j, k, Gmm);
            compute_christoffel_3D(i-1, j, k, Gm );
            compute_christoffel_3D(i+1, j, k, Gp );
            compute_christoffel_3D(i+2, j, k, Gpp);
        }
        else if (i >= 1 && i <= NX - 2) {
            compute_christoffel_3D(i-1, j, k, Gm );
            compute_christoffel_3D(i+1, j, k, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (i >= 2 && i <= NX - 3) {
						partialGamma[0][kk][aa][bb] = fourth_order_diff( 
							Gpp[kk][aa][bb],
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							Gmm[kk][aa][bb],
							DX
						);
					} else if (i >= 1 && i <= NX - 2) {
						partialGamma[0][kk][aa][bb] = second_order_diff(
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							DX 
						);
					} else if (i == 0) {
						partialGamma[0][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DX;
					} else if (i == NX - 1) {
						partialGamma[0][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DX;
                    } else {
                        partialGamma[0][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (j >= 2 && j <= NY - 3) {
            compute_christoffel_3D(i, j-2, k, Gmm);
            compute_christoffel_3D(i, j-1, k, Gm );
            compute_christoffel_3D(i, j+1, k, Gp );
            compute_christoffel_3D(i, j+2, k, Gpp);
        }
        else if (j >= 1 && j <= NY - 2) {
            compute_christoffel_3D(i, j-1, k, Gm );
            compute_christoffel_3D(i, j+1, k, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (j >= 2 && j <= NY - 3) {
						partialGamma[1][kk][aa][bb] = fourth_order_diff(
							Gpp[kk][aa][bb],
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							Gmm[kk][aa][bb],
							DY
						);
					}else if (j >= 1 && j <= NY - 2) {
						partialGamma[1][kk][aa][bb] = second_order_diff(
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							DY
						);
					} else if (j == 0) {
						partialGamma[1][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DY;
					} else if (j == NY - 1) {
						partialGamma[1][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DY;
                    } else {
                        partialGamma[1][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (k >= 2 && k <= NZ - 3) {
            compute_christoffel_3D(i, j, k-2, Gmm);
            compute_christoffel_3D(i, j, k-1, Gm );
            compute_christoffel_3D(i, j, k+1, Gp );
            compute_christoffel_3D(i, j, k+2, Gpp);
        }
        else if (k >= 1 && k <= NZ - 2) {
            compute_christoffel_3D(i, j, k-1, Gm );
            compute_christoffel_3D(i, j, k+1, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (k >= 2 && k <= NZ - 3) {
						partialGamma[2][kk][aa][bb] = fourth_order_diff( 
							Gpp[kk][aa][bb],
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							Gmm[kk][aa][bb],
							DZ
						);
					} else if (k >= 1 && k <= NZ - 2) {
						partialGamma[2][kk][aa][bb] = second_order_diff(
							Gp[kk][aa][bb],
							Gm[kk][aa][bb],
							DZ
						);
					} else if (k == 0) {
						partialGamma[2][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DZ;
					} else if (k == NZ - 1) {
						partialGamma[2][kk][aa][bb] = (Gp[kk][aa][bb] - Gm[kk][aa][bb]) / DZ;
                    } else {
                        partialGamma[2][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;

            for (int m = 0; m < 3; m++) {
                term1 += partialGamma[m][m][a][b];
                term2 += partialGamma[a][m][m][b];
            }

            for (int kk = 0; kk < 3; kk++) {
                for (int ll = 0; ll < 3; ll++) {
                    term3 += Gamma[kk][a][b] * Gamma[ll][kk][ll];
                    term4 += Gamma[ll][a][kk] * Gamma[kk][b][ll];
                }
            }

            Ricci[a][b] = term1 - term2 + term3 - term4;
        }
    }

    print_matrix_2D("Ricci(a,b)", Ricci);
}
