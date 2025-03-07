#include <Geodesics.h>

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


double GridTensor::partialX_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (i >= 2 && i <= NX - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i+2, j, k).gamma[a][b],
            grid_obj.getCell(i+1, j, k).gamma[a][b],
            grid_obj.getCell(i-1, j, k).gamma[a][b],
            grid_obj.getCell(i-2, j, k).gamma[a][b],
            DX 
        );
    } else if (i >= 1 && i <= NX - 2) {
        return second_order_diff(
            grid_obj.getCell(i+1, j, k).gamma[a][b],
            grid_obj.getCell(i-1, j, k).gamma[a][b],
            DX 
        );
    } else if (i == 0) {
        return (grid_obj.getCell(i+1, j, k).gamma[a][b] - 
                grid_obj.getCell(i, j, k).gamma[a][b]) / DX;
    } else if (i == NX - 1) {
        return (grid_obj.getCell(i, j, k).gamma[a][b] - 
                grid_obj.getCell(i-1, j, k).gamma[a][b]) / DX;
    }
    return 0.0;
}

double GridTensor::partialY_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (j >= 2 && j <= NY - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i, j+2, k).gamma[a][b],
            grid_obj.getCell(i, j+1, k).gamma[a][b],
            grid_obj.getCell(i, j-1, k).gamma[a][b],
            grid_obj.getCell(i, j-2, k).gamma[a][b],
            DY
        );
    } else if (j >= 1 && j <= NY - 2) {
        return second_order_diff(
            grid_obj.getCell(i, j+1, k).gamma[a][b],
            grid_obj.getCell(i, j-1, k).gamma[a][b],
            DY
        );
    } else if (j == 0) {
        return (grid_obj.getCell(i, j+1, k).gamma[a][b] - 
                grid_obj.getCell(i, j, k).gamma[a][b]) / DY;
    } else if (j == NY - 1) {
        return (grid_obj.getCell(i, j, k).gamma[a][b] - 
                grid_obj.getCell(i, j-1, k).gamma[a][b]) / DY;
    }
    return 0.0;
}

double GridTensor::partialZ_gamma(Grid &grid_obj, int i, int j, int k, int a, int b) {
    if (k >= 2 && k <= NZ - 3) {
        return fourth_order_diff(
            grid_obj.getCell(i, j, k+2).gamma[a][b],
            grid_obj.getCell(i, j, k+1).gamma[a][b],
            grid_obj.getCell(i, j, k-1).gamma[a][b],
            grid_obj.getCell(i, j, k-2).gamma[a][b],
            DZ
        );
    } else if (k >= 1 && k <= NZ - 2) {
        return second_order_diff(
            grid_obj.getCell(i, j, k+1).gamma[a][b],
            grid_obj.getCell(i, j, k-1).gamma[a][b],
            DZ
        );
    } else if (k == 0) {
        return (grid_obj.getCell(i, j, k+1).gamma[a][b] - 
                grid_obj.getCell(i, j, k).gamma[a][b]) / DZ;
    } else if (k == NZ - 1) {
        return (grid_obj.getCell(i, j, k).gamma[a][b] - 
                grid_obj.getCell(i, j, k-1).gamma[a][b]) / DZ;
    }
    return 0.0;
}

void GridTensor::compute_partial_christoffel(Grid &grid_obj, int i, int j, int k, int dim, double partialGamma[3][3][3][3], double d) {
    double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];
    double localPartialGamma[3][3][3] = {0.0}; 

    int a = (dim == 0) ? i : ((dim == 1) ? j : k);
    int max_a = (dim == 0) ? NX : ((dim == 1) ? NY : NZ);

    if (a >= 2 && a <= max_a - 3) {
        compute_christoffel_3D(grid_obj, i - 2 * (dim == 0), j - 2 * (dim == 1), k - 2 * (dim == 2), Gmm);
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
        compute_christoffel_3D(grid_obj, i + 2 * (dim == 0), j + 2 * (dim == 1), k + 2 * (dim == 2), Gpp);
    } else if (a >= 1 && a <= max_a - 2) {
        compute_christoffel_3D(grid_obj, i - 1 * (dim == 0), j - 1 * (dim == 1), k - 1 * (dim == 2), Gm);
        compute_christoffel_3D(grid_obj, i + 1 * (dim == 0), j + 1 * (dim == 1), k + 1 * (dim == 2), Gp);
    }

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                if (a >= 2 && a <= max_a - 3) {
                    localPartialGamma[kk][aa][bb] = fourth_order_diff(Gpp[kk][aa][bb],  \
																		Gp[kk][aa][bb], \
																		Gm[kk][aa][bb], \
																		Gmm[kk][aa][bb], d);
                } else if (a >= 1 && a <= max_a - 2) {
                    localPartialGamma[kk][aa][bb] = second_order_diff(Gp[kk][aa][bb], \
																		Gm[kk][aa][bb], d);
                } else {
                    localPartialGamma[kk][aa][bb] = 0.0;
                }
            }
        }
    }

    for (int kk = 0; kk < 3; kk++) {
        for (int aa = 0; aa < 3; aa++) {
            for (int bb = 0; bb < 3; bb++) {
                partialGamma[dim][kk][aa][bb] = localPartialGamma[kk][aa][bb];
            }
        }
    }
}


void Grid::compute_ricci_3D(Grid &grid_obj, int i, int j, int k, double Ricci[3][3]) {
    double Gamma[3][3][3];
    GridTensor gridTensor;
    gridTensor.compute_christoffel_3D(grid_obj, i, j, k, Gamma);

    double T[3][3]; 
    double partialGamma[3][3][3][3] = {};  
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 0, partialGamma, DX);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 1, partialGamma, DY);
    gridTensor.compute_partial_christoffel(grid_obj, i, j, k, 2, partialGamma, DZ);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;
            
            for (int m = 0; m < 3; m++) {
                term1 += partialGamma[m][m][a][b];  
                term2 += partialGamma[a][m][m][b];  
            }
			#pragma omp simd collapse(2)
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    term3 += Gamma[k][a][b] * Gamma[l][k][l]; // Γ^k_ab Γ^l_kl
                    term4 += Gamma[l][a][k] * Gamma[k][b][l]; // Γ^l_ak Γ^k_bl
                }
            }

            Ricci[a][b] = term1 - term2 - term3 + term4;
        }
    }
/*     Cell2D &cell = globalGrid[i][j][k]; */
/* #pragma omp for simd collapse(2) */
/* 	for (int a = 0; a < 3; a++) { */
/* 		for (int b = 0; b < 3; b++) { */
/* 			cell.T[a][b] = (cell.rho + cell.p) * cell.vx * cell.vy + cell.p * cell.gamma[a][b]; */
/* 			 */
/* 		} */
/* 	} */
/*  */
/* 	for (int a = 0; a < 3; a++) { */
/* 		for (int b = 0; b < 3; b++) { */
/* 			Ricci[a][b] += 8 * M_PI * T[a][b]; */
/* 		} */
/* 	} */
/* 	print_matrix_2D("T", cell.T); */
	/* print_matrix_2D("Ricci", Ricci); */
}
