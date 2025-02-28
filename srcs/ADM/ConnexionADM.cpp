#include "GridTensor.h"
#include <Geodesics.h>


void GridTensor::compute_dt_tildeGamma(int i, int j, int k, double dt_tildeGamma[3]) {
	Grid::Cell2D &cell = globalGrid[i][j][k];
	Grid grid_obj;	
    double gammaInv[3][3]; 


    double dtGamma[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            dtGamma[a][b] = cell.dgt[a][b]; 
			gammaInv[a][b] = cell.gamma_inv[a][b];
        }
    }

    double dtGammaInv[3][3];
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            dtGammaInv[a][b] = 0.0;
            for (int m = 0; m < 3; m++) {
                for (int n = 0; n < 3; n++) {
                    dtGammaInv[a][b] -= gammaInv[a][m] * gammaInv[b][n] * dtGamma[m][n];
                }
            }
        }
    }

    double dtChristoffel[3][3][3] = {0.0};  


	for (int a = 0; a < 3; a++) { 
		for (int b = 0; b < 3; b++) { 
			for (int c = 0; c < 3; c++) {
				double sum = 0.0;
				for (int m = 0; m < 3; m++) { 
					sum += gammaInv[a][m] * (
							partialX_gamma(i, j, k, b, c) +
							partialY_gamma(i, j, k, b, c) +
							partialZ_gamma(i, j, k, b, c)
							- partialX_gamma(i, j, k, m, c)  
							- partialY_gamma(i, j, k, m, c)
							- partialZ_gamma(i, j, k, m, c)
							);
				}
				dtChristoffel[a][b][c] = 0.5 * sum;
			}
		}
	}


	for (int i_comp = 0; i_comp < 3; i_comp++) {
		dt_tildeGamma[i_comp] = 0.0;
		for (int j_comp = 0; j_comp < 3; j_comp++) {
			for (int k_comp = 0; k_comp < 3; k_comp++) {
				dt_tildeGamma[i_comp] += gammaInv[j_comp][k_comp] * dtChristoffel[i_comp][j_comp][k_comp];
			}
		}
    }

	/* printf("d_tildeGamma at (%d, %d, %d): [%e, %e, %e]\n", */
	/*        i, j, k, dt_tildeGamma[0], dt_tildeGamma[1], dt_tildeGamma[2]); */
}

void GridTensor::compute_tildeGamma(int i, int j, int k, double tildeGamma[3]) {
    double gammaInv[3][3];
    double christof[3][3][3];

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = globalGrid[i][j][k].gamma_inv[a][b];
        }
    }

    compute_christoffel_3D(i, j, k, christof);

    for (int m = 0; m < 3; m++) {
        tildeGamma[m] = 0.0;
    }

    for (int i_comp = 0; i_comp < 3; i_comp++) { 
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma[i_comp] += gammaInv[j_comp][k_comp] * christof[i_comp][j_comp][k_comp];
            }
        }
    }
	/*  */
	/* printf("tildeGamma at (%d, %d, %d): [%e, %e, %e]\n", */
	/* 		i, j, k, tildeGamma[0], tildeGamma[1], tildeGamma[2]); */
}

void GridTensor::compute_christoffel_3D(int i, int j, int k, double christof[3][3][3]) {
	Matrix matrix_obj;
	double g[3][3];
	double invg[3][3];  
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			invg[a][b] = globalGrid[i][j][k].gamma_inv[a][b];
			g[a][b] = globalGrid[i][j][k].gamma[a][b];
		}
	}

	double dgamma[3][3][3];
	for(int a=0;a<3;a++){
		for(int b=0;b<3;b++){
			dgamma[0][a][b] = partialX_gamma(i,j,k,a,b);
			dgamma[1][a][b] = partialY_gamma(i,j,k,a,b);
			dgamma[2][a][b] = partialZ_gamma(i,j,k,a,b);
		}
	}
	for(int kk=0; kk<3; kk++){
		for(int aa=0; aa<3; aa++){
			for(int bb=0; bb<3; bb++){
				double sum=0.0;
				for(int ll=0; ll<3; ll++){

					double tmp = dgamma[aa][ll][bb] + dgamma[bb][ll][aa] - dgamma[ll][aa][bb];
					sum += invg[kk][ll] * tmp;
				}
				christof[kk][aa][bb] = 0.5 * sum;
			}
		}
	}
	/* printf("Christoffel symbols at (%d, %d, %d):\n", i, j, k); */
	/* for(int a=0; a<3; a++){ */
	/* 	for(int b=0; b<3; b++){ */
	/* 		for(int c=0; c<3; c++){ */
	/* 			printf("Gamma[%d][%d][%d] = %e\n", a, b, c, christof[a][b][c]); */
	/* 		} */
	/* 	} */
	/* } */
}

