#include "GridTensor.h"
#include <Geodesics.h>


void GridTensor::compute_dt_tildeGamma(int i, int j, int k, double dt_tildeGamma[3]) {
	Grid grid_obj;
	Grid::Cell2D &cell = grid_obj.getCell(i, j, k);
    
    double gammaInv[3][3], KLocal[3][3], Atilde[3][3], dBeta[3][3];
    double alpha = cell.alpha;
    double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };

    int iP = std::min(i + 1, NX - 1);
    int iM = std::max(i - 1, 0);
    int jP = std::min(j + 1, NY - 1);
    int jM = std::max(j - 1, 0);
    int kP = std::min(k + 1, NZ - 1);
    int kM = std::max(k - 1, 0);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = cell.gamma_inv[a][b];
            KLocal[a][b] = cell.K[a][b];
        }
    }

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * KLocal[a][b];
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Atilde[a][b] = KLocal[a][b] - (1.0 / 3.0) * Ktrace * cell.gamma[a][b];
        }
    }

    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            dBeta[m][n] = (grid_obj.getCell(iP, jP, kP).beta[n] - grid_obj.getCell(iM, jM, kM).beta[n]) / (2.0 * DX);
        }
    }

    for (int i_comp = 0; i_comp < 3; i_comp++) {
        double div_Atilde = 0.0;
        double tildeGamma_Atilde = 0.0;
        double beta_term = 0.0;

        for (int j_comp = 0; j_comp < 3; j_comp++) {
            div_Atilde += (grid_obj.getCell(iP, jP, kP).K[i_comp][j_comp] -
                           grid_obj.getCell(iM, jM, kM).K[i_comp][j_comp]) / (2.0 * DX);
        }

        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma_Atilde += gammaInv[j_comp][k_comp] * Atilde[j_comp][k_comp];
            }
        }

        for (int j_comp = 0; j_comp < 3; j_comp++) {
            beta_term += beta[j_comp] * (grid_obj.getCell(iP, jP, kP).beta[i_comp] -
                                          grid_obj.getCell(iM, jM, kM).beta[i_comp]) / (2.0 * DX);
        }

        dt_tildeGamma[i_comp] = -2.0 * alpha * div_Atilde + 2.0 * alpha * tildeGamma_Atilde + beta_term;
    }
}

void GridTensor::compute_tildeGamma(int i, int j, int k, double tildeGamma[3]) {
    double gammaInv[3][3];
    double christof[3][3][3];
	Grid grid_obj;
    Grid::Cell2D &cell = grid_obj.getCell(i, j, k);

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaInv[a][b] = cell.gamma_inv[a][b];
            for (int c = 0; c < 3; c++) {
                christof[a][b][c] = cell.Christoffel[a][b][c];
            }
        }
    }

    std::fill_n(tildeGamma, 3, 0.0);

    for (int i_comp = 0; i_comp < 3; i_comp++) { 
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            for (int k_comp = 0; k_comp < 3; k_comp++) {
                tildeGamma[i_comp] += gammaInv[j_comp][k_comp] * christof[i_comp][j_comp][k_comp];
            }
        }
    }
}

void GridTensor::compute_christoffel_3D(Grid &grid_obj, int i, int j, int k, double christof[3][3][3]) {
	Matrix matrix_obj;
	double g[3][3];
	double invg[3][3]; 
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			invg[a][b] = grid_obj.getCell(i, j, k).gamma_inv[a][b];
			g[a][b] = grid_obj.getCell(i, j, k).gamma[a][b];
			for (int c = 0; c < 3; c++) {
				christof[a][b][c] = grid_obj.getCell(i, j, k).Christoffel[a][b][c];
			}
		}
	}

	double dgamma[3][3][3];
	for(int a=0;a<3;a++){
		for(int b=0;b<3;b++){
			dgamma[0][a][b] = partialX_gamma(grid_obj, i, j, k, a, b);
			dgamma[1][a][b] = partialY_gamma(grid_obj, i, j, k, a, b);
			dgamma[2][a][b] = partialZ_gamma(grid_obj, i, j, k, a, b);
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
				grid_obj.getCell(i, j, k).Christoffel[kk][aa][bb] = 0.5 * sum;
				
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

