#include <Geodesics.h>

void Grid::compute_gauge_derivatives(int i, int j, int k, double &d_alpha_dt, double d_beta_dt[3]) {
    Grid::Cell2D &cell = globalGrid[i][j][k];
    double gammaLocal[3][3], KLocal[3][3];
	GridTensor gridTensor;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.gamma[a][b];
            KLocal[a][b]     = cell.K[a][b];
        }
    }

    double gammaInv[3][3];
    bool ok = invert_3x3(gammaLocal, gammaInv);
    if (!ok) {
        d_alpha_dt = 0.0;
        for (int m = 0; m < 3; m++) d_beta_dt[m] = 0.0;
        return;
	}

    double Ktrace = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * KLocal[a][b];
		}
    }

    double lambda = 1.0 / (1.0 + 2.0 * Ktrace * Ktrace);
    d_alpha_dt = -2.0 * cell.alpha * Ktrace * lambda;

    double eta = 2.0 / (1.0 + std::fabs(Ktrace));
    double d_Gamma_dt[3] = {0.0, 0.0, 0.0}; 

    double tildeGamma[3];
    gridTensor.compute_tildeGamma(i, j, k, tildeGamma); 
    gridTensor.compute_dt_tildeGamma(i, j, k, d_Gamma_dt);

    for (int m = 0; m < 3; m++) {
        d_beta_dt[m] = 3.0 / 4.0 * d_Gamma_dt[m] - eta * cell.beta[m];
    }
}

std::vector<std::vector<std::vector<double>>> hamiltonianGrid;
void Grid::initialize_grid() {
    globalGrid.resize(NX, std::vector<std::vector<Cell2D>>(NY, std::vector<Cell2D>(NZ)));
    hamiltonianGrid.resize(NX, std::vector<std::vector<double>>(NY, std::vector<double>(NZ, 0.0)));
}


void Grid::compute_constraints(Grid &grid_obj, int i, int j, int k, double &hamiltonian, double momentum[3]) {
    double Ricci[3][3];
    compute_ricci_3D(grid_obj, i, j, k, Ricci);
    
    double gammaLocal[3][3], gammaInv[3][3];
    Cell2D &cell = globalGrid[i][j][k];

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            gammaLocal[a][b] = cell.gamma[a][b];
        }
    }
    invert_3x3(gammaLocal, gammaInv);

    double R = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            R += gammaInv[a][b] * Ricci[a][b];
        }
    }

    double Ktrace = 0.0;
    double KK = 0.0;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            Ktrace += gammaInv[a][b] * cell.K[a][b];
            for (int c = 0; c < 3; c++) {
                KK += cell.K[a][b] * gammaInv[a][c] * cell.K[c][b];
            }
        }
    }
    
    hamiltonian = R + Ktrace * Ktrace - KK;
    hamiltonianGrid[i][j][k] = hamiltonian;
	/* if (hamiltonian > 5.0) { */
	/* 	printf("hamiltonian = %f, at i = %d, j = %d, k = %d\n", hamiltonian, i, j, k); */
	/* } */
    for (int i_comp = 0; i_comp < 3; i_comp++) {
        momentum[i_comp] = 0.0;
        for (int j_comp = 0; j_comp < 3; j_comp++) {
            double dKdx = (globalGrid[i+1][j][k].K[i_comp][j_comp] - globalGrid[i-1][j][k].K[i_comp][j_comp]) / (2.0 * DX);
            double dKdy = (globalGrid[i][j+1][k].K[i_comp][j_comp] - globalGrid[i][j-1][k].K[i_comp][j_comp]) / (2.0 * DY);
            double dKdz = (globalGrid[i][j][k+1].K[i_comp][j_comp] - globalGrid[i][j][k-1].K[i_comp][j_comp]) / (2.0 * DZ);
            momentum[i_comp] += (dKdx + dKdy + dKdz);
			/* if (momentum[i_comp] > 5.0) { */
			/* 	printf("momentum[%d] = %f, at i = %d, j = %d, k = %d\n", i_comp, momentum[i_comp], i, j, k); */
			/*          */
			/* } */
		}
    }
}
