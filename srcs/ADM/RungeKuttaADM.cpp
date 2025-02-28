#include <Geodesics.h>


void apply_boundary_conditions()
{
    for(int j=0; j<NY; j++){
        for(int k=0; k<NZ; k++){
            globalGrid[0][j][k] = globalGrid[1][j][k];        
            globalGrid[NX-1][j][k] = globalGrid[NX-2][j][k]; 
        }
    }

    for(int i=0; i<NX; i++){
        for(int k=0; k<NZ; k++){
            globalGrid[i][0][k] = globalGrid[i][1][k];        
            globalGrid[i][NY-1][k] = globalGrid[i][NY-2][k]; 
        }
    }

    for(int i=0; i<NX; i++){
        for(int j=0; j<NY; j++){
            globalGrid[i][j][0] = globalGrid[i][j][1];
            globalGrid[i][j][NZ-1] = globalGrid[i][j][NZ-2];
        }
    }
}

void Grid::copyInitialState(Cell2D &cell) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gamma0[a][b] = cell.gamma[a][b];
            cell.K0[a][b]     = cell.K[a][b];
        }
    }
    cell.alpha0 = cell.alpha;
    for (int m = 0; m < 3; m++) {
        cell.beta0[m] = cell.beta[m];
    }
}

void Grid::updateIntermediateState(Cell2D &cell, double dtCoeff, int stageIndex) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gamma[a][b] = cell.gamma0[a][b] + dtCoeff * cell.gammaStage[stageIndex][a][b];
            cell.K[a][b]     = cell.K0[a][b]     + dtCoeff * cell.KStage[stageIndex][a][b];
        }
    }
    cell.alpha = cell.alpha0 + dtCoeff * cell.alphaStage[stageIndex];
    for (int m = 0; m < 3; m++) {
        cell.beta[m] = cell.beta0[m] + dtCoeff * cell.betaStage[stageIndex][m];
    }
}

void Grid::storeStage(Cell2D &cell, int stage, double d_alpha_dt, double d_beta_dt[3]) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gammaStage[stage][a][b] = cell.dgt[a][b];
            cell.KStage[stage][a][b]     = cell.dKt[a][b];
        }
    }
    cell.alphaStage[stage] = d_alpha_dt;
    for (int m = 0; m < 3; m++) {
        cell.betaStage[stage][m] = d_beta_dt[m];
    }
}

void Grid::combineStages(Cell2D &cell, double dt) {
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            cell.gamma[a][b] = cell.gamma0[a][b] +
                (dt / 6.0) * (cell.gammaStage[0][a][b] +
                              2.0 * cell.gammaStage[1][a][b] +
                              2.0 * cell.gammaStage[2][a][b] +
                              cell.gammaStage[3][a][b]);
            cell.K[a][b] = cell.K0[a][b] +
                (dt / 6.0) * (cell.KStage[0][a][b] +
                              2.0 * cell.KStage[1][a][b] +
                              2.0 * cell.KStage[2][a][b] +
                              cell.KStage[3][a][b]);
        }
    }
    cell.alpha = cell.alpha0 +
        (dt / 6.0) * (cell.alphaStage[0] +
                      2.0 * cell.alphaStage[1] +
                      2.0 * cell.alphaStage[2] +
                      cell.alphaStage[3]);
    for (int m = 0; m < 3; m++) {
        cell.beta[m] = cell.beta0[m] +
            (dt / 6.0) * (cell.betaStage[0][m] +
                          2.0 * cell.betaStage[1][m] +
                          2.0 * cell.betaStage[2][m] +
                          cell.betaStage[3][m]);
    }
}

void Grid::evolve(double dtinitital, int nSteps) {
    initialize_grid(); 
	GridTensor gridTensor;
	double CFL = 0.5;
	double dt = dtinitital;
    for (int step = 0; step < nSteps; step++) {
		dt = computeCFL_dt(CFL);
        apply_boundary_conditions();
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    copyInitialState(globalGrid[i][j][k]);
                }
            }
        }
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    compute_time_derivatives(i, j, k);
                    double d_alpha_dt, d_beta_dt[3];
                    compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
                    storeStage(globalGrid[i][j][k], 0, d_alpha_dt, d_beta_dt);
                }
            }
        }
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 0);
                }
            }
        }
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    compute_time_derivatives(i, j, k);
                    double d_alpha_dt, d_beta_dt[3];
                    compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
                    storeStage(globalGrid[i][j][k], 1, d_alpha_dt, d_beta_dt);
                }
            }
        }
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    updateIntermediateState(globalGrid[i][j][k], 0.5 * dt, 1);
                }
            }
        }
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    compute_time_derivatives(i, j, k);
                    double d_alpha_dt, d_beta_dt[3];
                    compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
                    storeStage(globalGrid[i][j][k], 2, d_alpha_dt, d_beta_dt);
                }
            }
        }
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    updateIntermediateState(globalGrid[i][j][k], dt, 2);
                }
            }
        }
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    compute_time_derivatives(i, j, k);
                    double d_alpha_dt, d_beta_dt[3];
                    compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);
                    storeStage(globalGrid[i][j][k], 3, d_alpha_dt, d_beta_dt);
                }
            }
        }
        
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    combineStages(globalGrid[i][j][k], dt);
                }
            }
        }
        
        export_gamma_slice(NY / 2);
        export_gauge_slice(NY / 2);
		gridTensor.export_christoffel_slice(NY / 2);
    }
}
