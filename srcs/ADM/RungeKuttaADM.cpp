#include <Geodesics.h>



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
					update_fluid_velocity(i, j, k, dt);
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
					update_fluid_velocity(i, j, k, dt);
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
					update_fluid_velocity(i, j, k, dt);
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
					update_fluid_velocity(i, j, k, dt);
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
		export_fluid_slice(NY / 2);
    }
}
