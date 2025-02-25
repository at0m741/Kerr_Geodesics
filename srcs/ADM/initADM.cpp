#include <Geodesics.h>


void Grid::allocateGlobalGrid(){
	printf("Allocating global grid\n");
    globalGrid.resize(NX);
    for(int i=0;i<NX;i++){
        globalGrid[i].resize(NY);
        for(int j=0;j<NY;j++){
            globalGrid[i][j].resize(NZ);
        }
    }
}

void Grid::initializeData() {
	double x_min = -256.0, x_max = 256.0;
	double y_min = -256.0, y_max = 256.0;
	double z_min = -256.0, z_max = 256.0;
	double dx = (x_max - x_min)/(NX-1);
	double dy = (y_max - y_min)/(NY-1);
	double dz = (z_max - z_min)/(NZ-1);

	for(int i=0; i<NX; i++){
		double x = x_min + i*dx;
		for(int j=0; j<NY; j++){
			double y = y_min + j*dy;
			for(int k=0; k<NZ; k++){
				double z = z_min + k*dz;


				double epsilon = 1e-7;
				double rho = sqrt(x*x + y*y + z*z);
				rho = std::max(rho, epsilon); 
				double rho_eff = sqrt(rho * rho + epsilon * epsilon);

				double Phi = 1.0 + 0.5 * M / rho_eff;
				Cell2D &cell = globalGrid[i][j][k];
				cell.alpha = (1.0 - M/(2*rho_eff)) / (1.0 + M/(2*rho_eff));
				cell.beta[0]=0; cell.beta[1]=0; cell.beta[2]=0;
				for(int a=0;a<3;a++){
					for(int b=0;b<3;b++){
						cell.gamma[a][b] = (a==b)? Phi*Phi*Phi*Phi : 0.0;
					}
				}
				for(int a=0;a<3;a++){
					for(int b=0;b<3;b++){
						cell.K[a][b] = 0.0;
					}
				}
			}
		}
	}


	double test_radii[] = {1.0, 2.0, 5.0, 10.0, 20.0, 256.0};
	printf("\nVérification des valeurs analytiques aux points clés :\n");
	for (double test_r : test_radii) {
		double Phi_test = 1.0 + 0.5 * M / test_r;
		double alpha_test = (1.0 - M / (2 * test_r)) / (1.0 + M / (2 * test_r));
		printf("r = %f : Phi^4 = %e, alpha = %e\n", test_r, pow(Phi_test, 4), alpha_test);
	}

	printf("\nVérification de la diagonale de γ_{ij} :\n");
	for (int iCell = 0; iCell < 3; iCell++) {
		for (int jCell = 0; jCell < 3; jCell++) {
			for (int kCell = 0; kCell < 3; kCell++) {
				Cell2D &cell = globalGrid[iCell][jCell][kCell];
				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						if (a != b && fabs(cell.gamma[a][b]) > 1e-12) {
							printf("⚠️ γ[%d][%d] non diagonale en (%d,%d,%d) : %e\n", 
									a, b, iCell, jCell, kCell, cell.gamma[a][b]);
						}
					}
				}
			}
		}
	}
}




