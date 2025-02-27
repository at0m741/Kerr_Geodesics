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


void Grid::initializeData_Minkowski()
{
    printf("\n=== Initialisation Minkowski ===\n");

    double x_min = -128.0, x_max = 128.0;
    double y_min = -128.0, y_max = 128.0;
    double z_min = -128.0, z_max = 128.0;
    double dx = (x_max - x_min)/(NX-1);
    double dy = (y_max - y_min)/(NY-1);
    double dz = (z_max - z_min)/(NZ-1);

    for(int i=0; i<NX; i++)
    {
        double x = x_min + i*dx;
        for(int j=0; j<NY; j++)
        {
            double y = y_min + j*dy;
            for(int k=0; k<NZ; k++)
            {
                Cell2D &cell = globalGrid[i][j][k];

                cell.alpha = 1.0;
                cell.beta[0] = 0.0;
                cell.beta[1] = 0.0;
                cell.beta[2] = 0.0;

                for(int a=0; a<3; a++)
                {
                    for(int b=0; b<3; b++)
                    {
                        cell.gamma[a][b] = (a == b) ? 1.0 : 0.0;
                    }
                }

                for(int a=0; a<3; a++){
                    for(int b=0; b<3; b++){
                        cell.K[a][b] = 0.0;
                    }
                }
            }
        }
    }

    printf("Vérification Minkowski sur quelques points:\n");
    for(int test_i=0; test_i<3; test_i++)
    {
        for(int test_j=0; test_j<3; test_j++)
        {
            for(int test_k=0; test_k<3; test_k++)
            {
                Cell2D &cell = globalGrid[test_i][test_j][test_k];
                printf("Point (%d,%d,%d): alpha=%f, gamma=\n", test_i, test_j, test_k, cell.alpha);
                for(int a=0; a<3; a++)
                {
                    printf("  ");
                    for(int b=0; b<3; b++)
                    {
                        printf("%f ", cell.gamma[a][b]);
                    }
                    printf("\n");
                }
            }
        }
    }
    printf("=== Fin initialisation Minkowski ===\n");
}


void Grid::initializeData() {
	double L = 256.0; 
    double x_min = 2, x_max = L;
    double y_min = 2, y_max = L;
    double z_min = 2, z_max = L;
    double dx = (x_max - x_min)/(NX-1);
    double dy = (y_max - y_min)/(NY-1);
    double dz = (z_max - z_min)/(NZ-1);
	Matrix matrix;
    double epsilon = 1e-6;

    for(int i=0; i<NX; i++){
        double x = x_min + i*dx;
        for(int j=0; j<NY; j++){
            double y = y_min + j*dy;
            for(int k=0; k<NZ; k++){
                double z = z_min + k*dz;

                double rho = sqrt(x*x + y*y + z*z);
                Cell2D &cell = globalGrid[i][j][k];
                
				double Phi = 1.0 + 0.5 * M / rho;
				cell.alpha = (1.0 - M/(2*rho)) / (1.0 + M/(2*rho));
				for(int a=0; a<3; a++){
					for(int b=0; b<3; b++){
						cell.gamma[a][b] = (a == b) ? pow(Phi, 4) : 0.0;
					}
				}
				matrix.inverse_3x3(cell.gamma, cell.gamma_inv);
                for(int a=0; a<3; a++){
                    for(int b=0; b<3; b++){
                        cell.K[a][b] = 0.0;
                    }
                }
				
				if (j == NY / 2 && k == NZ / 2) { 
					printf("gamma[0][0] at (i=%d, j=%d, k=%d) = %e\n", i, j, k, globalGrid[i][j][k].gamma[0][0]);
				}
            }
        }
    }
	double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 128.0};
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
	int i_far = NX - 1;  
	int j_center = NY / 2;
	int k_center = NZ / 2;
	Cell2D &cell_far = globalGrid[i_far][j_center][k_center];

	std::cout << "Test 1 (rho -> ∞):\n";
	std::cout << "Alpha = " << cell_far.alpha << " (devrait être ~1)\n";
	std::cout << "Gamma_xx = " << cell_far.gamma[0][0] << " (devrait être ~1)\n";

	double rho_horizon = M / 2.0;  
	int i_horizon = static_cast<int>((rho_horizon - x_min) / dx);
	Cell2D &cell_horizon = globalGrid[i_horizon][j_center][k_center];

	std::cout << "\nTest 2 (rho = M/2):\n";
	std::cout << "Alpha = " << cell_horizon.alpha << " (devrait être ~0)\n";
	std::cout << "Position x = " << (x_min + i_horizon*dx) << "\n";
}




