#include <Geodesics.h>



void Grid::allocateGlobalGrid(){
    printf("Allocating global grid\n");

    globalGrid.resize(NX);
    for(int i = 0; i < NX; i++){
        globalGrid[i].resize(NY);
        for(int j = 0; j < NY; j++){
            globalGrid[i][j].resize(NZ);
        }
    }
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < NX; i++){
        for(int j = 0; j < NY; j++){
            for(int k = 0; k < NZ; k++){
                globalGrid[i][j][k] = Cell2D(); 
            }
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


double effective_potential(double x, double y, double a) {
    double r = sqrt(x*x + y*y);
    double cos_theta = y / (r + 1e-10);
    double sin2_theta = 1.0 - cos_theta * cos_theta;
    
    double delta = r*r - 2*r + a*a;
    double Sigma = r*r + a*a * cos_theta * cos_theta;
    double omega = 2.0 * r * a / (Sigma * delta);
    
    return 1.0 - sqrt(1.0 - 2.0 / r + a*a / (r*r));
}

void Grid::initializeData() {
    double L = 1.0; 
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    double dx = (x_max - x_min) / (NX - 1);
    double dy = (y_max - y_min) / (NY - 1);
    double dz = (z_max - z_min) / (NZ - 1);
    Matrix matrix;

    globalGrid.resize(NX);
    for (int i = 0; i < NX; i++) {
        globalGrid[i].resize(NY);
        for (int j = 0; j < NY; j++) {
            globalGrid[i][j].resize(NZ);
        }
    }

    #pragma omp parallel for collapse(3)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double z = z_min + k * dz;
                double r = sqrt(x * x + y * y + z * z);

                Cell2D cell;
                double Phi = 1.0 + 0.5 * M / r;
                cell.alpha = (1.0 - M / (2 * r)) / (1.0 + M / (2 * r));
                
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        cell.gamma[a][b] = (a == b) ? pow(Phi, 4) : 0.0;
                    }
                }

                matrix.inverse_3x3(cell.gamma, cell.gamma_inv);
                
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        cell.K[a][b] = 0.0;
                    }
                }

                cell.rho = exp(-r * r / 2.0);
                cell.p = 0.3 * cell.rho + 0.5 * cell.rho * cell.rho;

                double vr = 0.4;
                cell.vx = vr * x / r;
                cell.vy = vr * y / r;
                cell.vz = 0.0;

                if (r < 0.1) {
                    cell.vx = 0.0;
                    cell.vy = 0.0;
                    cell.vz = 0.0;
                    cell.rho = 0.0;
                    cell.p = 0.0;
                }

                globalGrid[i][j][k] = cell;

                if (j == NY / 2 && k == NZ / 2) { 
                    printf("gamma[0][0] at (i=%d, j=%d, k=%d) = %e\n", i, j, k, globalGrid[i][j][k].gamma[0][0]);
                }
            }
        }
    }
	double test_radii[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 128.0};
	for (double test_r : test_radii) {
		double Phi_test = 1.0 + 0.5 * M / test_r;
		double alpha_test = (1.0 - M / (2 * test_r)) / (1.0 + M / (2 * test_r));
		printf("r = %f : Phi^4 = %e, alpha = %e\n", test_r, pow(Phi_test, 4), alpha_test);
	}

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

	std::cout << "Test 1 (r -> ∞):\n";
	std::cout << "Alpha = " << cell_far.alpha << "\n";
	std::cout << "Gamma_xx = " << cell_far.gamma[0][0] << "\n";

	double rho_horizon = M / 2.0;  
	int i_horizon = static_cast<int>((rho_horizon - x_min) / dx);
	Cell2D &cell_horizon = globalGrid[i_horizon][j_center][k_center];

	std::cout << "\nTest 2 (r = M/2):\n";
	std::cout << "Alpha = " << cell_horizon.alpha << "\n";
	std::cout << "Position x = " << (x_min + i_horizon*dx) << "\n";
}




