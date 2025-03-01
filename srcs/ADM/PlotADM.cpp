
#include <Geodesics.h>

void export_gamma_slice(int j) {
    std::ofstream file("gamma_slice_full.csv");

    file << "x,z,"
         << "gamma_00,gamma_01,gamma_02,"
         << "gamma_10,gamma_11,gamma_12,"
         << "gamma_20,gamma_21,gamma_22\n";
    
    for(int i = 0; i < NX; i++) {
        for(int k = 0; k < NZ; k++) {
            double x = -256.0 + i * (512.0 / (NX - 1));
            double z = -256.0 + k * (512.0 / (NZ - 1));

            double g00 = globalGrid[i][j][k].gamma[0][0];
            double g01 = globalGrid[i][j][k].gamma[0][1];
            double g02 = globalGrid[i][j][k].gamma[0][2];
            double g10 = globalGrid[i][j][k].gamma[1][0];
            double g11 = globalGrid[i][j][k].gamma[1][1];
            double g12 = globalGrid[i][j][k].gamma[1][2];
            double g20 = globalGrid[i][j][k].gamma[2][0];
            double g21 = globalGrid[i][j][k].gamma[2][1];
            double g22 = globalGrid[i][j][k].gamma[2][2];

            file << x << "," << z << ","
                 << g00 << "," << g01 << "," << g02 << ","
                 << g10 << "," << g11 << "," << g12 << ","
                 << g20 << "," << g21 << "," << g22 << "\n";
        }
    }
    file.close();
    std::cout << "Gamma slice (all components) saved to gamma_slice_full.csv\n";
}

void export_K_slice(int j) {
    std::ofstream file("K_slice.csv");
    file << "x,z,K00,K01,K02,K10,K11,K12,K20,K21,K22\n";
    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            double x = -256.0 + i * (512.0 / (NX - 1));
            double z = -256.0 + k * (512.0 / (NZ - 1));
            double K00 = globalGrid[i][j][k].K[0][0];
            double K01 = globalGrid[i][j][k].K[0][1];
            double K02 = globalGrid[i][j][k].K[0][2];
            double K10 = globalGrid[i][j][k].K[1][0];
            double K11 = globalGrid[i][j][k].K[1][1];
            double K12 = globalGrid[i][j][k].K[1][2];
            double K20 = globalGrid[i][j][k].K[2][0];
            double K21 = globalGrid[i][j][k].K[2][1];
            double K22 = globalGrid[i][j][k].K[2][2];
            file << x << "," << z << ","
                 << K00 << "," << K01 << "," << K02 << ","
                 << K10 << "," << K11 << "," << K12 << ","
                 << K20 << "," << K21 << "," << K22 << "\n";
        }
    }
    file.close();
    std::cout << "K slice saved to K_slice.csv\n";
}

void export_alpha_slice(int j) {
    std::ofstream file("alpha_slice.csv");

    file << "x,z,alpha\n";
    for(int i = 0; i < NX; i++) {
        for(int k = 0; k < NZ; k++) {
            double x = -256.0 + i * (512.0 / (NX - 1));
            double z = -256.0 + k * (512.0 / (NZ - 1));

            double alpha = globalGrid[i][j][k].alpha;
            file << x << "," << z << "," << alpha << "\n";
        }
    }
    file.close();
    std::cout << "Alpha slice saved to alpha_slice.csv\n";
}



void export_gauge_slice(int j) {
	Grid grid;
    std::ofstream file("gauge_slice.csv");
    file << "x,z,alpha,beta0,beta1,beta2,d_alpha_dt,d_beta0_dt,d_beta1_dt,d_beta2_dt\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            double x = -256.0 + i * (512.0 / (NX - 1));
            double z = -256.0 + k * (512.0 / (NZ - 1));

			Grid::Cell2D &cell = globalGrid[i][j][k];

            double d_alpha_dt, d_beta_dt[3];
            grid.compute_gauge_derivatives(i, j, k, d_alpha_dt, d_beta_dt);

            file << x << "," << z << "," 
                 << cell.alpha << "," << cell.beta[0] << "," << cell.beta[1] << "," << cell.beta[2] << ","
                 << d_alpha_dt << "," << d_beta_dt[0] << "," << d_beta_dt[1] << "," << d_beta_dt[2] << "\n";
        }
    }
    file.close();
    std::cout << "Gauge slice saved to gauge_slice.csv\n";
}


void GridTensor::export_christoffel_slice(int j) {
	GridTensor grid;
    std::ofstream file("christoffel_slice.csv");
    file << "x,z,Gamma000,Gamma001,Gamma002\n";
    
    for(int i = 1; i < NX-1; i++) {
        for(int k = 1; k < NZ-1; k++) {
            double x = -256.0 + i * (512.0 / (NX-1));
            double z = -256.0 + k * (512.0 / (NZ-1));

            double christof[3][3][3];
            grid.compute_christoffel_3D(i, j, k, christof);

            file << x << "," << z << ","
                 << christof[0][0][0] << ","
                 << christof[0][0][1] << ","
                 << christof[0][0][2] << "\n";
        }
    }
    
    file.close();
    std::cout << "Christoffel slice saved to christoffel_slice.csv\n";
}


void Grid::export_hamiltonian_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Erreur: impossible d'ouvrir " << filename << std::endl;
        return;
    }

    file << "x,y,z,Hamiltonian\n";
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                file << i * DX << "," << j * DY << "," << k * DZ << "," << hamiltonianGrid[i][j][k] << "\n";
            }
        }
    }
    file.close();
    std::cout << "Hamiltonian exporté dans " << filename << std::endl;
}

void Grid::export_fluid_slice(int j_slice) {
	double L = 2.0; 
    double x_min = -L, x_max = L;
    double y_min = -L, y_max = L;
    double z_min = -L, z_max = L;
    std::ofstream file("fluid_slice.csv");
    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier fluid_slice.csv" << std::endl;
        return;
    }

    file << "x,z,rho,p,vx,vy,vz\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            double x = x_min + i * DX;
            double z = z_min + k * DZ;

            Cell2D &cell = globalGrid[i][j_slice][k];
            file << x << "," << z << ","
                 << cell.rho << "," << cell.p << ","
                 << cell.vx << "," << cell.vy << "," << cell.vz
                 << "\n";
        }
    }

    file.close();
    std::cout << "✅ Export du fluide terminé dans fluid_slice.csv" << std::endl;
}


void Grid::export_energy_momentum_tensor_slice(int slice_y) {
    std::ofstream file("T_energy_momentum.csv");
    if (!file.is_open()) {
        std::cerr << "Erreur lors de l'ouverture du fichier !\n";
        return;
    }

    // En-têtes du CSV
    file << "x,z,T_00,T_01,T_02,T_10,T_11,T_12,T_20,T_21,T_22\n";

    for (int i = 0; i < NX; i++) {
        for (int k = 0; k < NZ; k++) {
            Cell2D &cell = globalGrid[i][slice_y][k];

            file << i << "," << k;
            for (int a = 0; a < 3; a++) {
                for (int b = 0; b < 3; b++) {
                    file << "," << cell.T[a][b];
                }
            }
            file << "\n";
        }
    }

    file.close();
    std::cout << "Exportation de T_{ab} terminée : T_energy_momentum.csv\n";
}
