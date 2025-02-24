#include <Geodesics.h>
#include <chrono>
#include <fstream>
extern double (*geodesic_points)[5];
extern int num_points;
extern double a;



int Riemann_tensor(const char *metric) {
    Tensor tensor;
    Matrix matrix_obj;
    Connexion connexion;
    Metric metric_obj;
    
    double r0 = 20.0;
    std::array<double, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};

    if (strcmp(metric, "ds") == 0) {
        printf("KDS metric calculation\n");
        metric_obj.calculate_metric_kds(X, metric_obj.gcov, metric_obj.gcon);
    } else {
        metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
        matrix_obj.print_matrix("metric_obj.gcov", metric_obj.gcov);
        matrix_obj.print_matrix("gcon", metric_obj.gcon);
    }

    connexion.calculate_christoffel(X, DELTA, connexion.Gamma, metric_obj.gcov, metric_obj.gcon, "kerr");

    for (int d = 0; d < NDIM; d++) {
        tensor.calculate_Gamma_at_offset(X, d, DELTA, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_plus_h[d], "kerr");
        tensor.calculate_Gamma_at_offset(X, d, -DELTA, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_minus_h[d], "kerr");
    }

    for (int d = 0; d < NDIM; d++) {
        tensor.calculate_Gamma_at_offset(X, d, DELTA/2.0, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_plus_half_h[d], "kerr");
        tensor.calculate_Gamma_at_offset(X, d, -DELTA/2.0, DELTA, metric_obj.gcov, metric_obj.gcon, connexion.Gamma_minus_half_h[d], "kerr");
    }

    tensor.calculate_riemann(connexion.Gamma, connexion.Gamma_plus_h,
                             connexion.Gamma_minus_h, connexion.Gamma_plus_half_h,
                             connexion.Gamma_minus_half_h, tensor.Riemann, DELTA);

    std::cout << "Riemann tensor compute finished" << std::endl;
    tensor.print_riemann(tensor.Riemann);
    tensor.contract_riemann(tensor.Riemann, tensor.Ricci, metric_obj.gcon);
    
    return 0;
}



int Geodesics_prob() {
	Connexion connexion;
	Metric metric_obj;
	double r0 = 20.0;
    std::array<double, NDIM> X = {0.0, r0, M_PI/3.0, 0.0};;
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
	double g_tt = metric_obj.gcov[0][0];
	double g_tphi = metric_obj.gcov[0][3];
	double g_phiphi = metric_obj.gcov[3][3];
	double Omega = 1.0 / (pow(r0, 1.5) + a);
	double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
	double vt = 1.0 / sqrt(denom);
	double v[NDIM] = {vt, 0.0, 0.0, Omega * vt};
	double norm = g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
	printf("Norme quadrivecteur = %e (doit être proche de 0)\n", norm);
	double dt = 0.0910;
	double christoffel[NDIM][NDIM][NDIM];
	connexion.calculate_christoffel(X, DELTA, connexion.Gamma, metric_obj.gcov, metric_obj.gcon, "kerr");

	__m256d X_avx[NDIM], v_avx[NDIM];
	for (int i = 0; i < NDIM; i++) {
		X_avx[i] = _mm256_set1_pd(X[i]);
		v_avx[i] = _mm256_set1_pd(v[i]);
	}

	__m256d christoffel_avx[NDIM][NDIM][NDIM];
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				christoffel_avx[i][j][k] = _mm256_set1_pd(connexion.Gamma[i][j][k]);
				printf("christoffel_avx[%d][%d][%d] = %e\n", i, j, k, connexion.Gamma[i][j][k]);
			}
		}
	}

	printf("begin geodesic\n");
	auto start = std::chrono::high_resolution_clock::now();
	geodesic_AVX(X_avx, v_avx, max_dt + 4, ( __m256d (*)[NDIM][NDIM] )christoffel_avx, _mm256_set1_pd(dt));
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	printf("Elapsed time: %f\n", elapsed_seconds.count());
	write_vtk_file("output/geodesic.vtk");
	if (geodesic_points != NULL) {
		free(geodesic_points);
	}
	return 0;
}

int light_geodesics_prob() {
	Connexion connexion;
	Metric metric_obj;
	double r0 = 20.0;
    std::array<double, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};;
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);
	double g_tt = metric_obj.gcov[0][0];
	double g_tphi = metric_obj.gcov[0][3];
	double g_phiphi = metric_obj.gcov[3][3];
	double Omega = 1.0 / (pow(r0, 1.5) + a);
	double denom = -(g_tt + 2.0 * g_tphi * Omega + g_phiphi * Omega * Omega);
	double vt = 1.0 / sqrt(fabs(denom));
	double v[NDIM] = {vt, 0.0, 0.0, 3.5 * Omega * vt};
	double norm = g_tt * v[0] * v[0] + 2.0 * g_tphi * v[0] * v[3] + g_phiphi * v[3] * v[3];
	double dt = 0.0910;
	double christoffel[NDIM][NDIM][NDIM];
	connexion.calculate_christoffel(X, DELTA, connexion.Gamma, metric_obj.gcov, metric_obj.gcon, "kerr");
	__m256d X_avx[NDIM], v_avx[NDIM];
	for (int i = 0; i < NDIM; i++) {
		X_avx[i] = _mm256_set1_pd(X[i]);
		v_avx[i] = _mm256_set1_pd(v[i]);
	}
	__m256d christoffel_avx[NDIM][NDIM][NDIM];
	for (int i = 0; i < NDIM; i++) {
		for (int j = 0; j < NDIM; j++) {
			for (int k = 0; k < NDIM; k++) {
				christoffel_avx[i][j][k] = _mm256_set1_pd(connexion.Gamma[i][j][k]);
			}
		}
	}

	auto start = std::chrono::high_resolution_clock::now();
	geodesic_AVX(X_avx, v_avx, max_dt + 4, ( __m256d (*)[NDIM][NDIM] )christoffel_avx, _mm256_set1_pd(dt));
	auto end = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	printf("Elapsed time: %f\n", elapsed_seconds.count());
	write_vtk_file("output/light_geodesic.vtk");
	if (geodesic_points != NULL) {
		free(geodesic_points);
	}
	return 0;
}

int Metric_prob() {
	Metric metric_obj;
	double r0 = 20.0;
    std::array<double, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};;
	printf("=====================================================\n");
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);

	printf("=====================================================\n");
	std::array<std::array<double, NDIM>, NDIM> gcov_KN, gcon_KN;	
	metric_obj.calculate_metric_kerr_newman(X, gcov_KN, gcon_KN);

	printf("=====================================================\n");
	std::array<std::array<double, NDIM>, NDIM> gcov_kds, gcon_kds;
	metric_obj.calculate_metric_kds(X, gcov_kds, gcon_kds);
	return 0;
}



void Grid::save_christoffel_symbols(const std::vector<std::vector<Cell2D>>& grid,
                                      int Nx, int Ny,
                                      const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier " << filename << "\n";
        return;
    }
    
    file << "Cellule,i,j,k,l,m,Gamma3,Ricci,Kij\n";
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        file << i << "," << j << "," << k << "," << l << "," << m << ","
                             << std::setprecision(12) << grid[i][j].Gamma3[k][l][m] << ","
							 << std::setprecision(12) << grid[i][j].Ricci[k][l]  << ","
							 << std::setprecision(12) << grid[i][j].K[k][l] << "\n";
                    }
                }
            }
        }
    }
    
    file.close();
    std::cout << "Les symboles de Christoffel ont été enregistrés dans " << filename << "\n";
}

void Grid::initialize_grid(int Nr, int Ntheta, double r_min, double r_max, double theta_min, double theta_max) {
    this->Nr = Nr;
    this->Ntheta = Ntheta;
    this->dr = (r_max - r_min) / (Nr - 1);
    this->dtheta = (theta_max - theta_min) / (Ntheta - 1);
    grid.resize(Nr, std::vector<Cell2D>(Ntheta));
    Metric metric_obj;
    double h_riemann = 1e-5;

    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Ntheta; j++) {
            double r_i = r_min + i * dr;
            double th_j = theta_min + j * dtheta;
            std::array<double, 4> X4D_local = {0.0, r_i, th_j, 0.0};
            Vector3 X3D_local = {r_i, th_j, 0.0};

            Matrix4x4 gcov, gcon;
            metric_obj.calculate_metric(X4D_local, gcov, gcon);

            extract_3p1(gcov, gcon, &grid[i][j].alpha,
                        grid[i][j].beta_cov, grid[i][j].beta_con,
                        grid[i][j].gamma, grid[i][j].gamma_inv);

            calculate_christoffel_3D(X3D_local, grid[i][j].Gamma3, grid[i][j].gamma, grid[i][j].gamma_inv);

            Christoffel3D Gamma = grid[i][j].Gamma3;
            Christoffel3D Gamma_plus_h, Gamma_minus_h, Gamma_plus_half_h, Gamma_minus_half_h;
            {
                Vector3 Xp = X3D_local, Xm = X3D_local;
                Xp[0] += h_riemann;
                Xm[0] -= h_riemann;
                calculate_christoffel_3D(Xp, Gamma_plus_h, grid[i][j].gamma, grid[i][j].gamma_inv);
                calculate_christoffel_3D(Xm, Gamma_minus_h, grid[i][j].gamma, grid[i][j].gamma_inv);

                Vector3 Xp_half = X3D_local, Xm_half = X3D_local;
                Xp_half[0] += h_riemann / 2.0;
                Xm_half[0] -= h_riemann / 2.0;
                calculate_christoffel_3D(Xp_half, Gamma_plus_half_h, grid[i][j].gamma, grid[i][j].gamma_inv);
                calculate_christoffel_3D(Xm_half, Gamma_minus_half_h, grid[i][j].gamma, grid[i][j].gamma_inv);
            }

            Riemann3D Riemann;
            calculate_riemann_3d(Gamma, Gamma_plus_h, Gamma_minus_h,
                                 Gamma_plus_half_h, Gamma_minus_half_h,
                                 Riemann, h_riemann, 1.0);
			for (int k = 0; k < DIM3; k++) {
				for (int l = 0; l < DIM3; l++) {
					for (int m = 0; m < DIM3; m++) {
						printf("Riemann[%d][%d][%d][%d] = %e\n", k, l, m, m, Riemann[k][l][m][m]);
					}
				}
			}
            Matrix3x3 Ricci;
            calculate_ricci_3d_from_riemann(Riemann, Ricci);
            grid[i][j].Ricci = Ricci;
			print_ricci_tensor(Ricci);
            double H_constraint = compute_hamiltonian_constraint(grid[i][j].gamma_inv, grid[i][j].K, grid[i][j].Ricci);
            printf("r = %e, theta = %e, Hamiltonian constraint = %e\n", r_i, th_j, H_constraint);
        }
    }
    evolve_Kij(1e-8);
}

int grid_setup() {
    double r = 6.0;
    double theta = M_PI / 2.0;
    double phi = 0.0;

    Metric metric_obj;
    Grid grid_obj;
    Matrix matrix_obj; 

    Vector3 X3D = { r, theta, phi };  
    std::array<double, NDIM> X4D = { 0.0, r, theta, phi };

    metric_obj.calculate_metric(X4D, metric_obj.gcov, metric_obj.gcon);
    double alpha;
    Vector3 beta_cov, beta_con;
    Matrix3x3 gamma3, gamma3_inv;
    grid_obj.extract_3p1(metric_obj.gcov, metric_obj.gcon, &alpha, beta_cov, beta_con, gamma3, gamma3_inv);
    
    Tensor3D Gamma3;
    grid_obj.calculate_christoffel_3D(X3D, Gamma3, gamma3, gamma3_inv);
    
    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            for (int k = 0; k < DIM3; k++) {
                printf("Gamma3[%d][%d][%d] = %e\n", i, j, k, Gamma3[i][j][k]);
            }
        }
    }
    Matrix3x3 dbeta;
    grid_obj.calculate_dbeta(X3D, dbeta);

    Matrix3x3 K;
    grid_obj.compute_extrinsic_curvature_stationary_3D(X3D, alpha, beta_cov, Gamma3, dbeta, K);

    double K_trace = grid_obj.compute_K(gamma3_inv, K);
    double KijKij = grid_obj.compute_Kij_Kij(gamma3_inv, K);

    for (int i = 0; i < DIM3; i++) {
        for (int j = 0; j < DIM3; j++) {
            printf("K[%d][%d] = %e\n", i, j, K[i][j]);
        }
    }
    printf("Trace K = %e\n", K_trace);
    printf("Kij K^ij = %e\n", KijKij);
    

    double h_riemann = 1e-5;
    Riemann3D Riemann;
    Christoffel3D Gamma_plus_h, Gamma_minus_h, Gamma_plus_half_h, Gamma_minus_half_h;
    {
        Vector3 Xp = X3D, Xm = X3D, Xp_half = X3D, Xm_half = X3D;
        Xp[0] += h_riemann;
        Xm[0] -= h_riemann;
        Xp_half[0] += h_riemann / 2.0;
        Xm_half[0] -= h_riemann / 2.0;
        grid_obj.calculate_christoffel_3D(Xp, Gamma_plus_h, gamma3, gamma3_inv);
        grid_obj.calculate_christoffel_3D(Xm, Gamma_minus_h, gamma3, gamma3_inv);
        grid_obj.calculate_christoffel_3D(Xp_half, Gamma_plus_half_h, gamma3, gamma3_inv);
        grid_obj.calculate_christoffel_3D(Xm_half, Gamma_minus_half_h, gamma3, gamma3_inv);
    }
    grid_obj.calculate_riemann_3d(Gamma3, Gamma_plus_h, Gamma_minus_h,
                         Gamma_plus_half_h, Gamma_minus_half_h,
                         Riemann, h_riemann, DELTA);
	for (int k = 0; k < DIM3; k++) {
		for (int l = 0; l < DIM3; l++) {
			for (int m = 0; m < DIM3; m++) {
				printf("Riemann[%d][%d][%d][%d] = %e\n", k, l, m, m, Riemann[k][l][m][m]);
			}
		}
	}

    Matrix3x3 Ricci;
    grid_obj.calculate_ricci_3d_from_riemann(Riemann, Ricci);
    /*  */
    /* grid_obj.initialize_grid(20, 20, 6.0, 12.0, 0.0, M_PI); */
    /* 	 */
    /* grid_obj.calculate_christoffel_3D_grid(grid_obj.grid, grid_obj.Nr, grid_obj.Ntheta, grid_obj.dr, grid_obj.dtheta, 2.0, 0.0); */
    /* grid_obj.compute_ricci_3d_grid(grid_obj.grid, grid_obj.Nr, grid_obj.Ntheta, 2.0, 0.0, grid_obj.dr, grid_obj.dtheta, 1e-12); */
    /* grid_obj.save_christoffel_symbols(grid_obj.grid, grid_obj.Nr, grid_obj.Ntheta, "output/christoffel.csv"); */
    /*  */
    grid_obj.evolve_Kij(1e-8);
    
    return 0;
}
