#include <Geodesics.h>
extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

int grid_setup() {
    double r = 6.0;
    double theta = M_PI / 2.0;
    double phi = 0.0;

    Metric metric_obj;
    Grid grid_obj;

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

	Matrix3x3 Ricci2;
	grid_obj.compute_ricci_3d(X3D, Gamma3, Ricci2);
    
	printf("Hamiltonian constraint = %e\n", grid_obj.compute_hamiltonian_constraint(gamma3_inv, K, Ricci2));
    
	
	grid_obj.allocateGlobalGrid();
	grid_obj.initializeData();
	grid_obj.evolve(0.0000001, 1);
	printf("end of compute\n");
    return 0;
}
