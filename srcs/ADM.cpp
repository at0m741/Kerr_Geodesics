#include <Geodesics.h>


double second_derivative_alpha(int i, int j, double alpha, const Vector3& X) {
    double result = 0.0;
	if (i == j) {
		if (i == 0) {
			result = 0.0;
		} else if (i == 1) {
			result = 0.0;
		} else if (i == 2) {
			result = 0.0;
		}
	} else {
		result = 0.0;
	}
	return result;
}

Matrix3x3 RK4_update_gamma(const Matrix3x3& gamma, const Matrix3x3& dgamma_dt, double dt) {
    Matrix3x3 result;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result[i][j] = gamma[i][j] + dt * dgamma_dt[i][j];
    return result;
}

Matrix3x3 RK4_update_K(const Matrix3x3& K, const Matrix3x3& dK_dt, double dt) {
    Matrix3x3 result;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result[i][j] = K[i][j] + dt * dK_dt[i][j];
    return result;
}

double compute_rho(double r) {
    double R = 5.0;
    double rho0 = 1e-3;
    return (r < R) ? rho0 : 0.0;
}

void matter_sources(double r, double& rho, Vector3& j, Matrix3x3& S_ij, double& S) {
    rho = compute_rho(r);
    j = {0.0, 0.0, 0.0};
    for(int i=0; i<3; i++){
        for(int k=0; k<3; k++){
            S_ij[i][k] = 0.0;
        }
    }
    S = 0.0;
}

void evolveADM(const Matrix3x3& gamma, const Matrix3x3& K, double alpha, const Vector3& X, double dt, Matrix3x3& gamma_new, Matrix3x3& K_new) {
    Matrix3x3 R;
    Tensor3D Gamma3;
    Grid grid_obj;
    double rho, S;
	Vector3 j_i;
	Matrix3x3 S_ij;
	matter_sources(X[0], rho, j_i, S_ij, S);
    grid_obj.calculate_christoffel_3D(X, Gamma3);
    grid_obj.compute_ricci_3d(X, Gamma3, R);
    double K_trace = 0.0;
    for (int i = 0; i < 3; i++) {
        K_trace += K[i][i];
    }
    Matrix3x3 dgamma_dt;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            dgamma_dt[i][j] = -2.0 * alpha * K[i][j];
		}
	
    Matrix3x3 dK_dt;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double K2 = 0.0;
            for (int k = 0; k < 3; k++) {
                K2 += K[i][k] * K[k][j];
            }
            double D2alpha = second_derivative_alpha(i, j, alpha, X);
            dK_dt[i][j] = -D2alpha + alpha * (R[i][j] + K_trace * K[i][j] - 2.0 * K2) + 8.0 * M_PI * alpha * (S_ij[i][j] - 0.5 * gamma[i][j] * S);
		}
    }
    gamma_new = RK4_update_gamma(gamma, dgamma_dt, dt);
    K_new = RK4_update_K(K, dK_dt, dt);

}
