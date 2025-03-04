#include <Geodesics.h>

static double dgamma[NDIM3][NDIM3][NDIM3];

void calc_gamma_ij(const Vector3& X3D, Matrix3x3& gamma3, Matrix3x3& gamma3_inv) {
    Vector4 X4D = { 0.0, X3D[0], X3D[1], X3D[2] };
    Matrix4x4 g{};    
    Matrix4x4 g_inv{}; 
    Metric metric;
    Matrix matrix_obj;
    metric.calculate_metric(X4D, g, g_inv);

    Matrix3x3 gamma3_arr{};
    for (int i = 0; i < NDIM3; i++) {
        for (int j = 0; j < NDIM3; j++) {
            gamma3_arr[i][j] = g[i+1][j+1];
        }
    }

    Matrix3x3 gamma3_inv_arr{};
    if (!matrix_obj.inverse_3x3(gamma3_arr, gamma3_inv_arr)) {
        printf("Erreur: gamma_{ij} est singulière ou mal définie\n");
    }

    gamma3 = gamma3_arr;
    gamma3_inv = gamma3_inv_arr;
}


void Grid::calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3, 
                                     const Matrix3x3& gamma, Matrix3x3 gamma_inv) {
    Matrix3x3 gamma_p{}, gamma_m{};
    Tensor3D dgamma{};
    Gamma3.fill({});  

    for (int m = 0; m < NDIM3; m++) {
        Vector3 Xp = X;
        Vector3 Xm = X;
        Xp[m] += DELTA3;
        Xm[m] -= DELTA3;

        calc_gamma_ij(Xp, gamma_p, gamma_inv);
        calc_gamma_ij(Xm, gamma_m, gamma_inv);
        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                dgamma[m][i][j] = (gamma_p[i][j] - gamma_m[i][j]) / (2.0 * DELTA3);
            }
        }
    }

    for (int k = 0; k < NDIM3; k++) {
        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                double sum = 0.0;
                for (int l = 0; l < NDIM3; l++) {
                    sum += gamma_inv[k][l] * ( dgamma[i][l][j] + dgamma[j][l][i] - dgamma[l][i][j] );
                }
                Gamma3[k][i][j] = 0.5 * sum;
            }
        }
    }

    const double tol = 1e-12;
    for (int k = 0; k < NDIM3; k++) {
        for (int i = 0; i < NDIM3; i++) {
            for (int j = 0; j < NDIM3; j++) {
                if (fabs(Gamma3[k][i][j] - Gamma3[k][j][i]) > tol) {
                    printf("Erreur: Gamma[%d][%d][%d] != Gamma[%d][%d][%d]\n", 
                           k, i, j, k, j, i);
                    printf("Gamma[%d][%d][%d] = %f, Gamma[%d][%d][%d] = %f\n",
                           k, i, j, Gamma3[k][i][j], k, j, i, Gamma3[k][j][i]);
                }
            }
        }
    }
}



