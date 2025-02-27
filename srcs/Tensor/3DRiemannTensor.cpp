#include <Geodesics.h>

double Grid::richardson_derivative_ricci(
    const Tensor3D &Gamma_plus_h, 
    const Tensor3D &Gamma_minus_h,
    const Tensor3D &Gamma_plus_half_h,
    const Tensor3D &Gamma_minus_half_h,
    int mu, int nu, int sigma, 
    double h) 
{
    double diff_h = (Gamma_plus_h[mu][nu][sigma] - Gamma_minus_h[mu][nu][sigma]) / (2.0 * h);
    double diff_half_h = (Gamma_plus_half_h[mu][nu][sigma] - Gamma_minus_half_h[mu][nu][sigma]) / h;
    return (4.0 * diff_half_h - diff_h) / 3.0;
}

void Grid::calculate_riemann_3d(
    const Christoffel3D& Gamma, 
    const std::array<Christoffel3D, 3>& Gamma_plus_h,
    const std::array<Christoffel3D, 3>& Gamma_minus_h,
    const std::array<Christoffel3D, 3>& Gamma_plus_half_h,
    const std::array<Christoffel3D, 3>& Gamma_minus_half_h,
    Riemann3D& Riemann,
    double h,
    double scale  
) {
    double effective_h = h * scale;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double dGamma_k = richardson_derivative_ricci(
                        Gamma_plus_h[k], Gamma_minus_h[k],
                        Gamma_plus_half_h[k], Gamma_minus_half_h[k],
                        i, j, l, 
                        effective_h
                    );

                    double dGamma_l = richardson_derivative_ricci(
                        Gamma_plus_h[l], Gamma_minus_h[l],
                        Gamma_plus_half_h[l], Gamma_minus_half_h[l],
                        i, j, k, 
                        effective_h
                    );
                    double Gamma_terms = 0.0;
                    for (int m = 0; m < 3; m++) {
                        Gamma_terms += Gamma[i][k][m] * Gamma[m][l][j]
                                      - Gamma[i][l][m] * Gamma[m][k][j];
                    }
                    
                    Riemann[i][j][k][l] = dGamma_k - dGamma_l + Gamma_terms;
                }
            }
        }
    }
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				printf("Riemann[%d][%d][%d][%d] = %e\n", i, j, k, k, Riemann[i][j][k][k]);
			}
		}
	}
}

void Grid::calculate_riemann_4d_from_3d(
    const Riemann3D &Riemann3,
    const Matrix3x3 &K,     
    Riemann3D &Riemann4    
) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    Riemann4[i][j][k][l] = Riemann3[i][j][k][l]
                                          + K[i][k] * K[j][l]
                                          - K[i][l] * K[j][k];

                    if (i == j || k == l) {
                        Riemann4[i][j][k][l] = 0.0;  
                    }
                }
            }
        }
    }
    verify_riemann_symmetries(Riemann4);  
}
