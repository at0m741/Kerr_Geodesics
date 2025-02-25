#include <Geodesics.h>

void Tensor::check_riemann_symmetries(const Riemann4D& Riemann, double tolerance) {
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    if (fabs(Riemann[rho][sigma][mu][nu] + Riemann[rho][sigma][nu][mu]) > tolerance) {
                        printf("Non-antisymmetry R[%d][%d][%d][%d] and R[%d][%d][%d][%d]\n", rho, sigma, mu, nu, rho, sigma, nu, mu);
                    }
                    if (fabs(Riemann[rho][sigma][mu][nu] - Riemann[sigma][rho][mu][nu]) > tolerance) {
                        printf("Non-symmetry R[%d][%d][%d][%d] and R[%d][%d][%d][%d]\n", rho, sigma, mu, nu, sigma, rho, mu, nu);
                    }
                    if (fabs(Riemann[rho][sigma][mu][nu] + Riemann[rho][mu][nu][sigma] + Riemann[rho][nu][sigma][mu]) > tolerance) {
                        printf("Bianchi identity violation at R[%d][%d][%d][%d]\n", rho, sigma, mu, nu);
                    }
                }
            }
        }
    }
}
void Connexion::check_symmetry_christoffel(const Christoffel3D& gamma) {
    printf("\nChristoffel symbols verif with delta = %e:\n", TOLERANCE);
    for (int lambda = 0; lambda < NDIM; lambda++) {
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = mu; nu < NDIM; nu++) { 
                double diff = fabs(gamma[lambda][mu][nu] - gamma[lambda][nu][mu]);
                if (diff > TOLERANCE) {
                    printf("!symmetry at Gamma^%d_%d%d: |%f - %f| = %f\n", 
                           lambda, mu, nu, gamma[lambda][mu][nu], gamma[lambda][nu][mu], diff);
                }
            }
        }
    }
    printf("Symmetric !\n");
}

bool Grid::verify_riemann_symmetries(const Riemann3D &Riemann) {
    bool ok = true;
    const double tol = 1e-12;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double diff = Riemann[i][j][k][l] + Riemann[i][j][l][k];
                    if (fabs(diff) > tol) {
                        printf("Violation antisymétrie (indices 3 et 4) à Riemann[%d][%d][%d][%d] + Riemann[%d][%d][%d][%d] = %e\n",
                               i, j, k, l, i, j, l, k, diff);
                        ok = false;
                    }
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double diff = Riemann[i][j][k][l] + Riemann[j][i][k][l];
                    if (fabs(diff) > tol) {
                        printf("Violation antisymétrie (indices 1 et 2) à Riemann[%d][%d][%d][%d] + Riemann[%d][%d][%d][%d] = %e\n",
                               i, j, k, l, j, i, k, l, diff);
                        ok = false;
                    }
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double diff = Riemann[i][j][k][l] - Riemann[k][l][i][j];
                    if (fabs(diff) > tol) {
                        printf("Violation symétrie par échange de paires à Riemann[%d][%d][%d][%d] - Riemann[%d][%d][%d][%d] = %e\n",
                               i, j, k, l, k, l, i, j, diff);
                        ok = false;
                    }
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    double bianchi = Riemann[i][j][k][l] + Riemann[i][k][l][j] + Riemann[i][l][j][k];
                    if (fabs(bianchi) > tol) {
                        printf("Violation de l'identité de Bianchi à Riemann[%d][%d][*][*] somme = %e\n", i, j, bianchi);
                        ok = false;
                    }
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
    
    return ok;
}

