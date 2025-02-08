#include "Geodesics.h"


double richardson_derivative(double (*Gamma_plus_h)[NDIM][NDIM], 
                            double (*Gamma_minus_h)[NDIM][NDIM],
                            double (*Gamma_plus_half_h)[NDIM][NDIM],
                            double (*Gamma_minus_half_h)[NDIM][NDIM],
                            int rho, int mu, int nu, double h) {
    double diff_h = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / (2 * h);
    double diff_half_h = (Gamma_plus_half_h[rho][mu][nu] - Gamma_minus_half_h[rho][mu][nu]) / h;
    return (4 * diff_half_h - diff_h) / 3;
}

void calculate_riemann(double Gamma[NDIM][NDIM][NDIM], 
                       double Gamma_plus_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_minus_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM], 
                       double Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM],
                       double Riemann[NDIM][NDIM][NDIM][NDIM], 
                       double h) {
    memset(Riemann, 0, sizeof(double) * NDIM * NDIM * NDIM * NDIM);
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    double dGamma_mu = richardson_derivative(
                        Gamma_plus_h[mu], Gamma_minus_h[mu], 
                        Gamma_plus_half_h[mu], Gamma_minus_half_h[mu], 
                        rho, nu, sigma, h);
                    
                    double dGamma_nu = richardson_derivative(
                        Gamma_plus_h[nu], Gamma_minus_h[nu],
                        Gamma_plus_half_h[nu], Gamma_minus_half_h[nu],
                        rho, mu, sigma, h);

                    double Gamma_terms = 0.0;
                    for (int lambda = 0; lambda < NDIM; lambda++) {
                        Gamma_terms += Gamma[rho][mu][lambda] * Gamma[lambda][nu][sigma]
                                     - Gamma[rho][nu][lambda] * Gamma[lambda][mu][sigma];
                    }

                    Riemann[rho][sigma][mu][nu] = dGamma_mu - dGamma_nu + Gamma_terms;
                }
            }
        }
    }
	double Kretschmann_scalar = 0.0;
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    Kretschmann_scalar += Riemann[rho][sigma][mu][nu] * \
										  Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
    if (Kretschmann_scalar > 1e10) {
        printf("Kretschmann Scalar: INF (Singularity detected)\n");
    } else {
        printf("Kretschmann Scalar: %12.6f\n", Kretschmann_scalar);
    }
}

void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM],\
					  double Ricci[NDIM][NDIM], 
					  double g_inv[NDIM][NDIM]) {
    memset(Ricci, 0, sizeof(double) * NDIM * NDIM);
	 for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    Ricci[mu][nu] += g_inv[rho][sigma] * Riemann[rho][sigma][mu][nu];
                }
            }
        }
    }
    printf("\nRicci tensor:\n");
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            printf("%12.6f\t", Ricci[mu][nu]);
        }
        printf("\n");
    }
    
    double Ricci_scalar = 0.0;
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            Ricci_scalar += g_inv[mu][nu] * Ricci[mu][nu];
        }
    }
    printf("Ricci Scalar: %12.6f\n", Ricci_scalar);
}
