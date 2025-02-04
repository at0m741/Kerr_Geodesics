#include "Geodesics.h"

void initialize_riemann_tensor(double R[NDIM][NDIM][NDIM][NDIM]) {
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            for (int rho = 0; rho < NDIM; rho++) {
                for (int sigma = 0; sigma < NDIM; sigma++) {
                    R[mu][nu][rho][sigma] = 0.0; 
                }
            }
        }
    }
}


void print_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM]) {
    const double threshold = 1e-10; 
    printf("\nRiemann Tensor (Non-zero components):\n"); 
    for (int rho = 0; rho < NDIM; rho++) {
        for (int sigma = 0; sigma < NDIM; sigma++) {
            for (int mu = 0; mu < NDIM; mu++) {
                for (int nu = 0; nu < NDIM; nu++) {
                    double value = Riemann[rho][sigma][mu][nu];
                    if (fabs(value) > threshold) {
                        printf("Riemann[%d][%d][%d][%d] = %12.6f\n", rho, sigma, mu, nu, value);
                    }
                }
            }
        }
    }
}

void print_christoffel_matrix(double gamma[NDIM][NDIM][NDIM]) {
    printf("\nChristoffel Symbols:\n");
    for (int lambda = 0; lambda < NDIM; lambda++) {
        printf("\nGamma^%d:\n", lambda);
        for (int mu = 0; mu < NDIM; mu++) {
            for (int nu = 0; nu < NDIM; nu++) {
                printf("%12.6f\t", gamma[lambda][mu][nu]);
            }
            printf("\n");
        }
    }
}
