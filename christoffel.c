#include "geodesics.h"

void christoffel(long double g[4][4], long double christoffel[4][4][4])
{
    for (int mu = 0; mu < 4; mu++) {
        for (int beta = 0; beta < 4; beta++) {
            for (int nu = 0; nu < 4; nu++) {
                long double sum = 0;
                for (int sigma = 0; sigma < 4; sigma++) {
                    sum += 0.5 * (g[mu][sigma] * (g[sigma][beta] + g[beta][sigma] - g[beta][nu]));
                }
                christoffel[mu][beta][nu] = sum;
                printf("Christoffel[%d][%d][%d] = %Lf\n", mu, beta, nu, christoffel[mu][beta][nu]);
            }
        }
    }
}

void riemann(long double g[4][4], long double christoffel[4][4][4], long double riemann[4][4][4][4])
{
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            for (int k = 1; k < 4; k++) {
                for (int l = 1; l < 4; l++) {
                    long double sum = 0;
                    sum += (1 / (2 * g[0][0])) * (christoffel[k][i][l] * christoffel[l][j][k] - christoffel[k][j][l] * christoffel[l][i][k]);
                    sum += (1 / (2 * g[0][0])) * (christoffel[k][i][l] * (g[l][j] * g[k][k] - g[l][k] * g[j][k]) - christoffel[k][j][l] * (g[l][i] * g[k][k] - g[l][k] * g[i][k]));
                    sum += (1 / (2 * g[0][0])) * (christoffel[l][i][k] * (g[l][j] * g[k][k] - g[l][k] * g[j][k]) - christoffel[l][j][k] * (g[l][i] * g[k][k] - g[l][k] * g[i][k]));
                    riemann[i][j][k][l] = sum;
					printf("Riemann[%d][%d][%d][%d] = %Lf\n", i, j, k, l, riemann[i][j][k][l]);
                }
            }
        }
    }
}

