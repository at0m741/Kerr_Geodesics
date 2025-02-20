#include <Geodesics.h>

/* 
 * Richardson Extrapolation for the derivative of christoffel symbols
 * It is used to calculate the gamma derivative in the Riemann tensor
 */

double Tensor::richardson_derivative(
    const Tensor3D& Gamma_plus_h, 
    const Tensor3D& Gamma_minus_h,
    const Tensor3D& Gamma_plus_half_h,
    const Tensor3D& Gamma_minus_half_h,
    int rho, int mu, int nu, double h) 
{
    double diff_h = (Gamma_plus_h[rho][mu][nu] - Gamma_minus_h[rho][mu][nu]) / (2 * h);
    double diff_half_h = (Gamma_plus_half_h[rho][mu][nu] - Gamma_minus_half_h[rho][mu][nu]) / h;
    return (4 * diff_half_h - diff_h) / 3;
}

/* 
 * Calculate the Riemann tensor using the Christoffel symbols
 * The Riemann tensor is calculated using the formula:
 * R^rho_sigma_mu_nu = dGamma^rho_mu_nu/dx^sigma - dGamma^rho_nu/sigma
 * + Gamma^rho_mu_lambda * Gamma^lambda_nu_sigma - Gamma^rho_nu_lambda * Gamma^lambda_mu_sigma
 */

void Tensor::calculate_riemann(const Christoffel3D& Gamma, 
				const Christoffel4D& Gamma_plus_h, 
				const Christoffel4D& Gamma_minus_h, 
				const Christoffel4D& Gamma_plus_half_h, 
				const Christoffel4D& Gamma_minus_half_h,
				Riemann4D& Riemann, 
				double h) {
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

/* 
 * Contract the Riemann tensor to calculate the Ricci tensor
 * The Ricci tensor is calculated using the formula:
 * R_mu_nu = g^rho_sigma * R^sigma_rho_mu_nu
 */

void Tensor::contract_riemann(const Riemann4D& Riemann, MatrixNDIM& Ricci, const MatrixNDIM& g_inv) {
    for (auto &row : Ricci) {
        row.fill(0.0);
    }
    
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

void compute_partial_christoffel_3D(
    const double X[3],
    int m, 
    double dGamma[3][3][3], 
    double delta
)
{
    double Xp[3], Xm[3];
    memcpy(Xp, X, sizeof(Xp));
    memcpy(Xm, X, sizeof(Xm));

    Xp[m] += delta;
    Xm[m] -= delta;

    double Gamma_p[3][3][3], Gamma_m[3][3][3];
    memset(Gamma_p,0,sizeof(Gamma_p));
    memset(Gamma_m,0,sizeof(Gamma_m));

    calculate_christoffel_3D(Xp, Gamma_p);
    calculate_christoffel_3D(Xm, Gamma_m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                dGamma[i][j][k] = (Gamma_p[i][j][k] - Gamma_m[i][j][k]) / (2.0 * delta);
            }
        }
    }
}

void compute_ricci_3d(
    const double X[3],     
    double Gamma3[3][3][3], 
    double R3[3][3]         
)
{
    memset(R3,0,sizeof(double)*3*3);

    static double partialGamma[3][3][3][3];
    memset(partialGamma,0,sizeof(partialGamma));

    double delta = 1e-5;
    for(int m=0; m<3; m++){
        double dG[3][3][3];
        memset(dG,0,sizeof(dG));
        compute_partial_christoffel_3D(X,m,dG,delta);
        for(int k=0;k<3;k++)
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                    partialGamma[m][k][i][j] = dG[k][i][j];
    }

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            double term1=0, term2=0, term3=0, term4=0;
            for(int k=0;k<3;k++)
                term1 += partialGamma[k][k][i][j];
            for(int k=0;k<3;k++){
                term2 += partialGamma[j][k][i][k];
            for(int k=0;k<3;k++)
                for(int m=0;m<3;m++)
                    term3 += Gamma3[k][i][j]*Gamma3[m][k][m];
            for(int m=0;m<3;m++)
                for(int k=0;k<3;k++)
                    term4 += Gamma3[m][i][k]*Gamma3[k][j][m];
            }
            R3[i][j] = term1 - term2 + term3 - term4;
        }
    }

	print_ricci_tensor(R3);
}

void print_ricci_tensor(double R3[3][3]) {
	printf("\nRicci tensor:\n");
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%12.6f\t", R3[i][j]);
		}
		printf("\n");
	}
}

