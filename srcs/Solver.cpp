#include <Geodesics.h>
#include <cassert>

 std::vector<std::vector<std::vector<Grid::Cell2D>>> globalGrid;

void Grid::allocateGlobalGrid(){
	printf("Allocating global grid\n");
    globalGrid.resize(NX);
    for(int i=0;i<NX;i++){
        globalGrid[i].resize(NY);
        for(int j=0;j<NY;j++){
            globalGrid[i][j].resize(NZ);
        }
    }
}

void Grid::initializeData() {
	double x_min = -128.0, x_max = 128.0;
	double y_min = -128.0, y_max = 128.0;
	double z_min = -128.0, z_max = 128.0;
	double dx = (x_max - x_min)/(NX-1);
	double dy = (y_max - y_min)/(NY-1);
	double dz = (z_max - z_min)/(NZ-1);

	for(int i=0; i<NX; i++){
		double x = x_min + i*dx;
		for(int j=0; j<NY; j++){
			double y = y_min + j*dy;
			for(int k=0; k<NZ; k++){
				double z = z_min + k*dz;


				double epsilon = 1e-7;
				double rho = sqrt(x*x + y*y + z*z);
				// Smoothing type :  rho_effectif = sqrt(rho^2 + epsilon^2)
				double rho_eff = sqrt(rho * rho + epsilon * epsilon);

				double Phi = 1.0 + 0.5 * M / rho_eff;
				Cell2D &cell = globalGrid[i][j][k];
				cell.alpha = (1.0 - M/(2*rho)) / (1.0 + M/(2*rho));
				cell.beta[0]=0; cell.beta[1]=0; cell.beta[2]=0;
				for(int a=0;a<3;a++){
					for(int b=0;b<3;b++){
						cell.gamma[a][b] = (a==b)? Phi*Phi*Phi*Phi : 0.0;
					}
				}
				for(int a=0;a<3;a++){
					for(int b=0;b<3;b++){
						cell.K[a][b] = 0.0;
					}
				}
			}
		}
	}


	double test_radii[] = {1.0, 2.0, 5.0, 10.0, 20.0, 128.0};
	printf("\nVérification des valeurs analytiques aux points clés :\n");
	for (double test_r : test_radii) {
		double Phi_test = 1.0 + 0.5 * M / test_r;
		double alpha_test = (1.0 - M / (2 * test_r)) / (1.0 + M / (2 * test_r));
		printf("r = %f : Phi^4 = %e, alpha = %e\n", test_r, pow(Phi_test, 4), alpha_test);
	}

	printf("\nVérification de la diagonale de γ_{ij} :\n");
	for (int iCell = 0; iCell < 3; iCell++) {
		for (int jCell = 0; jCell < 3; jCell++) {
			for (int kCell = 0; kCell < 3; kCell++) {
				Cell2D &cell = globalGrid[iCell][jCell][kCell];
				for (int a = 0; a < 3; a++) {
					for (int b = 0; b < 3; b++) {
						if (a != b && fabs(cell.gamma[a][b]) > 1e-12) {
							printf("⚠️ γ[%d][%d] non diagonale en (%d,%d,%d) : %e\n", 
									a, b, iCell, jCell, kCell, cell.gamma[a][b]);
						}
					}
				}
			}
		}
	}

	/* printf("\nTest sur les dérivées de γ_{ij} :\n"); */
	/* for (int iCell = 1; iCell < NX-1; iCell++) { */
	/* 	for (int jCell = 1; jCell < NY-1; jCell++) { */
	/* 		for (int kCell = 1; kCell < NZ-1; kCell++) { */
	/* 			double dgammaX = (globalGrid[iCell+1][jCell][kCell].gamma[0][0] -  */
	/* 					globalGrid[iCell-1][jCell][kCell].gamma[0][0]) / (2.0*dx); */
	/* 			double rho = sqrt((iCell*dx)*(iCell*dx) + (jCell*dy)*(jCell*dy) + (kCell*dz)*(kCell*dz)); */
	/* 			double Phi_test = 1.0 + 0.5 * M / rho; */
	/* 			double dPhiX = (-0.5 * M / (rho*rho)) * (iCell*dx/rho); */
	/* 			double expected_dgammaX = 4 * pow(Phi_test, 3) * dPhiX; */
	/* 			printf("Cell (%d,%d,%d): dγ/dx = %e, attendu = %e\n", iCell, jCell, kCell, dgammaX, expected_dgammaX); */
	/* 		} */
	/* 	} */
	/* } */
}


bool invert_3x3(const double m[3][3], double inv[3][3]) {
	double det =
		m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
		- m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
		+ m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);

	if (std::fabs(det) < 1e-14) return false;
	double idet = 1.0/det;

	inv[0][0] =  (m[1][1]*m[2][2]-m[2][1]*m[1][2]) * idet;
	inv[0][1] = -(m[0][1]*m[2][2]-m[2][1]*m[0][2]) * idet;
	inv[0][2] =  (m[0][1]*m[1][2]-m[1][1]*m[0][2]) * idet;

	inv[1][0] = -(m[1][0]*m[2][2]-m[2][0]*m[1][2]) * idet;
	inv[1][1] =  (m[0][0]*m[2][2]-m[2][0]*m[0][2]) * idet;
	inv[1][2] = -(m[0][0]*m[1][2]-m[1][0]*m[0][2]) * idet;

	inv[2][0] =  (m[1][0]*m[2][1]-m[2][0]*m[1][1]) * idet;
	inv[2][1] = -(m[0][0]*m[2][1]-m[2][0]*m[0][1]) * idet;
	inv[2][2] =  (m[0][0]*m[1][1]-m[1][0]*m[0][1]) * idet;

	return true;
}


void compute_christoffel_3D(int i, int j, int k, double christof[3][3][3]) {
	Matrix matrix_obj;
	double g[3][3];
	for(int a=0;a<3;a++){
		for(int b=0;b<3;b++){
			g[a][b] = globalGrid[i][j][k].gamma[a][b];
		}
	}
	double invg[3][3];
	bool ok = invert_3x3(g, invg);
	double dgamma[3][3][3];
	for(int a=0;a<3;a++){
		for(int b=0;b<3;b++){
			dgamma[0][a][b] = partialX_gamma(i,j,k,a,b);
			dgamma[1][a][b] = partialY_gamma(i,j,k,a,b);
			dgamma[2][a][b] = partialZ_gamma(i,j,k,a,b);
		}
	}
	for(int kk=0; kk<3; kk++){
		for(int aa=0; aa<3; aa++){
			for(int bb=0; bb<3; bb++){
				double sum=0.0;
				for(int ll=0; ll<3; ll++){
					double tmp = dgamma[aa][ll][bb] + dgamma[bb][ll][aa] - dgamma[ll][aa][bb];
					sum += invg[kk][ll]*tmp;
				}
				christof[kk][aa][bb] = 0.5 * sum;
			}
		}
	}
	printf("\nChristoffel Symbols:\n");
	for (int lambda = 0; lambda < NDIM; lambda++) {
		printf("\nGamma^%d:\n", lambda);
		for (int mu = 0; mu < NDIM; mu++) {
			for (int nu = 0; nu < NDIM; nu++) {
				printf("%12.6f\t", christof[lambda][mu][nu]);
			}
			printf("\n");
		}
	}
}


void compute_ricci_3D(int i, int j, int k, double Ricci[3][3])
{
    double Gamma[3][3][3];
    compute_christoffel_3D(i, j, k, Gamma);

    static double partialGamma[3][3][3][3];

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (i >= 2 && i <= NX - 3) {
            compute_christoffel_3D(i-2, j, k, Gmm);
            compute_christoffel_3D(i-1, j, k, Gm );
            compute_christoffel_3D(i+1, j, k, Gp );
            compute_christoffel_3D(i+2, j, k, Gpp);
        }
        else if (i >= 1 && i <= NX - 2) {
            compute_christoffel_3D(i-1, j, k, Gm );
            compute_christoffel_3D(i+1, j, k, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (i >= 2 && i <= NX - 3) {
                        partialGamma[0][kk][aa][bb] = (
                            -   Gpp[kk][aa][bb]
                            + 8.*Gp [kk][aa][bb]
                            - 8.*Gm [kk][aa][bb]
                            +    Gmm[kk][aa][bb]
                        ) / (12.0 * DX);

                    } else if (i >= 1 && i <= NX - 2) {
                        partialGamma[0][kk][aa][bb] = (
                            Gp[kk][aa][bb] - Gm[kk][aa][bb]
                        ) / (2.0 * DX);

                    } else {
                        partialGamma[0][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (j >= 2 && j <= NY - 3) {
            compute_christoffel_3D(i, j-2, k, Gmm);
            compute_christoffel_3D(i, j-1, k, Gm );
            compute_christoffel_3D(i, j+1, k, Gp );
            compute_christoffel_3D(i, j+2, k, Gpp);
        }
        else if (j >= 1 && j <= NY - 2) {
            compute_christoffel_3D(i, j-1, k, Gm );
            compute_christoffel_3D(i, j+1, k, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (j >= 2 && j <= NY - 3) {
                        partialGamma[1][kk][aa][bb] = (
                            -   Gpp[kk][aa][bb]
                            + 8.*Gp [kk][aa][bb]
                            - 8.*Gm [kk][aa][bb]
                            +    Gmm[kk][aa][bb]
                        ) / (12.0 * DY);

                    } else if (j >= 1 && j <= NY - 2) {
                        partialGamma[1][kk][aa][bb] = (
                            Gp[kk][aa][bb] - Gm[kk][aa][bb]
                        ) / (2.0 * DY);

                    } else {
                        partialGamma[1][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    {
        double Gmm[3][3][3], Gm[3][3][3], Gp[3][3][3], Gpp[3][3][3];

        if (k >= 2 && k <= NZ - 3) {
            compute_christoffel_3D(i, j, k-2, Gmm);
            compute_christoffel_3D(i, j, k-1, Gm );
            compute_christoffel_3D(i, j, k+1, Gp );
            compute_christoffel_3D(i, j, k+2, Gpp);
        }
        else if (k >= 1 && k <= NZ - 2) {
            compute_christoffel_3D(i, j, k-1, Gm );
            compute_christoffel_3D(i, j, k+1, Gp );
        }

        for (int kk = 0; kk < 3; kk++) {
            for (int aa = 0; aa < 3; aa++) {
                for (int bb = 0; bb < 3; bb++) {

                    if (k >= 2 && k <= NZ - 3) {
                        partialGamma[2][kk][aa][bb] = (
                            -   Gpp[kk][aa][bb]
                            + 8.*Gp [kk][aa][bb]
                            - 8.*Gm [kk][aa][bb]
                            +    Gmm[kk][aa][bb]
                        ) / (12.0 * DZ);

                    } else if (k >= 1 && k <= NZ - 2) {
                        partialGamma[2][kk][aa][bb] = (
                            Gp[kk][aa][bb] - Gm[kk][aa][bb]
                        ) / (2.0 * DZ);

                    } else {
                        partialGamma[2][kk][aa][bb] = 0.0;
                    }
                }
            }
        }
    }

    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double term1 = 0.0, term2 = 0.0, term3 = 0.0, term4 = 0.0;

            for (int kk = 0; kk < 3; kk++) {
                term1 += partialGamma[kk][kk][a][b];
                term2 += partialGamma[b][kk][a][kk];
            }

            for (int kk = 0; kk < 3; kk++) {
                for (int ll = 0; ll < 3; ll++) {
                    term3 += Gamma[kk][a][b] * Gamma[ll][kk][ll];
                    term4 += Gamma[ll][a][kk] * Gamma[kk][b][ll];
                }
            }
            Ricci[a][b] = term1 - term2 + term3 - term4;
        }
    }

    /* printf("Tenseur de Ricci au point (%d,%d,%d) :\n", i, j, k); */
    /* for (int a = 0; a < 3; a++) { */
    /*     for (int b = 0; b < 3; b++) { */
    /*         printf(" % .6e", Ricci[a][b]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
}

void Grid::compute_time_derivatives(int i, int j, int k)
{
    Cell2D &cell = globalGrid[i][j][k];

    double alpha = cell.alpha; 
	double beta[3] = { cell.beta[0], cell.beta[1], cell.beta[2] };
    double gammaLocal[3][3], KLocal[3][3];
    for(int a=0; a<3; a++){
        for(int b=0; b<3; b++){
            gammaLocal[a][b] = cell.gamma[a][b];
            KLocal[a][b]     = cell.K[a][b];
        }
    }
    double gammaInv[3][3];
    invert_3x3(gammaLocal, gammaInv);

    double Ktrace = 0.0;
    for(int a=0;a<3;a++){
        Ktrace += gammaInv[a][a]*KLocal[a][a];
    }

    double Gamma[3][3][3];
    compute_christoffel_3D(i,j,k,Gamma);
    double Ricci[3][3];
    compute_ricci_3D(i,j,k,Ricci);
    double dtGamma[3][3];
    auto partial_m_gamma = [&](int a, int b, int dim){ 
        if(dim==0) return partialX_gamma(i,j,k,a,b);
        if(dim==1) return partialY_gamma(i,j,k,a,b);
        if(dim==2) return partialZ_gamma(i,j,k,a,b);
        return 0.0;
    };
    auto partial_j_beta = [&](int jaxis, int compBeta){
        if(jaxis==0) return partialX_betacomp(i,j,k, compBeta);
        if(jaxis==1) return partialY_betacomp(i,j,k, compBeta);
        if(jaxis==2) return partialZ_betacomp(i,j,k, compBeta);
        return 0.0;
    };

    for(int a=0;a<3;a++){
        for(int b=0;b<3;b++){
            double conv = 0.0;
            for(int m=0;m<3;m++){
                conv += beta[m] * partial_m_gamma(a,b,m);
            }
            double shiftTerm = 0.0;
            {
                for(int m=0; m<3; m++){
                    shiftTerm += gammaLocal[a][m] * partial_j_beta(b,m);
                }
            }
            {
                for(int m=0; m<3; m++){
                    shiftTerm += gammaLocal[b][m] * partial_j_beta(a,m);
                }
            }

            dtGamma[a][b] = -2.0*alpha*KLocal[a][b] + conv + shiftTerm;
        }
    }

    for(int a=0;a<3;a++){
        for(int b=0;b<3;b++){
            cell.dgt[a][b] = dtGamma[a][b];
        }
    }

    double dtK[3][3];
    auto partial_m_alpha = [&](int dim){
        if(dim==0) return partialX_alpha(i,j,k);
        if(dim==1) return partialY_alpha(i,j,k);
        if(dim==2) return partialZ_alpha(i,j,k);
        return 0.0;
    };
       auto second_partial_alpha = [&](int iaxis, int jaxis){
        if(iaxis==0 && jaxis==0) return partialXX_alpha(i,j,k);
        if(iaxis==1 && jaxis==1) return partialYY_alpha(i,j,k);
        if(iaxis==2 && jaxis==2) return partialZZ_alpha(i,j,k);
        return 0.0;
    };

    auto D_iD_j_alpha = [&](int a, int b){
        double p2 = second_partial_alpha(a,b);
        double sumG=0.0;
        for(int m=0;m<3;m++){
            sumG += Gamma[m][a][b]* partial_m_alpha(m);
        }
        return p2 - sumG;
    };

    double Kcontra[3][3];
    for(int p=0;p<3;p++){
        for(int q=0;q<3;q++){
            double s=0.0;
            for(int r=0;r<3;r++){
                for(int s2=0;s2<3;s2++){
                    s += gammaInv[p][r]*gammaInv[q][s2]*KLocal[r][s2];
                }
            }
            Kcontra[p][q] = s;
        }
    }
    auto KiKkj = [&](int ii, int jj){
        double s=0.0;
        for(int rr=0; rr<3; rr++){
            s += KLocal[ii][rr]*Kcontra[rr][jj];
        }
        return s;
    };
    auto partial_m_K = [&](int a, int b, int dim){
        if(dim==0){
            return (globalGrid[i+1][j][k].K[a][b] - globalGrid[i-1][j][k].K[a][b])/(2.0*DX-1);
        } else if(dim==1){
            return (globalGrid[i][j+1][k].K[a][b] - globalGrid[i][j-1][k].K[a][b])/(2.0*DY-1);
        } else {
            return (globalGrid[i][j][k+1].K[a][b] - globalGrid[i][j][k-1].K[a][b])/(2.0*DZ-1);
        }
    };

    auto Lie_beta_K = [&](int a, int b){
        double conv=0.0;
        for(int m=0;m<3;m++){
            conv += beta[m]*partial_m_K(a,b,m);
        }
        double shiftPart=0.0;
        for(int m=0;m<3;m++){
            shiftPart += KLocal[m][b]* partial_j_beta(a,m);
            shiftPart += KLocal[a][m]* partial_j_beta(b,m);
        }
        return conv + shiftPart;
    };

    for(int a=0;a<3;a++){
        for(int b=0;b<3;b++){
            double termA = - D_iD_j_alpha(a,b); 
            double termB = alpha*( Ricci[a][b] 
                                   + Ktrace*KLocal[a][b]
                                   - 2.0*KiKkj(a,b) );
            double termC = Lie_beta_K(a,b);
            dtK[a][b] = termA + termB + termC;
        }
    }

    for(int a=0;a<3;a++){
        for(int b=0;b<3;b++){
            cell.dKt[a][b] = dtK[a][b];
        }
    }
}


void Grid::evolve(double dt, int nSteps) {
	printf("Evolve\n");
    for(int step=0; step<nSteps; step++){
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
					compute_time_derivatives(i,j,k); 
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gammaStage[0][a][b] = cell.dgt[a][b];
                            cell.KStage[0][a][b]     = cell.dKt[a][b];
                        }
                    }
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0; a<3; a++){
                        for(int b=0; b<3; b++){
                            cell.gamma[a][b] += 0.5 * dt * cell.gammaStage[0][a][b];
                            cell.K[a][b]     += 0.5 * dt * cell.KStage[0][a][b];
                        }
                    }
                }
            }
        }

        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    compute_time_derivatives(i,j,k);
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gammaStage[1][a][b] = cell.dgt[a][b];
                            cell.KStage[1][a][b]     = cell.dKt[a][b];
                        }
                    }
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gamma[a][b] -= 0.5 * dt * cell.gammaStage[0][a][b];
                            cell.K[a][b]     -= 0.5 * dt * cell.KStage[0][a][b];

                            cell.gamma[a][b] += 0.5 * dt * cell.gammaStage[1][a][b];
                            cell.K[a][b]     += 0.5 * dt * cell.KStage[1][a][b];
                        }
                    }
                }
            }
        }

        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    compute_time_derivatives(i,j,k);
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gammaStage[2][a][b] = cell.dgt[a][b];
                            cell.KStage[2][a][b]     = cell.dKt[a][b];
                        }
                    }
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
							cell.gamma[a][b] -= 0.5 * dt * cell.gammaStage[1][a][b];
							cell.K[a][b]     -= 0.5 * dt * cell.KStage[1][a][b];
						}
					}
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gamma[a][b] += dt * cell.gammaStage[2][a][b];
                            cell.K[a][b]     += dt * cell.KStage[2][a][b];
                        }
                    }
                }
            }
        }

        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    compute_time_derivatives(i,j,k);
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];
                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gammaStage[3][a][b] = cell.dgt[a][b];
                            cell.KStage[3][a][b]     = cell.dKt[a][b];
                        }
                    }
                }
            }
        }
        for(int i=1; i<NX-1; i++){
            for(int j=1; j<NY-1; j++){
                for(int k=1; k<NZ-1; k++){
                    Cell2D &cell = globalGrid[i][j][k];

                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            cell.gamma[a][b] -= dt*cell.gammaStage[2][a][b];
                            cell.K[a][b]     -= dt*cell.KStage[2][a][b];
                        }
                    }

                    for(int a=0;a<3;a++){
                        for(int b=0;b<3;b++){
                            double K1 = cell.gammaStage[0][a][b];
                            double K2 = cell.gammaStage[1][a][b];
                            double K3 = cell.gammaStage[2][a][b];
                            double K4 = cell.gammaStage[3][a][b];

                            cell.gamma[a][b] += (dt/6.0)*(K1 + 2.0*K2 + 2.0*K3 + K4);

                            double K1K = cell.KStage[0][a][b];
                            double K2K = cell.KStage[1][a][b];
                            double K3K = cell.KStage[2][a][b];
                            double K4K = cell.KStage[3][a][b];
                            cell.K[a][b] += (dt/6.0)*(K1K + 2.0*K2K + 2.0*K3K + K4K);
                        }
                    }
                }
            }
        }
        // end step

    } // fin de la boucle sur nSteps
}
