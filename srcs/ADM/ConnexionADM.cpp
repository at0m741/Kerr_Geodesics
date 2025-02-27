#include <Geodesics.h>


void GridTensor::compute_christoffel_3D(int i, int j, int k, double christof[3][3][3]) {
	Matrix matrix_obj;
	double g[3][3];
	double invg[3][3];  
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			invg[a][b] = globalGrid[i][j][k].gamma_inv[a][b];
			g[a][b] = globalGrid[i][j][k].gamma[a][b];
		}
	}

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
					sum += invg[kk][ll] * tmp;
				}
				christof[kk][aa][bb] = 0.5 * sum;
			}
		}
	}
	/* printf("Christoffel symbols at (%d, %d, %d):\n", i, j, k); */
	/* for(int a=0; a<3; a++){ */
	/* 	for(int b=0; b<3; b++){ */
	/* 		for(int c=0; c<3; c++){ */
	/* 			printf("Gamma[%d][%d][%d] = %e\n", a, b, c, christof[a][b][c]); */
	/* 		} */
	/* 	} */
	/* } */
}

