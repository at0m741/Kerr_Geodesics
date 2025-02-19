#include <Geodesics.h>

double determinant2x2(double mat[2][2]) {
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

double determinant3x3(double mat[3][3]) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
         - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
         + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

double determinant4x4(double mat[NDIM][NDIM]) {
    double minor[3][3];
    double det = 0.0;
    for (int i = 0; i < NDIM; i++) {
        int subi = 0; 
        for (int j = 1; j < NDIM; j++) {
            int subj = 0;
            for (int k = 0; k < NDIM; k++) {
                if (k == i)
                    continue;
                minor[subi][subj] = mat[j][k];
                subj++;
            }
            subi++;
        }
        det += pow(-1, i) * mat[0][i] * determinant3x3(minor);
    }
    return det;
}

void cofactor(double mat[NDIM][NDIM], double cofactorMat[NDIM][NDIM]) {
    double minor[3][3];
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            int subi = 0;
            for (int x = 0; x < NDIM; x++) {
                if (x == i) continue;
                int subj = 0;
                for (int y = 0; y < NDIM; y++) {
                    if (y == j) continue;
                    minor[subi][subj] = mat[x][y];
                    subj++;
                }
                subi++;
            }
            cofactorMat[i][j] = pow(-1, i + j) * determinant3x3(minor);
        }
    }
}

void cofactor3x3(double mat[3][3], double cofactorMat[3][3]) {
    double minor[2][2];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int subi = 0;
            for (int x = 0; x < 3; x++) {
                if (x == i) continue; 
                int subj = 0;
                for (int y = 0; y < 3; y++) {
                    if (y == j) continue; 
                    minor[subi][subj] = mat[x][y];
                    subj++;
                }
                subi++;
            }
            cofactorMat[i][j] = pow(-1, i + j) * determinant2x2(minor);
        }
    }
}


void transpose(double mat[NDIM][NDIM], double transposed[NDIM][NDIM]) {
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            transposed[j][i] = mat[i][j];
        }
    }
}

void transpose3x3(double mat[3][3], double transposed[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            transposed[j][i] = mat[i][j];
        }
    }
}

int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]) {
    double det = determinant4x4(mat);
    double cofactorMat[NDIM][NDIM];
    cofactor(mat, cofactorMat);
    double adjugate[NDIM][NDIM];
    transpose(cofactorMat, adjugate);

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            inverse[i][j] = adjugate[i][j] / det;
        }
    }

    return 1;
}

int inverse_3x3(double mat[3][3], double inv[3][3]) {
    double det = determinant3x3(mat);
    if (fabs(det) < 1e-14) {
        printf("Matrix is singular or nearly singular!\n");
        return 0; 
    }

    double cofactorMat[3][3], adjugate[3][3];
    cofactor3x3(mat, cofactorMat);
    transpose3x3(cofactorMat, adjugate);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inv[i][j] = adjugate[i][j] / det;
        }
    }
    return 1;
}

void check_inverse_3x3(double mat[3][3], double inv[3][3]) {
    double product[3][3] = {{0.0}};
    double identity[3][3] = {{1.0, 0.0, 0.0},
                             {0.0, 1.0, 0.0},
                             {0.0, 0.0, 1.0}};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            product[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                product[i][j] += mat[i][k] * inv[k][j];
            }
        }
    }

    double TOL = 1e-10;  
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (fabs(product[i][j] - identity[i][j]) > TOL) {
                printf("Inversion 3x3 check failed at (%d, %d): %f != %f\n",
                       i, j, product[i][j], identity[i][j]);
            }
        }
    }
}


void check_inverse(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]) {
    double identity[NDIM][NDIM] = {0};
    for (int i = 0; i < NDIM; i++) {
        identity[i][i] = 1.0;
    }

    double product[NDIM][NDIM] = {0};

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            for (int k = 0; k < NDIM; k++) {
                product[i][j] += gcov[i][k] * gcon[k][j];
            }
        }
    }

    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            if (fabs(product[i][j] - identity[i][j]) > TOLERANCE) {
                printf("Matrix inversion check failed at element (%d, %d): %f\n", i, j, product[i][j]);
            }
        }
    }
}

void print_matrix(const char* name, double mat[NDIM][NDIM]) {
    printf("%s =\n", name);
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            printf("%.12f ", mat[i][j]);
        }
        printf("\n");
    }
}

void print_matrix_3x3(const char* name, double mat[3][3]) {
	printf("%s =\n", name);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%.12f ", mat[i][j]);
		}
		printf("\n");
	}
}
