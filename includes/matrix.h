#pragma once

#include "Geodesics.h"

class Matrix {
	public:
		double determinant3x3(double mat[3][3]);
		double determinant4x4(double mat[4][4]);
		void cofactor(double mat[NDIM][NDIM], double cofactorMat[NDIM][NDIM]);
		void transpose(double mat[NDIM][NDIM], double transposed[NDIM][NDIM]);
		int inverse_matrix(double mat[NDIM][NDIM], double inverse[NDIM][NDIM]);
		void check_inverse(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
		void print_matrix(const char *name, double matrix[NDIM][NDIM]);
		int inverse_3x3(double mat[3][3], double inv[3][3]);
		void check_inverse_3x3(double mat[3][3], double inv[3][3]);
		void print_matrix_3x3(const char* name, double mat[3][3]);
		void transpose3x3(double mat[3][3], double transposed[3][3]); 
		void cofactor3x3(double mat[3][3], double cofactorMat[3][3]);
		void check_cofactor3x3(double mat[3][3], double cofactorMat[3][3]);
		double determinant2x2(double mat[2][2]);
};
