#pragma once 

#include "Geodesics.h"

class Connexion {
	public:
		void calculate_christoffel(double X[NDIM], double h, \
							double gamma[NDIM][NDIM][NDIM],
							double g[NDIM][NDIM],
							double g_inv[NDIM][NDIM], const char *metric); 
		void check_symmetry_christoffel(double gamma[NDIM][NDIM][NDIM]); 
		void print_christoffel(double Gamma[NDIM][NDIM][NDIM]);
		void print_christoffel_matrix(double gamma[NDIM][NDIM][NDIM]);
};
