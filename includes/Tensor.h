#pragma once

#include "Geodesics.h"

class Tensor {
	public:
		double gamma[NDIM][NDIM][NDIM];
		double Gamma_plus_h[NDIM][NDIM][NDIM];
		double Gamma_minus_h[NDIM][NDIM][NDIM];
		double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM];
		double Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM];

		void calculate_Gamma_at_offset(double X[NDIM], int direction, 
				double offset, double delta,
				double gcov[NDIM][NDIM], 
				double gcon[NDIM][NDIM], 
				double Gamma_slice[NDIM][NDIM][NDIM], 
				const char *metric_type);
		void calculate_riemann(double Gamma[NDIM][NDIM][NDIM], 
				double Gamma_plus_h[NDIM][NDIM][NDIM][NDIM], 
				double Gamma_minus_h[NDIM][NDIM][NDIM][NDIM], 
				double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM], 
				double Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM],
				double Riemann[NDIM][NDIM][NDIM][NDIM], 
				double h);
		

 
		void initialize_riemann_tensor(double R[NDIM][NDIM][NDIM][NDIM]);
		void print_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM]);
		void check_riemann_symmetries(double Riemann[NDIM][NDIM][NDIM][NDIM], double tolerance);
		void contract_riemann(double Riemann[NDIM][NDIM][NDIM][NDIM], double Ricci[NDIM][NDIM], double g_inv[NDIM][NDIM]);

		Tensor() {
			memset(gamma, 0, sizeof(double) * NDIM * NDIM * NDIM);
			memset(Gamma_plus_h, 0, sizeof(double) * NDIM * NDIM * NDIM);
			memset(Gamma_minus_h, 0, sizeof(double) * NDIM * NDIM * NDIM);
			memset(Gamma_plus_half_h, 0, sizeof(double) * NDIM * NDIM * NDIM * NDIM);
			memset(Gamma_minus_half_h, 0, sizeof(double) * NDIM * NDIM * NDIM * NDIM);
		}
};


