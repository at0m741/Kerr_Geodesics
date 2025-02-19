#pragma once	

#include "Geodesics.h"

class Metric {
	public:
		double gcov[NDIM][NDIM];
		double gcon[NDIM][NDIM];
		double gcovK[NDIM][NDIM];
		double gconK[NDIM][NDIM];
		double gconMinkowski[NDIM][NDIM];
		double gcovMinkowski[NDIM][NDIM];

		Metric() {
			memset(gcov, 0, sizeof(double) * NDIM * NDIM);
			memset(gcon, 0, sizeof(double) * NDIM * NDIM);
			memset(gcovK, 0, sizeof(double) * NDIM * NDIM);
			memset(gconK, 0, sizeof(double) * NDIM * NDIM);
			memset(gcovMinkowski, 0, sizeof(double) * NDIM * NDIM);
			memset(gconMinkowski, 0, sizeof(double) * NDIM * NDIM);
		}
		
		~Metric() {}


		void calculate_metric(double x[4], double g[4][4], double g_inv[4][4]);
		void verify_metric(double g[4][4], double g_inv[4][4]);
		void calculate_metric_kds(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]); 
		void calculate_metric_kerr_newman(double x[NDIM], double g[NDIM][NDIM], double g_inv[NDIM][NDIM]);
};
