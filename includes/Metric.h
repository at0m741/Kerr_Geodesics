#pragma once	

#include "Geodesics.h"

class MetricTensor {
	public:
		double gcov[NDIM][NDIM];
		double gcon[NDIM][NDIM];
		double gcovK[NDIM][NDIM];
		double gconK[NDIM][NDIM];
		double gconMinkowski[NDIM][NDIM];
		double gcovMinkowski[NDIM][NDIM];

		MetricTensor() {
			memset(gcov, 0, sizeof(double) * NDIM * NDIM);
			memset(gcon, 0, sizeof(double) * NDIM * NDIM);
			memset(gcovK, 0, sizeof(double) * NDIM * NDIM);
			memset(gconK, 0, sizeof(double) * NDIM * NDIM);
			memset(gcovMinkowski, 0, sizeof(double) * NDIM * NDIM);
			memset(gconMinkowski, 0, sizeof(double) * NDIM * NDIM);
		}

		~MetricTensor() {}
};
