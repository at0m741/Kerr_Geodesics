#pragma once

#include "Geodesics.h"

class Tensor {
	public:
		double gamma[NDIM][NDIM][NDIM];
		double Gamma_plus_h[NDIM][NDIM][NDIM];
		double Gamma_minus_h[NDIM][NDIM][NDIM];
		double Gamma_plus_half_h[NDIM][NDIM][NDIM][NDIM];
		double Gamma_minus_half_h[NDIM][NDIM][NDIM][NDIM];
};
