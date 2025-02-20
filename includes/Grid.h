#pragma once

#include "Geodesics.h"

class Grid {
	public:
		using Vector3   = std::array<double, 3>;
		void extract_3p1(const Matrix4x4& g,
                 const Matrix4x4& /*g_inv*/,  
                 double* alpha,
                 Vector3& beta_cov,
                 Vector3& beta_con,
                 Matrix3x3& gamma,
                 Matrix3x3& gamma_inv) ;

};
