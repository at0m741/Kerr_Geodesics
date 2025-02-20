#pragma once

#include "Geodesics.h"

#define DIM3 3
using Matrix4x4 = std::array<std::array<double, NDIM>, NDIM>;
using Matrix3x3 = std::array<std::array<double, DIM3>, DIM3>;
using Vector3   = std::array<double, DIM3>;
using Vector4 = std::array<double, NDIM>;
using Tensor3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Tensor4D  = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;

class Grid {
    public:
        void extract_3p1(const Matrix4x4& g,
                         const Matrix4x4& g_inv,  
                         double* alpha,
                         Vector3& beta_cov,
                         Vector3& beta_con,
                         Matrix3x3& gamma,
                         Matrix3x3& gamma_inv);

        void calculeBeta(const Vector3& X, Vector3& beta_cov);
        void calculate_dbeta(const Vector3& X, Matrix3x3& dbeta);
        void compute_ricci_3d(const Vector3& X, const Tensor3D& Gamma3, Matrix3x3& R3);
        void print_ricci_tensor(const Matrix3x3& R3);
        double compute_K(const Matrix3x3& gamma_inv, const Matrix3x3& K);
        double compute_Kij_Kij(const Matrix3x3& gamma_inv, const Matrix3x3& K);
		void compute_extrinsic_curvature_stationary_3D(
				const Vector3& X,
				double alpha,
				const Vector3& beta_cov,
				std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>& Gamma3,
				const Matrix3x3& dbeta,
				Matrix3x3& K);
		double compute_hamiltonian_constraint(const Matrix3x3& gamma_inv, const Matrix3x3& K, const Matrix3x3& Ricci);
		void calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3);
};
