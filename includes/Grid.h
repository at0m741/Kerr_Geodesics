#pragma once

#include <Geodesics.h>

#define DIM3 3
using Matrix4x4 = std::array<std::array<double, NDIM>, NDIM>;
using Matrix3x3 = std::array<std::array<double, DIM3>, DIM3>;
using Vector3   = std::array<double, DIM3>;
using Vector4 = std::array<double, NDIM>;
using Tensor3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Tensor4D  = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;




class Grid {
    public:
		

		struct Cell2D {
			Matrix3x3 gamma;       // gamma_ij (métrique spatiale)
			Matrix3x3 gamma_inv;   // gamma^{ij} (inverse de gamma_ij)
			Tensor3D Gamma3;       // Gamma^k_{ij} (symboles de Christoffel en 3D)
			Matrix3x3 K;           // K_ij (courbure extrinsèque)
			Matrix3x3 Ricci;       // Ricci_ij (tenseur de Ricci)
			double H;              // Contrainte Hamiltonienne

			double alpha;          // Lapse function α
			Vector3 beta_cov;      // Shift vector β_i (composantes covariantes)
			Vector3 beta_con;      // Shift vector β^i (composantes contravariantes)

			Matrix3x3 dgamma_dt;   // ∂γ_ij / ∂t (évolution temporelle de gamma_ij)
			Matrix3x3 dK_dt;       // ∂K_ij / ∂t (évolution temporelle de K_ij)

			Vector3 dalpha_dx;     // Gradient de α : ∂α / ∂x^i
			Tensor3D dK_dx;        // Gradient de K_ij : ∂K_ij / ∂x^k
			Tensor3D dGamma3_dx;   // Gradient des symboles de Christoffel : ∂Γ^k_{ij} / ∂x^l

			double rho;            // Densité d'énergie de la matière
			Vector3 S_i;           // Flux de moment S_i
			Matrix3x3 S_ij;        // Tenseur de contraintes S_ij
			double S;              // Trace de S_ij : S = gamma^ij S_ij
		};

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
		void calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3);void calculate_christoffel_3D_grid(
				std::vector<std::vector<Cell2D>>& grid,
				int Nx, int Ny,
				double dx, double dy);
};
