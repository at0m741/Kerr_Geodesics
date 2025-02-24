#pragma once

#include <Geodesics.h>

#define DIM3 3
using Matrix4x4 = std::array<std::array<double, NDIM>, NDIM>;
using Matrix3x3 = std::array<std::array<double, DIM3>, DIM3>;
using Vector3   = std::array<double, DIM3>;
using Vector4 = std::array<double, NDIM>;
using Tensor3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Tensor4D  = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;
using Christoffel3D = std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>;
using Riemann3D = std::array<std::array<std::array<std::array<double, DIM3>, DIM3>, DIM3>, DIM3>;



class Grid {
    public:
		
		
		struct Cell2D {
			Matrix3x3 gamma;     
			Matrix3x3 gamma_inv; 
			Tensor3D Gamma3;    
			Matrix3x3 K;        
			Matrix3x3 Ricci;     
			double H;    

			double alpha;        
			Vector3 beta_cov;     
			Vector3 beta_con;    

			Matrix3x3 dgamma_dt;   
			Matrix3x3 dK_dt;      

			Vector3 dalpha_dx;    
			Tensor3D dK_dx;       
			Tensor3D dGamma3_dx;  

			double rho;          
			Vector3 S_i;         
			Matrix3x3 S_ij;     
			double S;            
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
		void compute_ricci_3d(
				Grid& grid_obj,  
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
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
		void calculate_christoffel_3D(const Vector3& X, Tensor3D& Gamma3, 
				const Matrix3x3& gamma, const Matrix3x3& gamma_inv);
		double compute_hamiltonian_constraint(const Matrix3x3& gamma_inv, const Matrix3x3& K, const Matrix3x3& Ricci);
		void calculate_christoffel_3D_grid(
				std::vector<std::vector<Grid::Cell2D>>& grid,
				int Nx, int Ny,
				double dr, double dtheta,
				double r_min,
				double theta_min);
		void save_christoffel_symbols(const std::vector<std::vector<Cell2D>>& grid,
				int Nx, int Ny,
				const std::string &filename);
		void compute_ricci_3d_grid(
				std::vector<std::vector<Cell2D>>& grid,
				int Nr, int Ntheta,
				double dr, double dtheta,
				double r_min, double theta_min,
				double delta);
		std::vector<std::vector<Cell2D>> grid; // Stocke la grille
		int Nr, Ntheta;
		double dr, dtheta;
		void save_extrinsic_curvature(const std::vector<std::vector<Cell2D>>& grid,
				int Nx, int Ny,
				const std::string &filename);
		void evolve_Kij(double dt);
		Matrix3x3 compute_second_derivative_alpha(int i, int j);

		Matrix3x3 compute_beta_gradient(int i, int j);
		double compute_KijKij_component(const Matrix3x3& gamma_inv, const Matrix3x3& K, int a, int b);
		void initialize_grid(int Nr, int Ntheta, double r_min, double r_max, double theta_min, double theta_max);
		void export_vtk(const std::string& filename);
		double richardson_derivative_ricci(
				const Tensor3D &Gamma_plus_h, 
				const Tensor3D &Gamma_minus_h,
				const Tensor3D &Gamma_plus_half_h,
				const Tensor3D &Gamma_minus_half_h,
				int mu, int nu, int sigma, 
				double h) ; 

		void calculate_riemann_3d(
				const Christoffel3D& Gamma, 
				const Christoffel3D& Gamma_plus_h, 
				const Christoffel3D& Gamma_minus_h, 
				const Christoffel3D& Gamma_plus_half_h,
				const Christoffel3D& Gamma_minus_half_h,
				Riemann3D& Riemann,
				double h,
				double scale);
		void calculate_ricci_3d_from_riemann(const Riemann3D& Riemann, Matrix3x3& Ricci);
};
