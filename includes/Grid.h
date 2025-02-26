#pragma once

#include <Geodesics.h>

#define DIM3 3
#define DX 1e-6
#define DY 1e-6
#define DZ 1e-6
#define NX 128
#define NY 128
#define NZ 128
#define GHOST 2  
#define NX_TOTAL (NX + 2*GHOST) 
#define NY_TOTAL (NY + 2*GHOST)
#define NZ_TOTAL (NZ + 2*GHOST)

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
			double beta[3];
			Matrix3x3 dgamma_dt;   
			Matrix3x3 dK_dt;      

			Vector3 dalpha_dx;    
			Tensor3D dK_dx;       
			Tensor3D dGamma3_dx;  

			double alphaStage[4];
			double betaStage[4][3];
			double rho;          
			Vector3 S_i;         
			Matrix3x3 S_ij;     
			double S; 
			double dgt[3][3]; 
			double dKt[3][3];  
			double gamma0[3][3];
			double K0[3][3];
			double alpha0;
			double beta0[3];
			double gammaStage[4][3][3];  
			double KStage[4][3][3];           
		};



		void evolve(double dt, int nSteps);
		void extract_3p1(const Matrix4x4& g,
				const Matrix4x4& g_inv,  
				double* alpha,
				Vector3& beta_cov,
				Vector3& beta_con,
				Matrix3x3& gamma,
				Matrix3x3& gamma_inv);
		void initializeData(); 
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
				const Matrix3x3& gamma, Matrix3x3 gamma_inv) ;
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
		void compute_ricci_3d(
				const Vector3& X,       
				const Tensor3D& Gamma3, 
				Matrix3x3& R3);
		void compute_ricci_3d_grid(
				std::vector<std::vector<Cell2D>>& grid,
				int Nr, int Ntheta,
				double dr, double dtheta,
				double r_min, double theta_min,
				double delta);
		std::vector<std::vector<Cell2D>> grid; 
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
				const std::array<Christoffel3D, 3>& Gamma_plus_h,
				const std::array<Christoffel3D, 3>& Gamma_minus_h,
				const std::array<Christoffel3D, 3>& Gamma_plus_half_h,
				const std::array<Christoffel3D, 3>& Gamma_minus_half_h,
				Riemann3D& Riemann,
				double h,
				double scale); 
		void calculate_ricci_3d_from_riemann(const Riemann3D& Riemann, Matrix3x3& Ricci);
		bool verify_riemann_symmetries(const Riemann3D &Riemann);
		void calculate_riemann_4d_from_3d(
				const Riemann3D &Riemann3,
				const Matrix3x3 &K,     
				Riemann3D &Riemann4) ;
		void compute_time_derivatives(int i, int j, int k);
		void allocateGlobalGrid();
		void initializeData_Minkowski();
};



extern std::vector<std::vector<std::vector<Grid::Cell2D>>> globalGrid;
double partialX_alpha(int i, int j, int k);
double partialY_alpha(int i, int j, int k);
double partialZ_alpha(int i, int j, int k);
double partialXX_alpha(int i, int j, int k);
double partialYY_alpha(int i, int j, int k);
double partialZZ_alpha(int i, int j, int k);
double partialX_gamma(int i, int j, int k, int a, int b);
double partialY_gamma(int i, int j, int k, int a, int b);
double partialZ_gamma(int i, int j, int k, int a, int b);
double partialX_betacomp(int i, int j, int k, int comp);
double partialY_betacomp(int i, int j, int k, int comp);
double partialZ_betacomp(int i, int j, int k, int comp);
double partialXY_alpha(int i, int j, int k);
double partialXZ_alpha(int i, int j, int k);
double partialYZ_alpha(int i, int j, int k);
double second_partial_alpha(int i, int j, int k, int a, int b);
void compute_gauge_derivatives(int i, int j, int k, double &d_alpha_dt, double d_beta_dt[3]);
void compute_christoffel_3D(int i, int j, int k, double christof[3][3][3]);
bool invert_3x3(const double m[3][3], double inv[3][3]);
void compute_ricci_3D(int i, int j, int k, double Ricci[3][3]);
