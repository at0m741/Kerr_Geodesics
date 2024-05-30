#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../headers/geodesics.h"


typedef struct {
    double rho;
    double vx;
    double vy;
    double P;
    double flux[4]; 
} __attribute__((aligned(32))) hydro;

#pragma omp declare simd
void christoffel(double g[4][4], double christoffel[4][4][4])
{
	#ifdef _OPENMP
		printf("OpenMP is supported and used\n");
	#else
		printf("OpenMP is not supported\n");
	#endif

    double (*g_aligned)[4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4]));
    double (*christoffel_aligned)[4][4] = aligned_alloc(ALIGNMENT, sizeof(double[4][4][4]));
    memcpy(g_aligned, g, sizeof(double[4][4]));
    memcpy(christoffel_aligned, christoffel, sizeof(double[4][4][4]));


    #pragma omp parallel for collapse(3)
    #pragma vector aligned
    for (int mu = 0; mu < 4; mu++) {
        for (int beta = 0; beta < 4; beta++) {
            for (int nu = 0; nu < 4; nu++) {
                double sum = 0;
                #pragma omp simd reduction(+:sum) aligned(g_aligned:ALIGNMENT)
                #pragma vector aligned
                for (int sigma = 0; sigma < 4; sigma++) {
                    sum += 0.5 * (g_aligned[mu][sigma] * (g_aligned[sigma][beta] + g_aligned[beta][sigma] - g_aligned[beta][nu]));
                }
                christoffel_aligned[mu][beta][nu] = sum;
            }
        }
    }

    memcpy(christoffel, christoffel_aligned, sizeof(double[4][4][4]));
    free(g_aligned);
    free(christoffel_aligned);
}

double random_perturbation(double min, double max) {
    srand(time(NULL));
    return min + ((double)rand() / RAND_MAX) * (max - min);
}

double max_wave_speed(double *ul, double *ur, double gamma) {
    double sl = sqrt(gamma * ul[2] / ul[0]);
    double sr = sqrt(gamma * ur[2] / ur[0]);
    return fmax(sl, sr);
}

double min_wave_speed(double *ul, double *ur, double gamma) {
    double sl = sqrt(gamma * ul[2] / ul[0]);
    double sr = sqrt(gamma * ur[2] / ur[0]);
    return fmin(sl, sr);
}

double stress_energie_tensor(double *u, double *T, double gamma, double christo[4][4][4]) {
    double v_rel[4];
    #pragma omp simd
    for (int i = 0; i < 4; ++i) {
        v_rel[i] = u[i];
    }
    #pragma omp simd
    for (int i = 1; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
            v_rel[i] += -christo[i][0][j] * u[0] * u[j];
        }
    }

    T[0] = u[0] * v_rel[1];
    T[1] = u[0] * v_rel[1] * v_rel[1] + u[3];
    T[2] = u[0] * v_rel[2];
    T[3] = v_rel[1] * (u[3] + u[1]);

    return T[0] + T[1] + T[2] + T[3];
}

#pragma omp declare simd
void hlle_solver(double *ul, double *ur, double *flux, double gamma, double christo[4][4][4]) {

    double T[4];
    double Tl[4] = {0};
    double ul_v_rel[4];
    double ur_v_rel[4];
    stress_energie_tensor(ul_v_rel, Tl, gamma, christo);
    #pragma omp simd
    for (int i = 0; i < 4; ++i) {
        ul_v_rel[i] = ul[i];
        ur_v_rel[i] = ur[i];
    }
    #pragma omp simd
    for (int i = 1; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
            ul_v_rel[i] += -christo[i][0][j] * ul[0] * ul[j];
            ur_v_rel[i] += -christo[i][0][j] * ur[0] * ur[j];
        }
    }

    double pL = (gamma - 1.0) * (ul[3] - 0.5 * ul[0] * (ul_v_rel[1] * ul_v_rel[1] + ul_v_rel[2] * ul_v_rel[2]));
    double pR = (gamma - 1.0) * (ur[3] - 0.5 * ur[0] * (ur_v_rel[1] * ur_v_rel[1] + ur_v_rel[2] * ur_v_rel[2]));

    double sL = fmin(ul_v_rel[1] - sqrt(gamma * pL / ul[0]), ur_v_rel[1] - sqrt(gamma * pR / ur[0]));
    double sR = fmax(ul_v_rel[1] + sqrt(gamma * pL / ul[0]), ur_v_rel[1] + sqrt(gamma * pR / ur[0]));
    if (sL >= 0) {
        flux[0] = ul[0] * ul_v_rel[1];
        flux[1] = ul[0] * ul_v_rel[1] * ul_v_rel[1] + pL;
        flux[2] = ul[0] * ul_v_rel[2];
        flux[3] = ul_v_rel[1] * (ul[3] + pL);
    } else if (sR <= 0) {
        flux[0] = ur[0] * ur_v_rel[1];
        flux[1] = ur[0] * ur_v_rel[1] * ur_v_rel[1] + pR;
        flux[2] = ur[0] * ur_v_rel[2];
        flux[3] = ur_v_rel[1] * (ur[3] + pR);
    } else {
        double invSL_SR = 1.0 / (sR - sL);
        flux[0] = (sR * ul[0] * ul_v_rel[1] - sL * ur[0] * ur_v_rel[1] + sL * sR * (ur[0] * ur_v_rel[1] - ul[0] * ul_v_rel[1])) * invSL_SR;
        flux[1] = (sR * (ul[0] * ul_v_rel[1] * ul_v_rel[1] + pL) - sL * (ur[0] * ur_v_rel[1] * ur_v_rel[1] + pR) + sL * sR * (ur[0] * ur_v_rel[1] - ul[0] * ul_v_rel[1])) * invSL_SR;
        flux[2] = (sR * ul[0] * ul_v_rel[2] - sL * ur[0] * ur_v_rel[2] + sL * sR * (ur[0] * ur_v_rel[2] - ul[0] * ul_v_rel[2])) * invSL_SR;
        flux[3] = (sR * ul_v_rel[1] * (ul[3] + pL) - sL * ur_v_rel[1] * (ur[3] + pR) + sL * sR * (ur[3] - ul[3])) * invSL_SR;
    }
}




void Riemann_solver(hydro *left, hydro *right, hydro *flux, double gamma, double christo[4][4][4]) {
    stress_energie_tensor(left->flux, left->flux, gamma, christo);
    double ul[4] = {left->rho, left->rho * left->vx, left->rho * left->vy, left->P / \
                    (gamma - 1.0) + 0.5 * left->rho * (left->vx * left->vx + left->vy * left->vy)};
    double ur[4] = {right->rho, right->rho * right->vx, right->rho * right->vy, right->P / \
                    (gamma - 1.0) + 0.5 * right->rho * (right->vx * right->vx + right->vy * right->vy)};

    hlle_solver(ul, ur, flux->flux, gamma, christo);
}

void initialize(hydro **grid, int nx, int ny, double gamma, double christo[4][4][4]) {
    stress_energie_tensor(grid[0][0].flux, grid[0][0].flux, gamma, christo);
    const double r_inner = -10.0;
    const double r_outer = 10.0;
    const double theta_inner = -1.0;
    const double theta_outer =  M_PI * 2.0;
    const double density_inner = 3.0;
    const double density_outer = 0.3;
    const double pressure_inner = 1.0;
    const double pressure_outer = 0.02;
    const double dx = 10.0 / nx;
    const double dy = 10.0 / ny;
    #pragma omp parallel for simd
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = i * dx;
            double y = j * dy;
            double Rs = 2 * G * M / c * c;
            double r = Rs + sqrt(x * x + y * y);
            double theta = atan2(y, x);

            if (r_inner <= r && r <= r_outer && theta_inner <= theta && theta <= theta_outer) {
                grid[i][j].rho = density_inner + cos(theta);
                grid[i][j].P = pressure_inner + sin(theta) + rand() / RAND_MAX;
            } else {
                grid[i][j].rho = density_outer;
                grid[i][j].P = pressure_outer;
            }

            grid[i][j].vx = 0.0;
            grid[i][j].vy = 0.0;
        }
    }
}
void momentum_flux(hydro **grid, hydro **flux_x, hydro **flux_y, int nx, int ny, double gamma, double christo[4][4][4]) {
    stress_energie_tensor(grid[0][0].flux, grid[0][0].flux, 1.4, christo);
    #pragma omp parallel for simd
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double pressure_term_x = grid[i][j].P / (grid[i][j].rho * (1.4 - 1.0)) * grid[i][j].vx;
            double pressure_term_y = grid[i][j].P / (grid[i][j].rho * (1.4 - 1.0)) * grid[i][j].vy;

            flux_x[i][j].flux[0] = grid[i][j].rho * grid[i][j].vx;
            flux_x[i][j].flux[1] = grid[i][j].rho * grid[i][j].vx * grid[i][j].vx + pressure_term_x;
            flux_x[i][j].flux[2] = grid[i][j].rho * grid[i][j].vx * grid[i][j].vy + pressure_term_x;
            flux_x[i][j].flux[3] = grid[i][j].vx * (pressure_term_x + 0.5 * (grid[i][j].vx * grid[i][j].vx + grid[i][j].vy * grid[i][j].vy));

            flux_y[i][j].flux[0] = grid[i][j].rho * grid[i][j].vy;
            flux_y[i][j].flux[1] = grid[i][j].rho * grid[i][j].vy * grid[i][j].vx;
            flux_y[i][j].flux[2] = grid[i][j].rho * grid[i][j].vy * grid[i][j].vy + pressure_term_y;
            flux_y[i][j].flux[3] = grid[i][j].vy * (pressure_term_y + 0.5 * (grid[i][j].vx * grid[i][j].vx + grid[i][j].vy * grid[i][j].vy));
        }
    }

    //conditions aux limites
    #pragma omp parallel for simd
    for (int j = 0; j < ny; ++j) {
        flux_x[nx-1][j].flux[0] = flux_x[0][j].flux[0];
        flux_x[nx-1][j].flux[1] = flux_x[0][j].flux[1];
        flux_x[nx-1][j].flux[2] = flux_x[0][j].flux[2];
        flux_x[nx-1][j].flux[3] = flux_x[0][j].flux[3];
    }
    #pragma omp parallel for simd
    for (int i = 0; i < nx; ++i) {
        flux_y[i][ny-1].flux[0] = flux_y[i][0].flux[0];
        flux_y[i][ny-1].flux[1] = flux_y[i][0].flux[1];
        flux_y[i][ny-1].flux[2] = flux_y[i][0].flux[2];
        flux_y[i][ny-1].flux[3] = flux_y[i][0].flux[3];
    }

}

void update(hydro **grid, hydro **flux_x, hydro **flux_y, int nx, int ny, double dt, double dx, double dy, double christo[4][4][4]) {
    momentum_flux(grid, flux_x, flux_y, nx, ny, 1.4, christo);
    #pragma omp parallel for simd
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            grid[i][j].rho -= dt / dx * (flux_x[i][j].flux[0] - flux_x[i-1][j].flux[0]) + dt / dy * (flux_y[i][j].flux[0] - flux_y[i][j-1].flux[0]);
            grid[i][j].vx -= dt / dx * (flux_x[i][j].flux[1] - flux_x[i-1][j].flux[1]) + dt / dy * (flux_y[i][j].flux[1] - flux_y[i][j-1].flux[1]);
            grid[i][j].vy -= dt / dx * (flux_x[i][j].flux[2] - flux_x[i-1][j].flux[2]) + dt / dy * (flux_y[i][j].flux[2] - flux_y[i][j-1].flux[2]);
            grid[i][j].P -= dt / dx * (flux_x[i][j].flux[3] - flux_x[i-1][j].flux[3]) + dt / dy * (flux_y[i][j].flux[3] - flux_y[i][j-1].flux[3]);
        }
    }
    //conditions aux limites
    #pragma omp parallel for simd
    for (int j = 0; j < ny; ++j) {
        grid[nx-1][j].rho = grid[0][j].rho;
        grid[nx-1][j].vx = grid[0][j].vx;
        grid[nx-1][j].vy = grid[0][j].vy;
        grid[nx-1][j].P = grid[0][j].P;
    }
    #pragma omp parallel for simd shared(grid)
    for (int i = 0; i < nx; ++i) {
        grid[i][ny-1].rho = grid[i][0].rho;
        grid[i][ny-1].vx = grid[i][0].vx;
        grid[i][ny-1].vy = grid[i][0].vy;
        grid[i][ny-1].P = grid[i][0].P;
    }

   // momentum_flux(grid, flux_x, flux_y, nx, ny);
}


// int main() {
//     double gamma = 0.4;
//     int nx = 256;
//     int ny = 256;
//     double Lx = 1.0;
//     double Ly = 1.0;
//     double dx = Lx;
//     double dy = Ly;
//     double dt = 0.0001;
//     int nt = 1; // Nombre de pas de temps

//     hydro **grid = (hydro**)malloc(nx * sizeof(hydro*));
//     hydro **flux_x = (hydro**)malloc(nx * sizeof(hydro*));
//     hydro **flux_y = (hydro**)malloc(nx * sizeof(hydro*));
//     for (int i = 0; i < nx; ++i) {
//         grid[i] = (hydro*)malloc(ny * sizeof(hydro));
//         flux_x[i] = (hydro*)malloc(ny * sizeof(hydro));
//         flux_y[i] = (hydro*)malloc(ny * sizeof(hydro));
//     }

// 	ldouble_a Q = 1.0;
// 	ldouble_a Rs = 2 * G * M / powf(c, 2);
//     double v[4] = {-7.0, 1.01, 1.0, 27.0};        // Vitesse initiale (dr/dλ, dθ/dλ, dφ/dλ, dt/dλ)
//     double x[4] = {Rs, M_PI / 2, 60 * M_PI, Rs}; // Position initiale (r, θ, φ, t)
//     ldouble_a r = sqrt(powf(x[1], 2) + powf(a, 2) * powf(cos(x[2]), 2));
// 	ldouble_a rho2_kn = powf(r, 2) + powf(a, 2) * powf(cos(x[1]), 2);
// 	ldouble_a delta_kn = powf(r, 2) - Rs * r + powf(a, 2) + powf(Q, 2);
// 	ldouble_a Sigma_kn = powf((powf(r, 2) + powf(a, 2)), 2) - powf(a, 2) * delta_kn *\
// 					     powf(sin(x[1]), 2);
// 	double g_kerr_newman[4][4] = {0};
//     double christoffel_sym[4][4][4] = {0};

// 	g_kerr_newman[0][0] = -(1 - (Rs * r - powf(Q, 2)) / rho2_kn);
// 	g_kerr_newman[0][3] = -(Rs * r - powf(Q, 2)) * a * powf(sin(x[1]), 2) / rho2_kn;
// 	g_kerr_newman[1][1] = rho2_kn / delta_kn;
// 	g_kerr_newman[2][2] = rho2_kn;
// 	g_kerr_newman[3][0] = g_kerr_newman[0][3];
// 	g_kerr_newman[3][3] = (powf(r, 2) + powf(a, 2) + (Rs * r - powf(Q, 2)) * powf(a, 2) * \
// 						   powf(sin(x[1]), 2) / rho2_kn) * powf(sin(x[1]), 2);

//     double riemann_tensor[4][4][4][4] = {0};
//     christoffel(g_kerr_newman, christoffel_sym);
//     initialize(grid, nx, ny, gamma, christoffel_sym);
//     #pragma omp parallel for simd
//     for (int t = 0; t < nt; ++t) {
//         for (int i = 0; i < nx - 1; ++i) {
//             for (int j = 0; j < ny - 1; ++j) {
//                 Riemann_solver(&grid[i][j], &grid[i+1][j], &flux_x[i][j], gamma, christoffel_sym);
//                 Riemann_solver(&grid[i][j], &grid[i][j+1], &flux_y[i][j], gamma, christoffel_sym);
//             }
//         }
//         update(grid, flux_x, flux_y, nx, ny, dt, dx, dy, christoffel_sym);
//     }
//     //write vts file for paraview
//     FILE *file = fopen("hydro.vts", "w");
//     fprintf(file, "<?xml version=\"1.0\"?>\n");
//     fprintf(file, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
//     fprintf(file, "<StructuredGrid WholeExtent=\"0 %d 0 %d 0 0\">\n", nx-1, ny-1);
//     fprintf(file, "<Piece Extent=\"0 %d 0 %d 0 0\">\n", nx-1, ny-1);
//     fprintf(file, "<Points>\n");
//     fprintf(file, "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             fprintf(file, "%f %f %f\n", i * dx, j * dy, 0.0);
//         }
//     }
//     fprintf(file, "</DataArray>\n");
//     fprintf(file, "</Points>\n");
//     fprintf(file, "<PointData Scalars=\"rho\">\n");
//     fprintf(file, "<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n");
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             fprintf(file, "%f ", grid[i][j].rho);
//         }
//         fprintf(file, "\n");
//     }
//     fprintf(file, "</DataArray>\n");
//     fprintf(file, "<DataArray type=\"Float64\" Name=\"vx\" format=\"ascii\">\n");
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             fprintf(file, "%f ", grid[i][j].vx);
//         }
//         fprintf(file, "\n");
//     }
//     fprintf(file, "</DataArray>\n");
//     fprintf(file, "<DataArray type=\"Float64\" Name=\"vy\" format=\"ascii\">\n");
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             fprintf(file, "%f ", grid[i][j].vy);
//         }
//         fprintf(file, "\n");
//     }
//     fprintf(file, "</DataArray>\n");
//     fprintf(file, "<DataArray type=\"Float64\" Name=\"P\" format=\"ascii\">\n");
//     for (int i = 0; i < nx; ++i) {
//         for (int j = 0; j < ny; ++j) {
//             fprintf(file, "%f ", grid[i][j].P);
//         }
//         fprintf(file, "\n");
//     }

//     fprintf(file, "</DataArray>\n");
//     fprintf(file, "</PointData>\n");
//     fprintf(file, "</Piece>\n");
//     fprintf(file, "</StructuredGrid>\n");
//     fprintf(file, "</VTKFile>\n");
//     fclose(file);



//     for (int i = 0; i < nx; ++i) {
//         free(grid[i]);
//         free(flux_x[i]);
//         free(flux_y[i]);
//     }
//     free(grid);
//     free(flux_x);
//     free(flux_y);

//     return 0;
// }