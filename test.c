#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NDIM 4
#define TT 0
#define RR 1
#define TH 2
#define PHI 3
#define ALIGNMENT 64

int Nr = 200;
int Nth = 10;
int Nphi = 10;

double r_min = 1.0;
double r_max = 2.0;
double theta_min = 0.0;
double theta_max = M_PI / 2.0;
double phi_min = 0.0;
double phi_max = 2.0 * M_PI;

#define SMALL 1e-10
double a = 0.9;
double viscosity = 0.1;

void Boyer_lindquist_coord(double *X, double *r, double *th)
{
    double R0 = 1.0, hslope = 0.3;
    *r = exp(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    if (fabs(*th) < SMALL) {
        if ((*th) >= 0)
            *th = SMALL;
        if ((*th) < 0)
            *th = -SMALL;
    }
    if (fabs(M_PI - (*th)) < SMALL) {
        if ((*th) >= M_PI)
            *th = M_PI + SMALL;
        if ((*th) < M_PI)
            *th = M_PI - SMALL;
    }
    return;
}

void gcov(double *X, double g_cov[][NDIM])
{
    double r, th;
    Boyer_lindquist_coord(X, &r, &th);
    double sth, cth, s2, rho2;
    cth = cos(th);
    sth = sin(th);
    s2 = sth * sth;
    rho2 = r * r + a * a * cth * cth;

    g_cov[TT][TT] = (-1. + 2. * r / rho2);
    g_cov[TT][RR] = (2. * r / rho2);
    g_cov[TT][PHI] = (-2. * a * r * s2 / rho2);
    g_cov[RR][TT] = g_cov[TT][RR];
    g_cov[RR][RR] = (1. + 2. * r / rho2);
    g_cov[RR][PHI] = (-a * s2 * (1. + 2. * r / rho2));
    g_cov[TH][TH] = rho2;
    g_cov[PHI][TT] = g_cov[TT][PHI];
    g_cov[PHI][RR] = g_cov[RR][PHI];
    g_cov[PHI][PHI] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            g_cov[mu][nu] *= -1.0;
        }
    }
}

void calculate_Tmunu(double *density, double *pressure, double *velocity, double g_cov[][NDIM], \
                     double Tmunu[][NDIM], int idx, int Nr, int Nth, int Nphi)
{
    double rho = density[idx];
    double p = pressure[idx];
    double u[NDIM];
    for (int d = 0; d < NDIM; d++) {
        u[d] = velocity[d * Nr * Nth * Nphi + idx];
    }
    #pragma omp simd aligned(Tmunu, g_cov: ALIGNMENT)
    for (int mu = 0; mu < NDIM; mu++) {
        for (int nu = 0; nu < NDIM; nu++) {
            Tmunu[mu][nu] = (rho + p) * u[mu] * u[nu] + p * g_cov[mu][nu];
        }
    }
}

void compute_fluxes(double *density, double *pressure, double *velocity, double *flux_density, double *flux_pressure, \
                    double *flux_velocity, double *r_grid, double *theta_grid, double *phi_grid, int Nr, int Nth, int Nphi)
{
    double g_cov[NDIM][NDIM];
    double Tmunu[NDIM][NDIM];
    double X[NDIM];

    #pragma omp parallel for private(g_cov, Tmunu, X)
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            #pragma omp simd aligned(flux_density, flux_pressure, flux_velocity: ALIGNMENT)
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                X[1] = r_grid[i];
                X[2] = theta_grid[j];
                X[3] = phi_grid[k];
                gcov(X, g_cov);
                calculate_Tmunu(density, pressure, velocity, g_cov, Tmunu, idx, Nr, Nth, Nphi);
                flux_density[idx] = Tmunu[TT][RR];
                flux_pressure[idx] = Tmunu[TT][RR];
                for (int d = 0; d < NDIM; d++) {
                    flux_velocity[d * Nr * Nth * Nphi + idx] = Tmunu[d][RR];
                }
            }
        }
    }
}

void apply_viscosity(double *velocity, int Nr, int Nth, int Nphi, double viscosity)
{
    #pragma omp parallel for
    for (int i = 1; i < Nr - 1; i++) {
        for (int j = 1; j < Nth - 1; j++) {
            for (int k = 1; k < Nphi - 1; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                #pragma omp simd aligned(velocity: ALIGNMENT)
                for (int d = 0; d < NDIM; d++) {
                    int idx_rp = (i + 1) * Nth * Nphi + j * Nphi + k;
                    int idx_rm = (i - 1) * Nth * Nphi + j * Nphi + k;
                    int idx_tp = i * Nth * Nphi + (j + 1) * Nphi + k;
                    int idx_tm = i * Nth * Nphi + (j - 1) * Nphi + k;
                    int idx_pp = i * Nth * Nphi + j * Nphi + (k + 1);
                    int idx_pm = i * Nth * Nphi + j * Nphi + (k - 1);
                    
                    double laplacian = (velocity[d * Nr * Nth * Nphi + idx_rp] +
                                        velocity[d * Nr * Nth * Nphi + idx_rm] +
                                        velocity[d * Nr * Nth * Nphi + idx_tp] +
                                        velocity[d * Nr * Nth * Nphi + idx_tm] +
                                        velocity[d * Nr * Nth * Nphi + idx_pp] +
                                        velocity[d * Nr * Nth * Nphi + idx_pm] -
                                        6.0 * velocity[d * Nr * Nth * Nphi + idx]);

                    velocity[d * Nr * Nth * Nphi + idx] += viscosity * laplacian;
                }
            }
        }
    }
}

void update_variables(double *density, double *pressure, double *velocity, double *flux_density, double *flux_pressure, \
                      double *flux_velocity, double dt, int Nr, int Nth, int Nphi)
{
    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            #pragma omp simd
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                density[idx] += dt * flux_density[idx];
                pressure[idx] += dt * flux_pressure[idx];
                #pragma omp simd
                for (int d = 0; d < NDIM; d++) {
                    velocity[d * Nr * Nth * Nphi + idx] += dt * flux_velocity[d * Nr * Nth * Nphi + idx];
                }
            }
        }
    }
}

void evolve_fluid(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid,\
                  double *phi_grid, int Nr, int Nth, int Nphi, double dt, int num_steps)
{
    double *flux_density = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *flux_pressure = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *flux_velocity = (double *)malloc(NDIM * Nr * Nth * Nphi * sizeof(double));

    for (int n = 0; n < num_steps; n++) {
        printf("Step %d\n", n);
        
        compute_fluxes(density, pressure, velocity, flux_density, flux_pressure, flux_velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi);
        apply_viscosity(velocity, Nr, Nth, Nphi, viscosity);
        printf("Before update - Density at [50, 25, 25]: %f\n", density[50 * Nth * Nphi + 25 * Nphi + 25]);
        update_variables(density, pressure, velocity, flux_density, flux_pressure, flux_velocity, dt, Nr, Nth, Nphi);
        printf("After update - Density at [50, 25, 25]: %f\n", density[50 * Nth * Nphi + 25 * Nphi + 25]);
    }

    free(flux_density);
    free(flux_pressure);
    free(flux_velocity);
}

void write_binary_file(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid,\
                       double *phi_grid, int Nr, int Nth, int Nphi, int step)
{
    char filename[256];
    sprintf(filename, "output_%d.bin", step);
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    fwrite(&Nr, sizeof(int), 1, fp);
    fwrite(&Nth, sizeof(int), 1, fp);
    fwrite(&Nphi, sizeof(int), 1, fp);
    fwrite(r_grid, sizeof(double), Nr, fp);
    fwrite(theta_grid, sizeof(double), Nth, fp);
    fwrite(phi_grid, sizeof(double), Nphi, fp);
    fwrite(density, sizeof(double), Nr*Nth*Nphi, fp);
    fwrite(pressure, sizeof(double), Nr*Nth*Nphi, fp);
    fwrite(velocity, sizeof(double), NDIM*Nr*Nth*Nphi, fp);

    fclose(fp);
}

int main()
{
    double *r_grid = (double *)malloc(Nr * sizeof(double));
    double *theta_grid = (double *)malloc(Nth * sizeof(double));
    double *phi_grid = (double *)malloc(Nphi * sizeof(double));

    for (int i = 0; i < Nr; i++) {
        r_grid[i] = r_min + i * (r_max - r_min) / (Nr - 1);
    }
    for (int j = 0; j < Nth; j++) {
        theta_grid[j] = theta_min + j * (theta_max - theta_min) / (Nth - 1);
    }
    for (int k = 0; k < Nphi; k++) {
        phi_grid[k] = phi_min + k * (phi_max - phi_min) / (Nphi - 1);
    }

    double *density = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *pressure = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *velocity = (double *)malloc(NDIM * Nr * Nth * Nphi * sizeof(double));

    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                density[idx] = exp(-r_grid[i] / 5.0) * pow(sin(theta_grid[j]), 2); 
                pressure[idx] = density[idx]; 
                velocity[TT * Nr * Nth * Nphi + idx] = 1.0; 
                velocity[RR * Nr * Nth * Nphi + idx] = sqrt(1.0 / r_grid[i]) * cos(theta_grid[j]);
                velocity[TH * Nr * Nth * Nphi + idx] = 0.0;
                velocity[PHI * Nr * Nth * Nphi + idx] = sqrt(1.0 / r_grid[i]) * sin(theta_grid[j]) + rand() / (double)RAND_MAX * 0.01;
            }
        }
    }

    double dt = 0.001;
    int num_steps = 50; 

    evolve_fluid(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi, dt, num_steps);
    write_binary_file(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi, num_steps);
    printf("Finished\n");

    free(r_grid);
    free(theta_grid);
    free(phi_grid);
    free(density);
    free(pressure);
    free(velocity);

    return 0;
}
