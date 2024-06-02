#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include <hdf5.h>

#define NDIM 4
#define TT 0
#define RR 1
#define TH 2
#define PHI 3
#define ALIGNMENT 64

int Nr = 100;
int Nth = 50;
int Nphi = 50;

double r_min = 1.0;
double r_max = 6.0;
double theta_min = 0.0;
double theta_max = M_PI;
double phi_min = 0.0;
double phi_max = 4.0 * M_PI ;

#define SMALL 1e-10
double a = 0.45;
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

double equation_of_state(double density) {
    double gamma = 1.5; 
    return density * pow(density, gamma - 1);
}

void compute_viscous_heating(double *density, double *pressure, double *velocity, double *viscous_heating, int Nr, int Nth, int Nphi) {
    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double rho = density[idx];
                double p = pressure[idx];
                
                double shear_viscosity = viscosity * rho; // Assuming a simple form of viscosity
                double dv_r = (velocity[RR * Nr * Nth * Nphi + idx + 1] - velocity[RR * Nr * Nth * Nphi + idx - 1]) / 2.0; // Simple finite difference
                
                viscous_heating[idx] = shear_viscosity * dv_r * dv_r;
            }
        }
    }
}

void update_energy(double *density, double *pressure, double *velocity, double *flux_energy, double *viscous_heating, double dt, int Nr, int Nth, int Nphi) {
    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            #pragma omp simd
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double energy_flux = flux_energy[idx] + pressure[idx] * velocity[RR * Nr * Nth * Nphi + idx];
                double heating = viscous_heating[idx];
                
                pressure[idx] += dt * (energy_flux + heating);
            }
        }
    }
}

void initialize_ADAF_torus(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid, double *phi_grid, int Nr, int Nth, int Nphi) {
    double r0 = 5.0; // Central radius of the torus
    double H = 2.0; // Scale height of the torus
    double K = 1.0; // Normalization constant for pressure
    double gamma = 1.5; // Adiabatic index for ADAF

    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;

                double r = r_grid[i];
                double theta = theta_grid[j];
                double phi = phi_grid[k];

                double rho = exp(-(r - r0) * (r - r0) / (2.0 * H * H)) * pow(sin(theta), 2);
                double p = K * pow(rho, gamma);

                density[idx] = rho;
                pressure[idx] = p;

                velocity[TT * Nr * Nth * Nphi + idx] = 1.0; // 4-velocity time component
                velocity[RR * Nr * Nth * Nphi + idx] = 0.0; // Radial velocity
                velocity[TH * Nr * Nth * Nphi + idx] = 0.0; // Theta velocity
                velocity[PHI * Nr * Nth * Nphi + idx] = sqrt(r) * sin(theta); // Keplerian azimuthal velocity
            }
        }
    }
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

void apply_local_viscosity(double *velocity, int Nr, int Nth, int Nphi, double *density)
{
    #pragma omp parallel for
    for (int i = 1; i < Nr - 1; i++) {
        for (int j = 1; j < Nth - 1; j++) {
            for (int k = 1; k < Nphi - 1; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double local_viscosity = viscosity * exp(-density[idx]);

                #pragma omp simd aligned(velocity: ALIGNMENT)
                for (int d = 0; d < NDIM; d++) {
                    int idx_rp = (i + 1) * Nth * Nphi + j * Nphi + k;
                    int idx_rm = (i - 1) * Nth * Nphi + j * Nphi + k;
                    int idx_tp = i * Nth * Nphi + (j + 1) * Nphi + k;
                    int idx_tm = i * Nth * Nphi + (j - 1) * Nth * Nphi + k;
                    int idx_pp = i * Nth * Nphi + j * Nphi + (k + 1);
                    int idx_pm = i * Nth * Nphi + j * Nphi + (k - 1);

                    double laplacian = (velocity[d * Nr * Nth * Nphi + idx_rp] +
                                        velocity[d * Nr * Nth * Nphi + idx_rm] +
                                        velocity[d * Nr * Nth * Nphi + idx_tp] +
                                        velocity[d * Nr * Nth * Nphi + idx_tm] +
                                        velocity[d * Nr * Nth * Nphi + idx_pp] +
                                        velocity[d * Nr * Nth * Nphi + idx_pm] -
                                        6.0 * velocity[d * Nr * Nth * Nphi + idx]);

                    velocity[d * Nr * Nth * Nphi + idx] += local_viscosity * laplacian;
                }
            }
        }
    }
}

void apply_boundary_conditions(double *density, double *pressure, double *velocity, int Nr, int Nth, int Nphi) {
    for (int j = 0; j < Nth; j++) {
        for (int k = 0; k < Nphi; k++) {
            int idx_inner = 0 * Nth * Nphi + j * Nphi + k;
            int idx_next = 1 * Nth * Nphi + j * Nphi + k;
            density[idx_inner] = density[idx_next];
            pressure[idx_inner] = pressure[idx_next];
            for (int d = 0; d < NDIM; d++) {
                velocity[d * Nr * Nth * Nphi + idx_inner] = -velocity[d * Nr * Nth * Nphi + idx_next];
            }
        }
    }

    for (int j = 0; j < Nth; j++) {
        for (int k = 0; k < Nphi; k++) {
            int idx_outer = (Nr-1) * Nth * Nphi + j * Nphi + k;
            int idx_prev = (Nr-2) * Nth * Nphi + j * Nphi + k;
            density[idx_outer] = density[idx_prev];
            pressure[idx_outer] = pressure[idx_prev];
            for (int d = 0; d < NDIM; d++) {
                velocity[d * Nr * Nth * Nphi + idx_outer] = -velocity[d * Nr * Nth * Nphi + idx_prev];
            }
        }
    }
}

// double equation_of_state(double density) {
//     double gamma = 5.0 / 3.0;
//     return density * pow(density, gamma - 1);
// }

void apply_forcing(double *density, double *pressure, double *velocity, int Nr, int Nth, int Nphi, double *r_grid) {
    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double forcing_term = 0.01 * sin(r_grid[i]);
                density[idx] += forcing_term;
                pressure[idx] += forcing_term;
                for (int d = 0; d < NDIM; d++) {
                    velocity[d * Nr * Nth * Nphi + idx] += forcing_term;
                }
            }
        }
    }
}

double compute_dt(double *velocity, int Nr, int Nth, int Nphi) {
    double max_velocity = 0.0;
    #pragma omp parallel for reduction(max:max_velocity)
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double v_local = sqrt(velocity[RR * Nr * Nth * Nphi + idx] * velocity[RR * Nr * Nth * Nphi + idx] +
                                      velocity[TH * Nr * Nth * Nphi + idx] * velocity[TH * Nr * Nth * Nphi + idx] +
                                      velocity[PHI * Nr * Nth * Nphi + idx] * velocity[PHI * Nr * Nth * Nphi + idx]);
                if (v_local > max_velocity) {
                    max_velocity = v_local;
                }
            }
        }
    }
    return 0.5 * fmin((r_max - r_min) / Nr, fmin((theta_max - theta_min) / Nth, (phi_max - phi_min) / Nphi)) / max_velocity;
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
void write_hdf5_file(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid,\
                     double *phi_grid, int Nr, int Nth, int Nphi, int step) {
    char filename[256];
    sprintf(filename, "output_%d.h5", step);

    hid_t file_id, dataset_id, dataspace_id;
    hsize_t dims[1];

    // Create a new HDF5 file
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Write the density array
    dims[0] = Nr * Nth * Nphi;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "density", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, density);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Write the pressure array
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "pressure", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Write the velocity array
    dims[0] = NDIM * Nr * Nth * Nphi;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "velocity", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocity);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Write the r_grid array
    dims[0] = Nr;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "r_grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, r_grid);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Write the theta_grid array
    dims[0] = Nth;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "theta_grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theta_grid);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Write the phi_grid array
    dims[0] = Nphi;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    dataset_id = H5Dcreate(file_id, "phi_grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi_grid);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Close the file
    H5Fclose(file_id);
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

void evolve_fluid(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid,\
                  double *phi_grid, int Nr, int Nth, int Nphi, int num_steps) {
    double *flux_density = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *flux_pressure = (double *)malloc(Nr * Nth * Nphi * sizeof(double));
    double *flux_velocity = (double *)malloc(NDIM * Nr * Nth * Nphi * sizeof(double));
    double *viscous_heating = (double *)malloc(Nr * Nth * Nphi * sizeof(double));

    for (int n = 0; n < num_steps; n++) {
        printf("Step %d\n", n);

        double dt = compute_dt(velocity, Nr, Nth, Nphi);
        compute_fluxes(density, pressure, velocity, flux_density, flux_pressure, flux_velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi);
        compute_viscous_heating(density, pressure, velocity, viscous_heating, Nr, Nth, Nphi);
        apply_local_viscosity(velocity, Nr, Nth, Nphi, density);
        apply_forcing(density, pressure, velocity, Nr, Nth, Nphi, r_grid);
        apply_boundary_conditions(density, pressure, velocity, Nr, Nth, Nphi);
        printf("Before update - Density at [50, 25, 25]: %f\n", density[50 * Nth * Nphi + 25 * Nphi + 25]);
        update_variables(density, pressure, velocity, flux_density, flux_pressure, flux_velocity, dt, Nr, Nth, Nphi);
        update_energy(density, pressure, velocity, flux_pressure, viscous_heating, dt, Nr, Nth, Nphi);
        printf("After update - Density at [50, 25, 25]: %f\n", density[50 * Nth * Nphi + 25 * Nphi + 25]);
        write_hdf5_file(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi, n);
    }

    free(flux_density);
    free(flux_pressure);
    free(flux_velocity);
    free(viscous_heating);
}



// Fonction pour initialiser le tore
void initialize_torus(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid, double *phi_grid, int Nr, int Nth, int Nphi) {
    double r0 = 10.0; // Rayon du tore
    double H = 5.0; // Échelle de hauteur du tore
    double K = 2.0; // Constante de normalisation pour la pression
    double gamma = 5.0 / 3.0; // Indice adiabatique

    #pragma omp parallel for
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;

                double r = r_grid[i];
                double theta = theta_grid[j];
                double phi = phi_grid[k];

                double rho = exp(-(r - r0) * (r - r0) / (2.0 * H * H)) * pow(sin(theta), 2);
                double p = K * pow(rho, gamma);

                density[idx] = rho;
                pressure[idx] = p;

                velocity[TT * Nr * Nth * Nphi + idx] = 1.0; // Vitesse temporelle (4-velocity component)
                velocity[RR * Nr * Nth * Nphi + idx] = sqrt(r); // Vitesse radiale
                velocity[TH * Nr * Nth * Nphi + idx] = 0.0; // Vitesse polaire
                velocity[PHI * Nr * Nth * Nphi + idx] = sqrt(r) * sin(theta); // Vitesse Keplerienne en azimutal
            }
        }
    }
}

//vtk file

void write_vtk_file(double *density, double *pressure, double *velocity, double *r_grid, double *theta_grid,\
                       double *phi_grid, int Nr, int Nth, int Nphi, int step)
{
    char filename[256];
    sprintf(filename, "output_%d.vtk", step);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening file");
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "vtk output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d %d\n", Nr, Nth, Nphi);
    fprintf(fp, "ORIGIN %f %f %f\n", r_min, theta_min, phi_min);
    fprintf(fp, "SPACING %f %f %f\n", (r_max - r_min) / (Nr - 1), (theta_max - theta_min) / (Nth - 1), (phi_max - phi_min) / (Nphi - 1));
    fprintf(fp, "POINT_DATA %d\n", Nr * Nth * Nphi);
    fprintf(fp, "SCALARS density float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                fprintf(fp, "%f\n", density[idx]);
            }
        }
    }
    fprintf(fp, "SCALARS pressure float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                fprintf(fp, "%f\n", pressure[idx]);
            }
        }
    }
    fprintf(fp, "VECTORS velocity float\n");
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                fprintf(fp, "%f %f %f\n", velocity[RR * Nr * Nth * Nphi + idx], velocity[TH * Nr * Nth * Nphi + idx], velocity[PHI * Nr * Nth * Nphi + idx]);
            }
        }
    }
    fprintf(fp, "SCALARS speed float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nth; j++) {
            for (int k = 0; k < Nphi; k++) {
                int idx = i * Nth * Nphi + j * Nphi + k;
                double v = sqrt(velocity[RR * Nr * Nth * Nphi + idx] * velocity[RR * Nr * Nth * Nphi + idx] +
                                velocity[TH * Nr * Nth * Nphi + idx] * velocity[TH * Nr * Nth * Nphi + idx] +
                                velocity[PHI * Nr * Nth * Nphi + idx] * velocity[PHI * Nr * Nth * Nphi + idx]);
                fprintf(fp, "%f\n", v);
            }
        }
    }
    printf("VTK file written\n");
    fclose(fp);
}

int main() {
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

    // Initialize the torus with ADAF conditions
    initialize_ADAF_torus(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi);

    int num_steps = 300;

    evolve_fluid(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi, num_steps);
    write_vtk_file(density, pressure, velocity, r_grid, theta_grid, phi_grid, Nr, Nth, Nphi, num_steps);
    printf("Finished\n");

    free(r_grid);
    free(theta_grid);
    free(phi_grid);
    free(density);
    free(pressure);
    free(velocity);

    return 0;
}
