#include "geodesics.h"

void christoffel(long double g[4][4], long double christoffel[4][4][4])
{
    for (int mu = 0; mu < 4; mu++) {
        for (int beta = 0; beta < 4; beta++) {
            for (int nu = 0; nu < 4; nu++) {
                long double sum = 0;
                for (int sigma = 0; sigma < 4; sigma++) {
                    sum += 0.5 * (g[mu][sigma] * (g[sigma][beta] + g[beta][sigma] - g[beta][nu]));
                }
                christoffel[mu][beta][nu] = sum;
                printf("Christoffel[%d][%d][%d] = %Lf\n", mu, beta, nu, christoffel[mu][beta][nu]);
            }
        }
    }
}

void riemann(long double g[4][4], long double christoffel[4][4][4], long double riemann[4][4][4][4])
{
    for (int i = 1; i < 4; i++) {
        for (int j = 1; j < 4; j++) {
            for (int k = 1; k < 4; k++) {
                for (int l = 1; l < 4; l++) {
                    long double sum = 0;
                    sum += (1 / (2 * g[0][0])) * (christoffel[k][i][l] * christoffel[l][j][k] - christoffel[k][j][l] * christoffel[l][i][k]);
                    sum += (1 / (2 * g[0][0])) * (christoffel[k][i][l] * (g[l][j] * g[k][k] - g[l][k] * g[j][k]) - christoffel[k][j][l] * (g[l][i] * g[k][k] - g[l][k] * g[i][k]));
                    sum += (1 / (2 * g[0][0])) * (christoffel[l][i][k] * (g[l][j] * g[k][k] - g[l][k] * g[j][k]) - christoffel[l][j][k] * (g[l][i] * g[k][k] - g[l][k] * g[i][k]));
                    riemann[i][j][k][l] = sum;
					printf("Riemann[%d][%d][%d][%d] = %Lf\n", i, j, k, l, riemann[i][j][k][l]);
                }
            }
        }
    }
}

// Ajoutez ces lignes juste avant la fonction main()
long double (*geodesic_points)[4] = NULL;
int num_points = 0;

void store_geodesic_point(long double x[4], long double lambda)
{
    long double r = x[1];
    long double theta = x[2];
    long double phi = x[3];
    long double sin_theta = sin(theta);
    long double cos_theta = cos(theta);
    long double sin_phi = sin(phi);
    long double cos_phi = cos(phi);

    geodesic_points = realloc(geodesic_points, (num_points + 1000) * sizeof(*geodesic_points));
    if (geodesic_points == NULL)
    {
        fprintf(stderr, "Error: failed to allocate memory for geodesic_points\n");
        exit(1);
    }

    geodesic_points[num_points][0] = r * sin_theta * cos_phi;
    geodesic_points[num_points][1] = r * sin_theta * sin_phi;
    geodesic_points[num_points][2] = r * cos_theta;
	geodesic_points[num_points][3] = lambda; // Stocker la valeur de lambda comme scalaire

    num_points++;
}

void geodesic(long double x[4], long double v[4], long double lambda_max, long double christoffel[4][4][4], long double step_size, void (*store_point)(long double[], long double))
{
    long double k1[4], k2[4], k3[4], k4[4];
    long double temp_x[4], temp_v[4];
    ldouble_a32 lambda = 0;
	int step = 0;
    while (lambda < lambda_max) 
	{
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k1[mu] = v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k1[mu] -= christoffel[mu][alpha][beta] * v[alpha] * v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + 0.5 * step_size * k1[mu];
            temp_v[mu] = v[mu] + 0.5 * step_size * k1[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k2[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k2[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + 0.5 * step_size * k2[mu];
            temp_v[mu] = v[mu] + 0.5 * step_size * k2[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k3[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k3[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            temp_x[mu] = x[mu] + step_size * k3[mu];
            temp_v[mu] = v[mu] + step_size * k3[mu];
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            k4[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int beta = 0; beta < 4; beta++) {
                    k4[mu] -= christoffel[mu][alpha][beta] * temp_v[alpha] * temp_v[beta];
                }
            }
        }
		#pragma omp parallel for
        for (int mu = 0; mu < 4; mu++) {
            x[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
            v[mu] += step_size * (k1[mu] + 2 * k2[mu] + 2 * k3[mu] + k4[mu]) / 6.0;
        }
		store_geodesic_point(x, lambda);
		if (step % 1000 == 0)
			printf("Step %d: x = %Lf, y = %Lf, z = %Lf\n", step, x[0], x[1], x[2]);
		step++;
        lambda += step_size;
	}
}



void write_vtk_file(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return;
    }

    fprintf(file, "# vtk DataFile Version 3.0\n");
    fprintf(file, "Geodesic Points\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET POLYDATA\n");
    fprintf(file, "POINTS %d double\n", num_points);

    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%Lf %Lf %Lf\n", geodesic_points[i][0], geodesic_points[i][1], geodesic_points[i][2]);
    }

    fprintf(file, "LINES %d %d\n", num_points - 1, 3 * (num_points - 1));

    for (int i = 0; i < num_points - 1; ++i)
    {
        fprintf(file, "2 %d %d\n", i, i + 1);
    }

    fprintf(file, "POINT_DATA %d\n", num_points);
    fprintf(file, "SCALARS lambda double\n");
    fprintf(file, "LOOKUP_TABLE default\n");

    for (int i = 0; i < num_points; ++i)
    {
        fprintf(file, "%Lf\n", geodesic_points[i][3]);
    }
    fclose(file);
}

static inline double sqrt_asm(double n)
{
    double result;
    asm("fld %1;"
        "fsqrt;"
        "fstp %0;"
        : "=m" (result)
        : "m" (n));
    return result;
}


int main()
{
	ldouble_a32 Rs = 2 * G * M / powf(c, 2);

	
    long double x[4] = {Rs * 3, M_PI / 2, M_PI * 2, 1.0}; // Position initiale (r, θ, φ, t)
    long double v[4] = {100, 10.0, 0.0, 0.0};        // Vitesse initiale (dr/dλ, dθ/dλ, dφ/dλ, dt/dλ)

    ldouble_a32 H = 1.0;

    ldouble_a32 r = sqrt_asm(powf(x[1], 2) + powf(a, 2) * powf(cos(x[2]), 2));
    ldouble_a32 f = 2 * r * H / (powf(r, 2) + powf(a, 2) * powf(cos(x[2]), 2));
	ldouble_a32 rho2 = powf(r, 2) + powf(a, 2) * powf(cos(M_PI / 2), 2);
	ldouble_a32 delta = powf(r, 2) - Rs * r + powf(a, 2);
	ldouble_a32 Sigma = powf(r, 2) + powf(a, 2);
	ldouble_a32 A = (powf(r, 2) + powf(a, 2)) * Sigma + Rs * r * powf(a, 2) * powf(sin(M_PI / 2), 2);

	long double g_kerr[4][4] = {
		{-(1 - Rs * r / rho2), 0, 0, -Rs * r * a * powf(sinf(M_PI / 2), 2) / rho2},
		{0, rho2 / delta, 0, 0},
		{0, 0, rho2, 0},
		{-Rs * r * a * powf(sinf(M_PI / 2), 2) / rho2, 0, 0, (Sigma * powf(sinf(M_PI / 2), 2)) / rho2}
	};

    ldouble_a32 l0 = -1.0;
    ldouble_a32 l1 = (r * x[1] + a * x[3] * sinf(x[2])) / (powf(r, 2) + powf(a, 2));
    ldouble_a32 l2 = 0.0;
    ldouble_a32 l3 = (a * sin(x[2]) * (powf(r, 2) + powf(a, 2) - 2 * r)) / (powf(r, 2) + powf(a, 2));
    
	long double g_schwarzschild[4][4] = {0};
	ldouble_a32 denominator = 1 - 2 * G * M / (powf(c, 2));
	if (denominator != 0.0) {
		g_schwarzschild[0][0] = -denominator;
		g_schwarzschild[1][1] = 1 / denominator;
		g_schwarzschild[2][2] = powf(r, 2);
		g_schwarzschild[3][3] = powf(r, 2) * powf(sinf(M_PI / 2), 2);
	}

    long double g_kerr_schild[4][4] = {
        {-1 + f * l0 * l0, f * l0 * l1, f * l0 * l2, f * l0 * l3},
        {f * l1 * l0, 1 + f * l1 * l1, f * l1 * l2, f * l1 * l3},
        {f * l2 * l0, f * l2 * l1, 1 + f * l2 * l2, f * l2 * l3},
        {f * l3 * l0, f * l3 * l1, f * l3 * l2, 1 + f * l3 * l3}
    };

    long double christoffel_sym[4][4][4] = {0};
    long double riemann_tensor[4][4][4][4] = {0};

    christoffel(g_kerr, christoffel_sym);
    riemann(g_kerr, christoffel_sym, riemann_tensor);
    geodesic(x, v, 20.0, christoffel_sym, 0.0001, store_geodesic_point);
    write_vtk_file("geodesic.vtk");
    free(geodesic_points);

    return 0;
}