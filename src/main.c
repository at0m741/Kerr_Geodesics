#include "../headers/geodesics.h"

extern double (*geodesic_points)[4];
extern int num_points;

void print_arch()
{
	printf("Architecture: %s\n", ARCH);
	printf("OS: %s\n", Plateform);
	if (ALIGNMENT == 64)
	{
		printf("AVX-512 is supported\n");
		printf("AVX2 is supported\n");
		printf("SSE is supported\n");
		printf("SSE2 is supported\n");
	}
	else if (ALIGNMENT == 32)
	{
		printf("AVX2 is supported\n");
		printf("SSE2 is supported\n");
		printf("SSE is supported\n");
	}
	else
	{
		printf("SSE is supported\n");
	}
	printf("Alignment: %d bits\n", ALIGNMENT);
}

int main(int argc, char **argv)
{
	(void)argc;
	(void)argv;
	ldouble_a Q = 1.0;
	ldouble_a Rs = 2 * G * M / powf(c, 2);
    double v[4] = {-7.0, 1.01, 1.0, 27.0};        // Vitesse initiale (dr/dλ, dθ/dλ, dφ/dλ, dt/dλ)
    double x[4] = {Rs, M_PI / 2, 60 * M_PI, 0.0}; // Position initiale (r, θ, φ, t)
    ldouble_a r = sqrt_asm(powf(x[1], 2) + powf(a, 2) * powf(cos(x[2]), 2));
	ldouble_a rho2_kn = powf(r, 2) + powf(a, 2) * powf(cos(x[1]), 2);
	ldouble_a delta_kn = powf(r, 2) - Rs * r + powf(a, 2) + powf(Q, 2);
	ldouble_a Sigma_kn = powf((powf(r, 2) + powf(a, 2)), 2) - powf(a, 2) * delta_kn *\
					     powf(sin(x[1]), 2);
	double g_kerr_newman[4][4] = {0};
    double christoffel_sym[4][4][4] = {0};

	g_kerr_newman[0][0] = -(1 - (Rs * r - powf(Q, 2)) / rho2_kn);
	g_kerr_newman[0][3] = -(Rs * r - powf(Q, 2)) * a * powf(sin(x[1]), 2) / rho2_kn;
	g_kerr_newman[1][1] = rho2_kn / delta_kn;
	g_kerr_newman[2][2] = rho2_kn;
	g_kerr_newman[3][0] = g_kerr_newman[0][3];
	g_kerr_newman[3][3] = (powf(r, 2) + powf(a, 2) + (Rs * r - powf(Q, 2)) * powf(a, 2) * \
						   powf(sin(x[1]), 2) / rho2_kn) * powf(sin(x[1]), 2);

    double riemann_tensor[4][4][4][4] = {0};

    christoffel(g_kerr_newman, christoffel_sym);
    riemann(g_kerr_newman, christoffel_sym, riemann_tensor);
    geodesic(x, v, max_dt, christoffel_sym, DT, store_geodesic_point);
    write_vtk_file("geodesic.vtk");
	print_arch();
	if (geodesic_points != NULL)
	{
		free(geodesic_points);
	}
	return 0;
}
