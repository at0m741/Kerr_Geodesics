#include "geodesics.h"

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
	printf("Block size: %d bits\n", BLOCK_SIZE);
}

int main(int argc, char **argv)
{
	ldouble_a Rs = 2 * G * M / powf(c, 2);

	
    // double x[4] = {Rs * 4, M_PI / 2, M_PI * 2, Rs}; // Position initiale (r, θ, φ, t)
    // double v[4] = {100.0, 1.01, 1.0, 1.0};        // Vitesse initiale (dr/dλ, dθ/dλ, dφ/dλ, dt/dλ)
	double x[4] = {Rs * 4, M_PI / 2, 0.0, 0.0}; // Position initiale (t, r, θ, φ)
	double v[4] = {1.0, 0.0, 0.0, 0.1};          // Vitesse initiale (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
    ldouble_a H = 1.0;

    ldouble_a r = sqrt_asm(powf(x[1], 2) + powf(a, 2) * powf(cos(x[2]), 2));
    ldouble_a f = 2 * r * H / (powf(r, 2) + powf(a, 2) * powf(cos(x[2]), 2));
	ldouble_a rho2 = powf(r, 2) + powf(a, 2) * powf(cos(M_PI / 2), 2);
	ldouble_a delta = powf(r, 2) - Rs * r + powf(a, 2);
	ldouble_a Sigma = powf(r, 2) + powf(a, 2);
	ldouble_a A = (powf(r, 2) + powf(a, 2)) * Sigma + Rs * r * powf(a, 2) * powf(sin(M_PI / 2), 2);

	double g_kerr[4][4] = {
		{-(1 - Rs * r / rho2), 0, 0, -Rs * r * a * powf(sinf(M_PI / 2), 2) / rho2},
		{0, rho2 / delta, 0, 0},
		{0, 0, rho2, 0},
		{-Rs * r * a * powf(sinf(M_PI / 2), 2) / rho2, 0, 0, (Sigma * powf(sinf(M_PI / 2), 2)) / rho2}
	};

    ldouble_a l0 = -1.0;
    ldouble_a l1 = (r * x[1] + a * x[3] * sinf(x[2])) / (powf(r, 2) + powf(a, 2));
    ldouble_a l2 = 0.0;
    ldouble_a l3 = (a * sin(x[2]) * (powf(r, 2) + powf(a, 2) - 2 * r)) / (powf(r, 2) + powf(a, 2));
    double g_kerr_schild[4][4] = {
        {-1 + f * l0 * l0, f * l0 * l1, f * l0 * l2, f * l0 * l3},
        {f * l1 * l0, 1 + f * l1 * l1, f * l1 * l2, f * l1 * l3},
        {f * l2 * l0, f * l2 * l1, 1 + f * l2 * l2, f * l2 * l3},
        {f * l3 * l0, f * l3 * l1, f * l3 * l2, 1 + f * l3 * l3}
    };

	ldouble_a Q = 1.03;
	double g_kerr_newman[4][4] = {0};

	ldouble_a rho2_kn = powf(r, 2) + powf(a, 2) * powf(cos(x[1]), 2);
	ldouble_a delta_kn = powf(r, 2) - Rs * r + powf(a, 2) + powf(Q, 2);
	ldouble_a Sigma_kn = powf((powf(r, 2) + powf(a, 2)), 2) - powf(a, 2) * delta_kn * powf(sin(x[1]), 2);

	g_kerr_newman[0][0] = -(1 - (Rs * r - powf(Q, 2)) / rho2_kn);
	g_kerr_newman[0][3] = -(Rs * r - powf(Q, 2)) * a * powf(sin(x[1]), 2) / rho2_kn;
	g_kerr_newman[1][1] = rho2_kn / delta_kn;
	g_kerr_newman[2][2] = rho2_kn;
	g_kerr_newman[3][0] = g_kerr_newman[0][3];
	g_kerr_newman[3][3] = (powf(r, 2) + powf(a, 2) + (Rs * r - powf(Q, 2)) * powf(a, 2) * powf(sin(x[1]), 2) / rho2_kn) * powf(sin(x[1]), 2);

    double christoffel_sym[4][4][4] = {0};
    double riemann_tensor[4][4][4][4] = {0};

    christoffel(g_kerr, christoffel_sym);
    riemann(g_kerr, christoffel_sym, riemann_tensor);
    geodesic(x, v, 10.4, christoffel_sym, 0.00001, store_geodesic_point);
	printf("Rs = %f\n", Rs * 1000);
	print_arch();
    write_vtk_file("geodesic.vtk");
	free(geodesic_points);

    return 0;
}