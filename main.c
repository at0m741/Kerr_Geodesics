#include "geodesics.h"

extern long double (*geodesic_points)[4];
extern int num_points;


int main(int argc, char **argv)
{
	ldouble_a32 Rs = 2 * G * M / powf(c, 2);

	
    // long double x[4] = {Rs * 4, M_PI / 2, M_PI * 2, Rs}; // Position initiale (r, θ, φ, t)
    // long double v[4] = {100.0, 1.01, 1.0, 1.0};        // Vitesse initiale (dr/dλ, dθ/dλ, dφ/dλ, dt/dλ)
	long double x[4] = {Rs * 4, M_PI / 2, 0.0, 0.0}; // Position initiale (t, r, θ, φ)
	long double v[4] = {1.0, 0.0, 0.0, 0.1};          // Vitesse initiale (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
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
    long double g_kerr_schild[4][4] = {
        {-1 + f * l0 * l0, f * l0 * l1, f * l0 * l2, f * l0 * l3},
        {f * l1 * l0, 1 + f * l1 * l1, f * l1 * l2, f * l1 * l3},
        {f * l2 * l0, f * l2 * l1, 1 + f * l2 * l2, f * l2 * l3},
        {f * l3 * l0, f * l3 * l1, f * l3 * l2, 1 + f * l3 * l3}
    };

	ldouble_a32 Q = 1.03;
	long double g_kerr_newman[4][4] = {0};

	ldouble_a32 rho2_kn = powf(r, 2) + powf(a, 2) * powf(cos(x[1]), 2);
	ldouble_a32 delta_kn = powf(r, 2) - Rs * r + powf(a, 2) + powf(Q, 2);
	ldouble_a32 Sigma_kn = powf((powf(r, 2) + powf(a, 2)), 2) - powf(a, 2) * delta_kn * powf(sin(x[1]), 2);

	g_kerr_newman[0][0] = -(1 - (Rs * r - powf(Q, 2)) / rho2_kn);
	g_kerr_newman[0][3] = -(Rs * r - powf(Q, 2)) * a * powf(sin(x[1]), 2) / rho2_kn;
	g_kerr_newman[1][1] = rho2_kn / delta_kn;
	g_kerr_newman[2][2] = rho2_kn;
	g_kerr_newman[3][0] = g_kerr_newman[0][3];
	g_kerr_newman[3][3] = (powf(r, 2) + powf(a, 2) + (Rs * r - powf(Q, 2)) * powf(a, 2) * powf(sin(x[1]), 2) / rho2_kn) * powf(sin(x[1]), 2);

    long double christoffel_sym[4][4][4] = {0};
    long double riemann_tensor[4][4][4][4] = {0};

    christoffel(g_kerr, christoffel_sym);
    riemann(g_kerr, christoffel_sym, riemann_tensor);
    geodesic(x, v, 20.4, christoffel_sym, 0.00001, store_geodesic_point);
	printf("Rs = %Lf\n", Rs * 1000);
    write_vtk_file("geodesic.vtk");
	free(geodesic_points);

    return 0;
}