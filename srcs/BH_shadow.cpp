#include <Geodesics.h>

extern double a;

double calculate_impact_parameter(double p_t, double p_phi, double g_tt, double g_tphi, double g_phiphi) {
    return - (g_tphi * p_t + g_phiphi * p_phi) / (g_tt * p_t + g_tphi * p_phi);
}

double calculate_emission_angle(double p_r, double p_phi, double g_rr, double g_phiphi) {
    return atan2(sqrt(g_phiphi) * p_phi, sqrt(g_rr) * p_r) * 180.0 / M_PI; 
}

double b_critique_kerr(double a, int sense) {
    double r_ph = 2.0 * M;
    for (int i = 0; i < 10; i++) {  
        double f = 2 * r_ph - 3 * M + 4 * a * sqrt(M / r_ph) * sense;
        double df = 2 - 2 * a * sqrt(M / (r_ph * r_ph * r_ph)) * sense;
        r_ph -= f / df;
    }
    return (r_ph * r_ph + a * a) / (a + sense * sqrt(r_ph * r_ph * r_ph / M));
}

void compute_photon_properties(double g[4][4], double p[4]) {
    double b = calculate_impact_parameter(p[0], p[3], g[0][0], g[0][3], g[3][3]);
    double alpha = calculate_emission_angle(p[1], p[3], g[1][1], g[3][3]);

    printf("Impact parameter b = %f\n", b);
    printf("Emission angle alpha = %f\n", alpha);
	if (b < b_critique_kerr(a, 1) && b > b_critique_kerr(a, -1)) {
		printf("Photon is captured by the black hole\n");
	} else {
		printf("Photon is not captured by the black hole\n");
	}
}
