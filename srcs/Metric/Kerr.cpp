
#include <Geodesics.h>

extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

void Metric::calculate_metric(const std::array<double, NDIM>& x, 
                              std::array<std::array<double, NDIM>, NDIM>& g,
                              std::array<std::array<double, NDIM>, NDIM>& g_inv) {
    Matrix matrix_obj;
    Metric metric_obj;
    double r = x[1];
    double theta = x[2];
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_theta2 = sin_theta * sin_theta;
    double cos_theta2 = cos_theta * cos_theta;
    
    double Sigma = r * r + a * a * cos_theta2;
    double Delta = r * r - 2.0 * M * r + a * a;
    
    g[0][0] = -(1.0 - (2.0 * M * r) / Sigma);
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = sin_theta2 * (r * r + a * a + (2.0 * M * r * a * a * sin_theta2) / Sigma);
    g[0][3] = - (2.0 * M * r * a * sin_theta2) / Sigma;
    g[3][0] = g[0][3];

    matrix_obj.inverse_matrix(g, g_inv);
}


