
#include <Geodesics.h>
extern double (*geodesic_points)[5];
extern int num_points;
extern double a;

int Metric_prob() {
	Metric metric_obj;
	double r0 = 20.0;
    std::array<double, NDIM> X = {0.0, r0, M_PI/4.0, 0.0};;
	printf("=====================================================\n");
	metric_obj.calculate_metric(X, metric_obj.gcov, metric_obj.gcon);

	printf("=====================================================\n");
	std::array<std::array<double, NDIM>, NDIM> gcov_KN, gcon_KN;	
	metric_obj.calculate_metric_kerr_newman(X, gcov_KN, gcon_KN);

	printf("=====================================================\n");
	std::array<std::array<double, NDIM>, NDIM> gcov_kds, gcon_kds;
	metric_obj.calculate_metric_kds(X, gcov_kds, gcon_kds);
	return 0;
}

