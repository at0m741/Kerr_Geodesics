#include <Geodesics.h>

double (*geodesic_points)[5] = NULL;
int num_points = 0;
double a = 0.935;

int main(int argc, char **argv)
{
	double r0 = 20.0;
	double X[NDIM] = {0.0, r0, M_PI/4.0, 0.0};

	if (argc < 3) {
		printf("Usage: <options>\n");
		printf("Options:\n");
		printf("       -G <Spin value a> - Geodesic calculation\n");
		printf("       -L <Spin value a> - Light geodesics calculation\n");
		printf("       -R <Spin value a> - Riemann tensor calculation\n");
		printf("       -M <Spin value a> - Metric tensor calculation (a = 0 -> Schwarzschild or a > 0 -> Kerr)\n");
		
		return 0;
	}

	a = atof(argv[2]);
	if (strcmp(argv[1], "-R") == 0) {
		if (argc < 3) {
			printf("Usage: -R <Spin value a>\n");
			return 0;
		} 
		Riemann_tensor(); 
	} else if (strcmp(argv[1], "-G") == 0) {
		if (argc < 3) {
			printf("Usage: -G <Spin value a>\n");
			return 0;
		} 
		Geodesics_prob();	
	}else if (strncmp(argv[1], "-L", 2) == 0) {
		if (argc < 3) {
			printf("Usage: -L <Spin value a>\n");
			return 0;
		}
		light_geodesics_prob();
	} else if (strncmp(argv[1], "-M", 2) == 0) {
		if (argc < 3) {
			printf("Usage: -M <Spin value a>\n");
			return 0;
		}
		Metric_prob();
	} else {
		printf("Invalid option\n");
		return 0;
	}
}
