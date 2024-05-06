#include "geodesics.h"

#ifdef __INTEL_COMPILER
	__assume_aligned(geodesic_points, ALIGNMENT);
	__assume_aligned(num_points, ALIGNMENT);
	void __init_KNL_MPI(int argc, char **argv)
	{
		#ifdef __INTEL_COMPILER
			__assume_aligned(geodesic_points, ALIGNMENT);
			__assume_aligned(num_points, ALIGNMENT);
		#endif
		#ifdef __ICC
			#pragma vector aligned
			#pragma simd aligned
			#pragma ivdep
			#pragma omp simd aligned(geodesic_points: ALIGNMENT)
		#endif
		MPI_init(&argc, &argv);
		printf("MPI has been initialized over %s\n", ARCH);
	}

#elif USE_MPI
	void __init_KNL_MPI(int argc, char **argv)
	{
		MPI_init(&argc, &argv);
		printf("MPI has been initialized over %s\n", ARCH);
	}
#endif