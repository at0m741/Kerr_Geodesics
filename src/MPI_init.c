/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   MPI_init.c                                         :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/17 15:51:37 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/17 15:51:37 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

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
		printf("MPI has been initialized over %s\n", ARCH);
	}

#elif USE_MPI
	void __init_KNL_MPI(int argc, char **argv)
	{
		(void)argc;
		(void)argv;
	}

	void mpi_preliminary_task() {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0) {
			printf("MPI first %d\n", rank);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
#endif	