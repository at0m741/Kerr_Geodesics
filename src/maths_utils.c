/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   maths_utils.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/21 19:12:53 by ltouzali          #+#    #+#             */
/*   Updated: 2024/09/01 02:40:17 by at0m             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"


VEC_TYPE _mm256_sqrt_simd_asm(VEC_TYPE x) {
    __asm__ __volatile__("vsqrtpd %ymm0, %ymm0");
    return x;
}

#ifdef AVX512
    double _mm512_sqrt_simd_asm(double x) {__asm__ __volatile__("sqrtpd %zmm0, %zmm0");}
#endif

VEC_TYPE _mm256_exp_pd(VEC_TYPE x) {
    double vals[4];
    VEC_STORE_PD(vals, x);
    for (int i = 0; i < 4; i++)
        vals[i] = expf(vals[i]);
    return VEC_LOADU_PD(vals);
}

VEC_TYPE _mm256_sin_pd(VEC_TYPE x) {
    double vals[4];
    VEC_STORE_PD(vals, x);
    for (int i = 0; i < 4; i++)
        vals[i] = sinf(vals[i]);
    return VEC_LOADU_PD(vals);
}

VEC_TYPE _mm256_cos_pd(VEC_TYPE x) {
    double vals[4];
    VEC_STORE_PD(vals, x);
    for (int i = 0; i < 4; i++)
        vals[i] = cosf(vals[i]);
    return VEC_LOADU_PD(vals);
}
