/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   maths_utils.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/21 19:12:53 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/26 18:16:51 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

double _mm256_sqrt_simd_asm(double x) {__asm__ __volatile__("sqrtpd %ymm0, %ymm0");}
#ifdef __AVX512__
    double _mm512_sqrt_simd_asm(double x) {__asm__ __volatile__("sqrtpd %zmm0, %zmm0");}
#endif
__m256d _mm256_exp_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = expf(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

__m256d _mm256_sin_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = sinf(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

__m256d _mm256_cos_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = cosf(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}
