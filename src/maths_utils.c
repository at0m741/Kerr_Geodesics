/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   maths_utils.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/21 19:12:53 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/21 19:12:53 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

__m256d _mm256_exp_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = exp(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

__m256d _mm256_sin_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = sin(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

__m256d _mm256_cos_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = cos(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}