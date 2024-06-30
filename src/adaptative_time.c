/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   adaptative_time.c                                  :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ltouzali <ltouzali@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2024/06/25 16:59:23 by ltouzali          #+#    #+#             */
/*   Updated: 2024/06/30 14:53:38 by ltouzali         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../headers/geodesics.h"

/* 
    TODO : Implement adaptative time step and change for a rk45 method for better precision
            then compare k4 and k5 to get the error and adjust the step size
*/
double calculate_max_speed(__m256d christoffel[4][4][4], __m256d x[4], __m256d v[4]) {return 1.0;}

void adjust_step_size(__m256d christoffel[4][4][4], __m256d x[4], __m256d v[4], __m256d step_size, double delta_x, double safety_factor) {
    double max_speed = calculate_max_speed(christoffel, x, v);
    double new_step_size = safety_factor * delta_x / max_speed;
    if (new_step_size < _mm256_cvtsd_f64(step_size)) {
        step_size = VEC_SET1_PD(new_step_size);
    }
}