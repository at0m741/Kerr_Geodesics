#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//include for avx maths
#define ALIGNMENT 32
#define DLOOP  for(int j=0;j<NDIM;j++) for(int k=0;k<NDIM;k++)
#define NDIM 4
#define TT 0
#define RR 1
#define TH 2
#define PH 3
#define SMALL 1.e-40
#define DT 0.0000005
#define max_dt 2.9
#define a 0.9375
 

void Boyer_lindquist_coord(double *X, double *r, double *th)
{
	double  R0 = 1.0, Rin = 2.0, Rout = 10.0, hslope = 0.3;
    *r = exp(X[1]) + R0;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
    if (fabs(*th) < SMALL)
    {
        if ((*th) >= 0)
            *th = SMALL;
        if ((*th) < 0)
            *th = -SMALL;
    }
    if (fabs(M_PI - (*th)) < SMALL)
    {
        if ((*th) >= M_PI)
            *th = M_PI + SMALL;
        if ((*th) < M_PI)
            *th = M_PI - SMALL;
    }
    return;
}
void gcov(double *X, double gcov[][NDIM])
{
    DLOOP gcov[j][k] = 0.;
    double dXpdX[NDIM];

    double r, th;
	double R0 = 1.0;
	double hslope = 0.3;
    Boyer_lindquist_coord(X, &r, &th);
    double sth, cth, s2, rho2;
    cth = cos(th);
    sth = sin(th);
    s2 = sth * sth;
    rho2 = r * r + a * a * cth * cth;

    gcov[TT][TT] = (-1. + 2. * r / rho2);
    gcov[TT][1] = (2. * r / rho2);
    gcov[TT][3] = (-2. * a * r * s2 / rho2);
    gcov[1][TT] = gcov[TT][1];
    gcov[1][1] = (1. + 2. * r / rho2);
    gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2));
    gcov[2][2] = rho2;
    gcov[3][TT] = gcov[TT][3];
    gcov[3][1] = gcov[1][3];
    gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));

}


void christoffel_AVX(__m256d g[4][4], __m256d christoffel[4][4][4]) 
{
    #ifdef USE_MPI
        printf("MPI is defined for Christoffel\n");
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #elif _OPENMP
        printf("OpenMP is supported and used\n");
    #else
        printf("OpenMP is not supported\n");
    #endif

    __m256d (*g_aligned)[4] = aligned_alloc(ALIGNMENT, sizeof(__m256d[4][4]));
    __m256d (*christoffel_aligned)[4][4] = aligned_alloc(ALIGNMENT, sizeof(__m256d[4][4][4]));
    memcpy(g_aligned, g, sizeof(__m256d[4][4]));
    memcpy(christoffel_aligned, christoffel, sizeof(__m256d[4][4][4]));

    #pragma omp parallel for collapse(3)
    for (int mu = 0; mu < 4; mu++) {
        for (int beta = 0; beta < 4; beta++) {
            for (int nu = 0; nu < 4; nu++) {
                __m256d sum = _mm256_setzero_pd();
                for (int sigma = 0; sigma < 4; sigma++) {
                    __m256d g_mu_sigma = _mm256_broadcast_sd((const double*)&g_aligned[mu][sigma]);
                    __m256d g_sigma_beta = _mm256_broadcast_sd((const double*)&g_aligned[sigma][beta]);
                    __m256d g_beta_sigma = _mm256_broadcast_sd((const double*)&g_aligned[beta][sigma]);
                    __m256d g_beta_nu = _mm256_broadcast_sd((const double*)&g_aligned[beta][nu]);

                    __m256d term1 = _mm256_mul_pd(g_mu_sigma, g_sigma_beta);
                    __m256d term2 = _mm256_mul_pd(g_mu_sigma, g_beta_sigma);
                    __m256d term3 = _mm256_mul_pd(g_mu_sigma, g_beta_nu);

                    sum = _mm256_add_pd(sum, _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_sub_pd(_mm256_add_pd(term1, term2), term3)));
                }
                christoffel_aligned[mu][beta][nu] = sum;
            }
        }
    }


    memcpy(christoffel, christoffel_aligned, sizeof(__m256d[4][4][4]));
    free(g_aligned);
    free(christoffel_aligned);
}

void store_geodesic_point(__m256d x[4], __m256d lambda);

void geodesic(__m256d x[4], __m256d v[4], float lambda_max, __m256d christoffel[4][4][4], __m256d step_size)
{
    __m256d k1[4], k2[4], k3[4], k4[4], temp_x[4], temp_v[4];
    float lambda = 0.0;
    __m256i step = _mm256_set1_epi64x(0); 

    while (lambda < lambda_max)
    {
        #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++) {
            k1[mu] = v[mu];
            for (int nu = 0; nu < 4; nu++) {
                for (int beta = 0; beta < 4; beta++) {
                    __m256d product = _mm256_mul_pd(christoffel[mu][nu][beta], 
                                        _mm256_mul_pd(v[nu], v[beta]));
                    k1[mu] = _mm256_sub_pd(k1[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k1[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k1[mu]));
        }
        #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++){
            k2[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++){
                for (int beta = 0; beta < 4; beta++){
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], \
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k2[mu] = _mm256_sub_pd(k2[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k2[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k2[mu]));
        }
        #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++){
            k3[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++){
                for (int beta = 0; beta < 4; beta++){
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], \
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k3[mu] = _mm256_sub_pd(k3[mu], product);
                }
            }
            temp_x[mu] = _mm256_add_pd(x[mu], _mm256_mul_pd(step_size, k3[mu]));
            temp_v[mu] = _mm256_add_pd(v[mu], _mm256_mul_pd(step_size, k3[mu]));
        }
        #pragma omp simd aligned (temp_v, christoffel: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++){
            k4[mu] = temp_v[mu];
            for (int alpha = 0; alpha < 4; alpha++){
                for (int beta = 0; beta < 4; beta++){
                    __m256d product = _mm256_mul_pd(christoffel[mu][alpha][beta], \
                                        _mm256_mul_pd(temp_v[alpha], temp_v[beta]));
                    k4[mu] = _mm256_sub_pd(k4[mu], product);
                }
            }
        }
        #pragma omp simd aligned (x, v, k1, k2, k3, k4: ALIGNMENT)
        for (int mu = 0; mu < 4; mu++){
            x[mu] = _mm256_add_pd(x[mu], _mm256_div_pd(_mm256_mul_pd(step_size, _mm256_add_pd(k1[mu], \
                                            _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), k2[mu]), 
                                            _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), k3[mu]), k4[mu])))), 
                                            _mm256_set1_pd(6.0)));

            v[mu] = _mm256_add_pd(v[mu], _mm256_div_pd(_mm256_mul_pd(step_size, _mm256_add_pd(k1[mu], \
                                            _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), k2[mu]), 
                                            _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(2.0), k3[mu]), k4[mu])))), 
                                            _mm256_set1_pd(6.0)));
        }
        #pragma omp single
        {
            __m256d lambda_cast = _mm256_set1_pd(lambda);
            store_geodesic_point(x, lambda_cast);
            step = _mm256_add_epi64(step, _mm256_set1_epi64x(1));
            lambda += step_size[0];
        }
    }
}


__m256d sin_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = sin(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

__m256d cos_pd(__m256d x) {
    double vals[4];
    _mm256_storeu_pd(vals, x);
    for (int i = 0; i < 4; i++) {
        vals[i] = cos(vals[i]);
    }
    return _mm256_loadu_pd(vals);
}

void store_geodesic_point(__m256d x[4], __m256d lambda)
{
    __m256d (*x_points)[4];
    __m256i capacity = _mm256_set1_epi64x(0);
    __m256i num_points = _mm256_set1_epi64x(0);
    if (capacity[0] == 0) {
        capacity = _mm256_set1_epi64x(1000);
        __m256d (*x_points)[4] = aligned_alloc(ALIGNMENT, sizeof(__m256d[1000][4]));
        __m256d (*lambda_points) = aligned_alloc(ALIGNMENT, sizeof(__m256d[1000]));
        if (x_points == NULL || lambda_points == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for geodesic points\n");
            exit(1);
        }

        if (num_points[0] < capacity[0]) {
            memcpy(x_points, x_points, sizeof(__m256d[1000][4]));
            memcpy(lambda_points, lambda_points, sizeof(__m256d[1000]));
        }
    }

    __m256d r = x[1];
    __m256d th = x[2];
    __m256d ph = x[3];
    __m256d sin_th = sin_pd(th);
    __m256d cos_th = cos_pd(th);
    __m256d sin_ph = sin_pd(ph);
    __m256d cos_ph = cos_pd(ph);
    __m256d dt = x[4];

    x_points[num_points[0]][0] = r * sin_th * cos_ph;
    x_points[num_points[0]][1] = r * sin_th * sin_ph;
    x_points[num_points[0]][2] = r * cos_th;
    x_points[num_points[0]][3] = lambda;

    //print
    double x_val[4];
    for (int i = 0; i < 4; i++) {
        _mm256_storeu_pd(x_val, x_points[num_points[0]][i]);
        fprintf(stdout, "x[%d]: %f, y[%d]: %f, z[%d]: %f, t[%d]: %f\n",
                i, x_val[0], i, x_val[1], i, x_val[2], i, x_val[3]);
    }
}

int main() {
    double x_vals[4] = {50.0, M_PI / 2, M_PI / 2, 20.0};
    double v_vals[4] = {-20.2, 10.0, 12.0, 27.0};
    double g_vals[NDIM][NDIM] = {0};

    __m256d x[4], v[4], g[4][4], christoffel_avx[4][4][4];
    for (int i = 0; i < NDIM; i++) {
        x[i] = _mm256_set1_pd(x_vals[i]);
        v[i] = _mm256_set1_pd(v_vals[i]);
        for (int j = 0; j < NDIM; j++) {
            g[i][j] = _mm256_set1_pd(g_vals[i][j]);
        }
    }
    gcov(x_vals, g_vals); 
    for (int i = 0; i < NDIM; i++) {
        for (int j = 0; j < NDIM; j++) {
            g[i][j] = _mm256_set1_pd(g_vals[i][j]);
        }
    }
    christoffel_AVX(g, christoffel_avx);
    geodesic(x, v, max_dt, christoffel_avx, _mm256_set1_pd(DT));

    return 0;
}