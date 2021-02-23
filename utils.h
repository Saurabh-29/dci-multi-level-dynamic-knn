#include <iostream> 
#include <set>
#include <iterator> 
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdio>

using namespace std; 


// struct Val
// {
//     double val;
//     // double *data_ptr;
//     int global_id;
//     int local_id;
// };

extern "C" {
#ifdef USE_MKL
#define DGEMM dgemm
#else
#define DGEMM dgemm_
#endif
extern void DGEMM(const char* const transa, const char* const transb, const int* const m, const int* const n, const int* const k, const double* const alpha, const double* const A, const int* const lda, const double* const B, const int* const ldb, const double* const beta, double* const C, const int* const ldc);

void matmul(const int M, const int N, const int K, const double* const A, const double* const B, double* const C) {
    const char TRANSA = 'T';
    const char TRANSB = 'N';
    const double ALPHA = 1.; 
    const double BETA = 0.; 
    const int LDA = K;
    const int LDB = K;
    const int LDC = M;
    DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
    // dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}
}
void gen_data(double* const data, const int ambient_dim, const int intrinsic_dim, const int num_points) {
    int i;
    
    double* latent_data;
    double* transformation;
    assert(posix_memalign((void **)&latent_data, 64, sizeof(double)*intrinsic_dim*num_points) == 0);
    assert(posix_memalign((void **)&transformation, 64, sizeof(double)*intrinsic_dim*ambient_dim) == 0);

    //double* latent_data = (double *)memalign(64, sizeof(double)*intrinsic_dim*num_points);
    //double* transformation = (double *)memalign(64, sizeof(double)*intrinsic_dim*ambient_dim);
    for (i = 0; i < intrinsic_dim*num_points; i++) {
        //latent_data[i] = 2 * ((double)rand() / RAND_MAX) - 1;
        latent_data[i] = 2 * drand48() - 1;
    }
    for (i = 0; i < intrinsic_dim*ambient_dim; i++) {
        //transformation[i] = 2 * ((double)rand() / RAND_MAX) - 1;
        transformation[i] = 2 * drand48() - 1;
    }
    // Assuming column-major layout, transformation is intrisic_dim x ambient_dim, 
    // latent_data is intrinsic_dim x num_points, data is ambient_dim x num_points
    matmul(ambient_dim, num_points, intrinsic_dim, transformation, latent_data, data);
    free(latent_data);
    free(transformation);
}

double rand_normal() {
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            //double U1 = (double)rand() / RAND_MAX;
            //double U2 = (double)rand() / RAND_MAX;
            double U1 = drand48();
            double U2 = drand48();
            
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

double compute_dist(const double* const vec1, const double* const vec2, const int dim) {
    int i;
    double sq_dist = 0.0;
    for (i = 0; i < dim; i++) {
        sq_dist += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    return sqrt(sq_dist);
}
