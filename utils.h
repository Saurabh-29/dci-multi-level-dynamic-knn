#ifndef UTILS_H
#define UTILS_H

#include <iostream> 
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <cstdio>

using namespace std; 

extern "C" {
#ifdef USE_MKL
#define DGEMM dgemm
#else
#define DGEMM dgemm_
#endif
extern void DGEMM(const char* const transa, const char* const transb, const int* const m, const int* const n, const int* const k, const double* const alpha, const double* const A, const int* const lda, const double* const B, const int* const ldb, const double* const beta, double* const C, const int* const ldc);

void matmul(const int M, const int N, const int K, const double* const A, const double* const B, double* const C);
}

void gen_data(double* const data, const int ambient_dim, const int intrinsic_dim, const int num_points);

double rand_normal();

double compute_dist(const double* const vec1, const double* const vec2, const int dim);

#endif // UTILS_H
