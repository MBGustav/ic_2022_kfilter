#ifndef _MKL_OMP_HPP_
#define _MKL_OMP_HPP_


#include "mkl_common.h"
#include "mkl.h"

static struct mkl_matrix global_FP, global_GQ, global_FPF;
static gemm_conf conf1, conf2, conf3, conf4;

void omp_AllocCovMatrix(struct mkl_matrix &F,struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q)
{
    global_FP.allocate( F.row, P.col);
    global_GQ.allocate( G.row, Q.col);
    global_FPF.allocate(F.row, F.col);

    //gemm1(F,P)
    conf1 = set_gemm_conf(F, P, false, false);
    //gemm2(G,Q)
    conf2 = set_gemm_conf( G, Q, false, false);
    //gemm3(FP, Ft)
    conf3 = set_gemm_conf(global_FP, F, false, true);    
    //gemm4(GQ, Gt)
    conf4 = set_gemm_conf(global_GQ, G, false, true);
}


void  omp_FreeCovMatrix()
{
    global_FP.deallocate();
    global_GQ.deallocate();
    global_FPF.deallocate();
}

void cblas_CovMatrix(struct mkl_matrix &F, struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q, struct mkl_matrix &propP)
{    
    __data_type alpha = 1.0, beta =  0.0f;

    #pragma omp parallel sections shared(F, P, global_FP, global_FPF, alpha, beta, conf1, conf3)
    {
        #pragma omp section
        #if DATA_PRECISION == 1
        cblas_dgemm(CblasRowMajor, conf1.trA, conf1.trB, 
                    conf1.m, conf1.n, conf1.k, alpha, F.data, conf1.lda, 
                    P.data, conf1.ldb, beta, global_FP.data, conf1.ldc);
        #else
        cblas_sgemm(CblasRowMajor, conf1.trA, conf1.trB, 
                    conf1.m, conf1.n, conf1.k, alpha, F.data, conf1.lda, 
                    P.data, conf1.ldb, beta, global_FP.data, conf1.ldc);
        #endif

        #pragma omp section
        #if DATA_PRECISION == 1
        cblas_dgemm(CblasRowMajor, conf3.trA, conf3.trB, 
                    conf3.m, conf3.n, conf3.k, alpha, global_FP.data, conf3.lda, 
                    F.data, conf3.ldb, beta, global_FPF.data, conf3.ldc);
        #else
        cblas_sgemm(CblasRowMajor, conf3.trA, conf3.trB, 
                    conf3.m, conf3.n, conf3.k, alpha, global_FP.data, conf3.lda, 
                    F.data, conf3.ldb, beta, global_FPF.data, conf3.ldc);
        #endif
    }

    #pragma omp parallel sections shared(G, Q, global_GQ, propP, alpha, beta, conf2, conf4)
    {
        #pragma omp section
        #if DATA_PRECISION == 1
        cblas_dgemm(CblasRowMajor, conf2.trA, conf2.trB, 
                    conf2.m, conf2.n, conf2.k, alpha, G.data, conf2.lda, 
                    Q.data, conf2.ldb, beta, global_GQ.data, conf2.ldc);
        #else
        cblas_sgemm(CblasRowMajor, conf2.trA, conf2.trB, 
                    conf2.m, conf2.n, conf2.k, alpha, G.data, conf2.lda, 
                    Q.data, conf2.ldb, beta, global_GQ.data, conf2.ldc);
        #endif

        #pragma omp section
        #if DATA_PRECISION == 1
        cblas_dgemm(CblasRowMajor, conf4.trA, conf4.trB, 
                    conf4.m, conf4.n, conf4.k, alpha, global_GQ.data, conf4.lda, 
                    G.data, conf4.ldb, beta, propP.data, conf4.ldc);
        #else
        cblas_sgemm(CblasRowMajor, conf4.trA, conf4.trB, 
                    conf4.m, conf4.n, conf4.k, alpha, global_GQ.data, conf4.lda, 
                    G.data, conf4.ldb, beta, propP.data, conf4.ldc);
        #endif
    }

    #pragma omp taskwait
    #if DATA_PRECISION == 1
    cblas_daxpy(propP.size, alpha, global_FPF.data, 1, propP.data, 1);
    #else
    cblas_saxpy(propP.size, alpha, global_FPF.data, 1, propP.data, 1);
    #endif
}

void omp_CovMatrix(mkl_matrix &F,
                   mkl_matrix &P,
                   mkl_matrix &G, 
                   mkl_matrix &Q,
                   mkl_matrix &propP)
{    
    __data_type alpha = 1.0, beta =  0.0f;

    #pragma omp parallel for
    for (int i = 0; i < conf1.m; ++i) {
        #pragma omp simd
        for (int j = 0; j < conf1.n; ++j) {
            global_FP.data[i * conf1.lda + j] = 0.0f;
            for (int k = 0; k < conf1.k; ++k) {
                global_FP.data[i * conf1.lda + j] += alpha * F.data[i * conf1.lda + k] * P.data[k * conf1.ldb + j];
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < conf3.m; ++i) {
        #pragma omp simd
        for (int j = 0; j < conf3.n; ++j) {
            global_FPF.data[i * conf3.lda + j] = 0.0f;
            for (int k = 0; k < conf3.k; ++k) {
                global_FPF.data[i * conf3.lda + j] += alpha * global_FP.data[i * conf3.lda + k] * F.data[k * conf3.ldb + j];
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < conf2.m; ++i) {
        #pragma omp simd
        for (int j = 0; j < conf2.n; ++j) {
            global_GQ.data[i * conf2.lda + j] = 0.0f;
            for (int k = 0; k < conf2.k; ++k) {
                global_GQ.data[i * conf2.lda + j] += alpha * G.data[i * conf2.lda + k] * Q.data[k * conf2.ldb + j];
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < conf4.m; ++i) {
        #pragma omp simd
        for (int j = 0; j < conf4.n; ++j) {
            propP.data[i * conf4.lda + j] = 0.0f;
            for (int k = 0; k < conf4.k; ++k) {
                propP.data[i * conf4.lda + j] += alpha * global_GQ.data[i * conf4.lda + k] * G.data[k * conf4.ldb + j];
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < propP.size; ++i) {
        propP.data[i] += alpha * global_FPF.data[i];
    }
}

// Solução, Implementar via kernel os codigos de openMP 5.0 em vez da chamada de função cblas..
void omp_CovMatrix_offload(mkl_matrix &F,
                           mkl_matrix &P,
                           mkl_matrix &G, 
                           mkl_matrix &Q,
                           mkl_matrix &propP)
{    
    __data_type alpha = 1.0, beta =  0.0f;

    #pragma omp parallel sections shared(F,         \
                                P, global_FP,       \
                                global_FPF, alpha,  \
                                beta, conf1, conf3)
    {
    #pragma omp section
    #pragma omp target teams distribute parallel for\
        map(to: F.data[0:F.size], P.data[0:P.size]) \
        map(alloc: global_FP.data[0:global_FP.size])
        for (int i = 0; i < conf1.m; ++i) 
            #pragma omp simd
            for (int j = 0; j < conf1.n; ++j) {
                global_FP.data[i * conf1.lda + j] = 0.0f;
                for (int k = 0; k < conf1.k; ++k) {
                    global_FP.data[i * conf1.lda + j] += alpha * 
                            F.data[i * conf1.lda + k] * 
                            P.data[k * conf1.ldb + j];
                }
            }
    
    #pragma omp section
    #pragma omp target teams distribute parallel for \
        map(to: global_FP.data[0:global_FP.size])    \
        map(tofrom: global_FPF.data[0:global_FPF.size])
        for (int i = 0; i < conf3.m; ++i) {

            #pragma omp simd
            for (int j = 0; j < conf3.n; ++j) {
                global_FPF.data[i * conf3.lda + j] = 0.0f;
                for (int k = 0; k < conf3.k; ++k){

                    global_FPF.data[i * conf3.lda + j] += alpha * 
                     global_FP.data[i * conf3.lda + k] * 
                             F.data[k * conf3.ldb + j];
                } 
            }
        }
    }

        #pragma omp target teams distribute parallel for map(to: G.data[0:G.size], Q.data[0:Q.size]) \
                                                         map(tofrom: global_GQ.data[0:global_GQ.size])
        for (int i = 0; i < conf2.m; ++i) {
            #pragma omp simd
            for (int j = 0; j < conf2.n; ++j) {
                global_GQ.data[i * conf2.lda + j] = 0.0f;
                for (int k = 0; k < conf2.k; ++k) {
                    global_GQ.data[i * conf2.lda + j] += alpha * G.data[i * conf2.lda + k] * Q.data[k * conf2.ldb + j];
                }
            }
        }

        #pragma omp target teams distribute parallel for map(to: global_GQ.data[0:global_GQ.size]) \
                                                         map(alloc: propP.data[0:propP.size])
        for (int i = 0; i < conf4.m; ++i) {
            #pragma omp simd
            for (int j = 0; j < conf4.n; ++j) {
                propP.data[i * conf4.lda + j] = 0.0f;
                for (int k = 0; k < conf4.k; ++k) {
                    propP.data[i * conf4.lda + j] += alpha * global_GQ.data[i * conf4.lda + k] * G.data[k * conf4.ldb + j];
                }
            }
        }
    

    for (int i = 0; i < propP.size; ++i) {
        propP.data[i] += alpha * global_FPF.data[i];
    }
}

#endif //_MKL_OMP_HPP_
