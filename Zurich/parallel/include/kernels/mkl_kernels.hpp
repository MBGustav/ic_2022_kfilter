#ifndef _MKL_KERNELS_HPP_
#define _MKL_KERNELS_HPP_
//builtin libs
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "mkl_common.h"

using namespace sycl;
//global pointers for cov_matrix_prediction
static struct mkl_matrix global_FP, global_GQ, global_FPF;
static gemm_conf conf1, conf2, conf3, conf4;

#ifdef MKL_ONEAPI
//Global arrays for reduce allocation latency -> CovMatrix Prediction
void oneapi_AllocCovMatrix(struct mkl_matrix &F, struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q)
{
    global_FP.q = F.q;
    global_GQ.q = F.q;
    global_FPF.q = F.q;
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

void oneapi_FreeCovMatrix()
{
    global_FP.deallocate();
    global_GQ.deallocate();
    global_FPF.deallocate();
}


void oneapi_CovMatrix(sycl::queue &queue, struct mkl_matrix &F, struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q, struct mkl_matrix propP)
{
    
    sycl::event gemm1_done, gemm2_done, gemm3_done, gemm4_done, axpy_done;
    __data_type alpha = 1.0, beta =  0.0f;


    //send to execution
    gemm1_done = oneapi::mkl::blas::row_major::gemm(queue, conf1.trA, conf1.trB, 
                                                    conf1.m, conf1.n, conf1.k, alpha, F.data, conf1.lda,
                                                    P.data, conf1.ldb, beta, global_FP.data, conf1.ldc);
    
    
    gemm2_done = oneapi::mkl::blas::row_major::gemm(queue, conf2.trA, conf2.trB, 
                                                    conf2.m, conf2.n, conf2.k, alpha, G.data, conf2.lda, 
                                                    Q.data, conf2.ldb, beta, global_GQ.data, conf2.ldc);

    gemm1_done.wait();
    gemm2_done.wait();
    gemm3_done = oneapi::mkl::blas::row_major::gemm(queue, conf3.trA, conf3.trB, 
                                                    conf3.m, conf3.n, conf3.k, alpha, global_FP.data, conf3.lda,
                                                    F.data, conf3.ldb, beta, global_FPF.data, conf3.ldc);

    gemm4_done = oneapi::mkl::blas::row_major::gemm(queue, conf4.trA, conf4.trB, 
                                                    conf4.m, conf4.n, conf4.k, alpha, global_GQ.data, conf4.lda,
                                                    G.data, conf4.ldb, beta, propP.data, conf4.ldc);

    gemm3_done.wait();
    gemm4_done.wait();
    axpy_done = oneapi::mkl::blas::row_major::axpy(queue, propP.size, alpha, global_FPF.data, 1, propP.data, 1);
    // axpy_done.wait();
}
#endif

// struct matrix IPIV, scratchpad;
// void oneapi_allocInverse(sycl::queue &queue, struct mkl_matrix &A)
// {
//     if(A.row != A.col){
//         fprintf(stderr, "Error, Not a square Matrix\n");
//         exit(EXIT_FAILURE);
//     }

//     __data_type size_scratch =  oneapi::lapack::getrf_scratchpad_size<__data_type>(queue, A.row, A.row, A.ldw);
//     scratchpad.q = queue;
//     IPIV.q = queue;
//     scratchpad.allocate(size_scratch, 1);
//     IPIV.allocate(N,N);
// }


// void oneapi_inverse(sycl::queue &queue, struct mkl_matrix)
// {
//     sycl::event potri_done, potrf_done;
//     vector<sycl::event> event_list;


//     potrf_done = lapack::getrf(queue, N, N, A, N, IPIV, scratchpad, scratch_size, event_list);
//     potri_done = lapack::getri(queue, N, A, N, IPIV, scratchpad, scratch_size, event_list);
    

// }

#endif /*_MKL_KERNELS_HPP_*/