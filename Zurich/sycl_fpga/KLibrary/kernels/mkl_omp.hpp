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

void omp_CovMatrix(struct mkl_matrix &F, struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q, struct mkl_matrix propP)
{    

    double alpha = 1.0, beta =  0.0f;
    #pragma omp parallel sections shared(F, P, global_FP, global_FPF, alpha, beta, conf1, conf3)
    {

        #pragma omp section
        cblas_dgemm(CblasRowMajor, conf1.trA, conf1.trB, 
                    conf1.m, conf1.n, conf1.k, alpha, F.data, conf1.lda, 
                    P.data, conf1.ldb, beta, global_FP.data, conf1.ldc);

        #pragma omp section
        cblas_dgemm(CblasRowMajor, conf3.trA, conf3.trB, 
                    conf3.m, conf3.n, conf3.k, alpha, global_FP.data, conf3.lda, 
                    F.data, conf3.ldb, beta, global_FPF.data, conf3.ldc);
    }

     #pragma omp parallel sections shared(G, Q, global_GQ, propP, alpha, beta, conf2, conf4)
    {
    #pragma omp section
    cblas_dgemm(CblasRowMajor, conf2.trA, conf2.trB, 
                conf2.m, conf2.n, conf2.k, alpha, G.data, conf2.lda, 
                Q.data, conf2.ldb, beta, global_GQ.data, conf2.ldc);
    
    #pragma omp section
    cblas_dgemm(CblasRowMajor, conf4.trA, conf4.trB, 
                conf4.m, conf4.n, conf4.k, alpha, global_GQ.data, conf4.lda, 
                G.data, conf4.ldb, beta, propP.data, conf4.ldc);
    }
    #pragma omp taskwait
    cblas_daxpy(propP.size, alpha, global_FPF.data, 1, propP.data, 1);

}



//Nao vai rolar, cblas_dgemm eh feito para CPU ---> cuBLAS ou naive eh o mais indicado!
void offload_CovMatrix(struct mkl_matrix &F, struct mkl_matrix &P, struct mkl_matrix &G, struct mkl_matrix &Q, struct mkl_matrix propP)
{
    pragma omp target data map(to:    F.data[0:F.size], P.data[0:P.size],     \
                                      G.data[0:G.size], Q.data[0:Q.size])     \
                           map(alloc: global_FP.data[0:global_FP.size],       \
                                      global_FPF.data[0:global_FPF.size])     \
                           map(from:  propP.data[0:propP.size\])
    {
        #pragma omp dispatch target
        cblas_dgemm(CblasRowMajor, conf1.trA, conf1.trB, 
                    conf1.m, conf1.n, conf1.k, 1.0f, F.data, conf1.lda, 
                    P.data, conf1.ldb, 0.0f, global_FP.data, conf1.ldc);
        
        #pragma omp dispatch target
        cblas_dgemm(CblasRowMajor, conf3.trA, conf3.trB, 
                    conf3.m, conf3.n, conf3.k, 1.0f, global_FP.data, conf3.lda, 
                    F.data, conf3.ldb, 0.0f, global_FPF.data, conf3.ldc)
    
        #pragma omp dispatch target
        cblas_dgemm(CblasRowMajor, conf2.trA, conf2.trB, 
                    conf2.m, conf2.n, conf2.k, 1.0f, G.data, conf2.lda, 
                    Q.data, conf2.ldb, 0.0f, global_GQ.data, conf2.ldc);
        
        #pragma omp dispatch target
        cblas_dgemm(CblasRowMajor, conf4.trA, conf4.trB, 
                    conf4.m, conf4.n, conf4.k, 1.0f, global_GQ.data, conf4.lda, 
                    G.data, conf4.ldb, 0.0f, propP.data, conf4.ldc);

        #pragma omp dispatch target
        cblas_daxpy(propP.size, alpha, global_FPF.data, 1, propP.data, 1);
    }

}



#endif //_MKL_OMP_HPP_