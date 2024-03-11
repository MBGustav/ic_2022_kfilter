#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#include "libmat.h"
#include "Timer.hpp"



#ifndef _SERIAL
#include "kernels/mkl_kernels.hpp"
#ifdef MKL_ONEAPI
#include <sycl/sycl.hpp>
#include "oneapi/mkl.hpp"
#else
#include "mkl.h"
#endif
#endif

Timer T;


struct matrix* cov_matrix_serial(struct matrix* F, struct matrix* P,
                                 struct matrix* G, struct matrix* Q)
{
    int Gr = G->n_row;
    int Gc = G->n_col;
    int Qr = Q->n_row;
    int Qc = Q->n_col;
    int Fr = F->n_row;
    int Fc = F->n_col;
    int Pr = P->n_row;
    int Pc = P->n_col;
    
    struct matrix* Gt = matrix_m(Gc, Gr);
    struct matrix* Ft = matrix_m(Fc, Fr);
    struct matrix* FP = matrix_m(Fr, Pc);
    struct matrix* GQ = matrix_m(Gr, Qc);
    struct matrix* FPF = matrix_m(Fr, Pc);
    struct matrix* GQG = matrix_m(Fr, Pc);

    struct matrix* propP = matrix_m(Fr, Pc);
    
    // FPF  = F * P * F'
    times_m(F, P, FP);      
    transp_m(F, Ft);
    times_m(FP, Ft, FPF);
    
    // GQG = G x Q x G'
    times_m(G, Q, GQ);
    transp_m(G, Gt);
    times_m(GQ, Gt, GQG);

    sum_m(FPF, GQG, propP); 

    delete_m (Gt);
    delete_m (Ft);
    delete_m (FP);
    delete_m (GQ);
    delete_m (FPF);
    delete_m (GQG);

    return propP;
}

void cov_matrix_test(int m)
{
    
    struct matrix *F, *serial_F;
    struct matrix *P, *serial_P;
    struct matrix *G, *serial_G;
    struct matrix *Q, *serial_Q;
    struct matrix *propP, *serial_propP;

    char simd_name[50], cblas_name[50];

    int i,j, nr = 4, nc = 4; 
    // sprintf(simd_name, "simd-vector,%d",m);
    #ifndef _SERIAL
    #ifdef MKL_ONEAPI

    sycl::queue q; //select any device
    sprintf(cblas_name, "oneapi_mkl,%d",m);
    F = rnd_matrix_mkl(q ,m, m);
    P = rnd_matrix_mkl(q ,m, m);
    G = rnd_matrix_mkl(q ,m, m);
    Q = rnd_matrix_mkl(q ,m, m);
    propP = matrix_m_mkl(q, m, m);
    // printf("F\n");
    // display_matrix_mkl(F);
    // printf("P\n");
    // display_matrix_mkl(P);
    // printf("G\n");
    // display_matrix_mkl(G);
    // printf("Q\n");
    // display_matrix_mkl(Q);

    alloc_CovMatrixPredAux(q, F, P, G, Q);
    timer_start(&T);
        oneapi_cov_matrix_pred(q, F, P, G, Q, propP);
    timer_stop(&T, cblas_name);
    
    free_CovMatrixPredAux(q);

    printf("Printing first block of 4 of propP:\n");
    for (i=0; i<nr; i++){
        for (j=0; j<nc; j++)  {
            printf("  %05.4lf", propP->elements[i*propP->n_col + j]);
        }
            printf(" ;\n"); 
    }

    
    delete_m_mkl(q, F);
    delete_m_mkl(q, P);
    delete_m_mkl(q, G);
    delete_m_mkl(q, Q);
    delete_m_mkl(q, propP);
#else
    
    sprintf(cblas_name, "cblas,%d",m);

    F = rnd_matrix_mkl(m, m);
    P = rnd_matrix_mkl(m, m);
    G = rnd_matrix_mkl(m, m);
    Q = rnd_matrix_mkl(m, m);
    propP = matrix_m_mkl(m, m);

    alloc_CovMatrixPredAux(F, P, G, Q);
    // printf("P\n");
    // display_matrix_mkl(P);
    // printf("F\n");
    // display_matrix_mkl(F);
    // printf("G\n");
    // display_matrix_mkl(G);
    // printf("Q\n");
    // display_matrix_mkl(Q);
    

    timer_start(&T);
        mkl_cov_pred(F, P, G, Q, propP);
    timer_stop(&T, cblas_name);
    
    free_CovMatrixPredAux();



    printf("Printing first block of 4 of propP:\n");
    for (i=0; i<nr; i++){
        for (j=0; j<nc; j++)  {
            printf("  %05.4lf", propP->elements[i*propP->n_row + j]);                                            
        }
            printf(" ;\n"); 
    }

    delete_m_mkl(F);
    delete_m_mkl(P);
    delete_m_mkl(G);
    delete_m_mkl(Q);
    delete_m_mkl(propP);

#endif
#else //from serial
    F = rnd_matrix(m, m);
    P = rnd_matrix(m, m);
    G = rnd_matrix(m, m);
    Q = rnd_matrix(m, m);
    // propP = matrix_m(m, m);

    
    printf("P\n");
    print_m(P,1);
    printf("F\n");
    print_m(F,1);
    printf("G\n");
    print_m(G,1);
    printf("Q\n");
    print_m(Q,1);
    
    sprintf(simd_name, "simd,%d",m);


    timer_start(&T);
    propP = cov_matrix_serial(F, P, G, Q);
    timer_stop(&T, simd_name);
    



    printf("Printing first block of 4 of propP:\n");
    for (i=0; i<nr; i++){
        for (j=0; j<nc; j++)  {
            printf("  %05.4lf", propP->elements[i*propP->n_row + j]);                                            
        }
            printf(" ;\n"); 
    }


    delete_m(F);
    delete_m(P);
    delete_m(G);
    delete_m(Q);
    delete_m(propP);

#endif
}


int main(int argc, char **argv){

    timer_setconf("compare_cvm.csv", &T, "function,size");
    srand(777);
    int i;
    if(argc != 2)
        exit(EXIT_FAILURE);

    i = atoi(argv[1]);

    
    cov_matrix_test(i);

    printf("Finished\n");
    return 0;
}