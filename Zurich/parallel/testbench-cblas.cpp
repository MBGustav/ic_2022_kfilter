#include <stdio.h>
#include <stdlib.h>

// #include "structs.h"
// #include "libmat.h"
// #include "kalman_filter.h"
#include "Timer.hpp"
#include "libmat.h"
#include "kalman_filter.h"
#include "defines.h"
#include "kernels/mkl_cblas.hpp"
#include "Timer_omp.hpp"

int main(int argc, char **argv)
{

    int N, verify;
    if(argc <1)
    {
        printf("Benchmark tester that shows a efficience curve to offload or mkl use\n");
        printf("Here we considere row-major matrices operations");
        printf("How to use: ./testbench-omp  <N> [<verifify]\n");
        exit(EXIT_FAILURE);
    }

    if(argc == 2)
        N = std::max(3,atoi(argv[1])); 

    if(argc == 3)
        verify = atoi(argv[2]);
    else 
        verify = 0;
    
    printf("Testing Time for CovMatrix Prediction\n");
    std::string label = "CovPred," +         /*function*/
                        std::to_string(N) +","+  /*size*/
                        "CPU,"+ /*dev_type*/
    #if defined(PAD_LD)
                        "Y,"+
    #else
                        "N,"+
    #endif
                        "cblas"; /*type_exec*/

    //Matrices decl
    struct mkl_matrix F, G, Q, P, propP;
    G.allocate(N, N);
    Q.allocate(N, N);
    P.allocate(N, N);
    F.allocate(N, N);
    propP.allocate(N,N);

    
    ParallelTimer Timer("testbench_CovMatrix.csv", "device_type,padding,type_execution");
    
    //shuffle matrix
    G.rnd_matrix();
    Q.rnd_matrix();
    P.rnd_matrix();

    //Allocate intermediary matrices
    AllocCovMatrix(F,P,G,Q);
    //warm up..
    cblas_CovMatrix(F, P, G, Q, propP);
    
    G.data[0] = 5.4f;

    //run sample
    Timer.start();
    cblas_CovMatrix(F,P,G, Q, propP);
    Timer.stop(label.c_str());





    if(verify){
        //create serial matrices
        struct matrix *G_serial = matrix_m(G.row,G.col);free(G_serial->elements);
        struct matrix *Q_serial = matrix_m(Q.row,Q.col);free(Q_serial->elements);
        struct matrix *P_serial = matrix_m(P.row,P.col);free(P_serial->elements);
        struct matrix *F_serial = matrix_m(F.row,F.col);free(F_serial->elements);
        struct matrix *propP_omp = matrix_m(propP.row, propP.col);
        struct matrix *exp = matrix_m(propP.row, propP.col);

        //changing data
        G_serial->elements  = G.data_NoPadding();
        Q_serial->elements  = Q.data_NoPadding();
        P_serial->elements  = P.data_NoPadding();
        F_serial->elements  = F.data_NoPadding();
        propP_omp->elements = propP.data_NoPadding();

        
        cov_matrix_serial(F_serial, P_serial, G_serial, Q_serial, exp);
        if(matrix_isEqual(exp, propP_omp))
            printf("Matrices are equal, OK!\n");
        else 
            printf("Matrices are NOT equal, NOT OK!\n");
    }
        


    propP.printMat();
    FreeCovMatrix();





    return 0;
}