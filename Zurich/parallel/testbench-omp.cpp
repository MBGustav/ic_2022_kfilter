#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


// #include "structs.h"
// #include "libmat.h"
// #include "kalman_filter.h"
#include "Timer.hpp"
#include "libmat.h"
#include "kalman_filter.h"
#include "defines.h"
#include "kernels/mkl_omp.hpp"
#include "Timer_omp.hpp"

enum class DeviceType {
    CPU,
    GPU
};

struct CommandLineOptions {
    int N;
    bool verify;
    DeviceType deviceType;
};



void print_banner()
{
    printf("Benchmark tester that shows a efficience curve to offload or mkl use\n");
    printf("Here we consider row-major matrices operations");
    printf("How to use: ./testbench-omp  <N> [<verifify] [device_specify]\n");
    printf("[verify]: Set a compare the results on CPU(optional).\n");
    printf("\t0 - No(default argument); \n\t1 -Yes;\n");
    printf("[device_specify]: Set a parameter to run on a specific device(optional).\n");
    printf("\t0 - cpu device(default argument); \n\t1 - gpu device;\n");

    
    exit(EXIT_FAILURE);
}




CommandLineOptions parseCommandLine(int argc, char **argv) {
    CommandLineOptions options;

    if (argc < 2 || argc > 4) print_banner();
    

    options.N = std::max(3, atoi(argv[1]));

    if (argc >= 3) {
        options.verify = (atoi(argv[2]) == 1);
    } else {
        options.verify = false;
    }

    if (argc == 4) {
        int deviceType = atoi(argv[3]);
        if (deviceType == 0) {
            options.deviceType = DeviceType::CPU;
        } else if (deviceType == 1) {
            options.deviceType = DeviceType::GPU;
        } else {
            std::cerr << "Invalid device type.\n";
            exit(EXIT_FAILURE);
        }
    } else {
        options.deviceType = DeviceType::CPU;
    }

    return options;
}

int main(int argc, char **argv)
{

    CommandLineOptions options = parseCommandLine(argc, argv);
    int N = options.N;
    printf("Testing Time for CovMatrix Prediction\n");
    printf("Range Defined for N = %i\n\n", N);

    
    #if defined(PAD_LD) 
    std::cout << "The test-case is considering padding\n";
    #else
    std::cout << "The test-case is NOT considering padding\n";
    #endif
    
    //Matrices decl
    struct mkl_matrix F, G, Q, P, propP, propP_cblas;
    
    G.allocate(N, N);
    Q.allocate(N, N);
    F.allocate(N, N);
    P.allocate(N, N);
    propP.allocate(N,N);
    propP_cblas.allocate(N,N);
    
    Timer T;
    // function,device_type,padding,type_execution,time_us
    ParallelTimer Timer("testbench_CovMatrix.csv", "size,device_type,padding,type_execution");
    // timer_setconf("testbench_CovMatrix.csv", &T, "size,device_type,padding,type_execution");
    //shuffle matrix
    G.rnd_matrix();
    Q.rnd_matrix();
    F.rnd_matrix();
    P.rnd_matrix();
    
    std::string label = "CovPred," +         /*function*/
                    std::to_string(N) +","+  /*size*/
        (options.deviceType == DeviceType::GPU ? "GPU,":"CPU,")+ /*dev_type*/
    #if defined(PAD_LD)
                "Y,"+
    #else
                "N,"+
    #endif
                "openMP"; /*type_exec*/


    std::string label_cblas = "CovPred,"+/*function*/
                std::to_string(N) +","+  /*size*/
        (options.deviceType == DeviceType::GPU ? "GPU,":"CPU,")+ /*dev_type*/
    #if defined(PAD_LD)
                "Y,"+
    #else
                "N,"+
    #endif
                "cblas"; /*type_exec*/


    //Allocate intermediary matrices
    omp_AllocCovMatrix(F,P,G,Q);
    if(options.deviceType == DeviceType::CPU){
        //warm up..
        std::cout << "CPU execution: \n";
        omp_CovMatrix(F, P, G, Q, propP);
        blas_CovMatrix(F,P, G, Q propP_cblas);
        G.data[0] = 5.4f;

        //run sample
        // timer_start(&T);
        Timer.start();
        omp_CovMatrix(F,P,G, Q, propP);
        Timer.stop(label.c_str());
        // timer_stop(&T, label.c_str());

    }else{
        std::cout << "GPU execution: \n";

        omp_CovMatrix_offload(F,P,G, Q, propP);
        G.data[0] = 5.4f;
        
        // timer_start(&T);

        Timer.start();
        omp_CovMatrix_offload(F,P,G, Q, propP);
        Timer.stop(label.c_str());

        // timer_stop(&T, label.c_str());
    }


    // if(verify){
        //create serial matrices
        // struct matrix *G_serial = matrix_m(G.row,G.c ol);free(G_serial->elements);
        // struct matrix *Q_serial = matrix_m(Q.row,Q.col);free(Q_serial->elements);
        // struct matrix *P_serial = matrix_m(P.row,P.col);free(P_serial->elements);
        // struct matrix *F_serial = matrix_m(F.row,F.col);free(F_serial->elements);
        // struct matrix *propP_omp = matrix_m(propP.row, propP.col);
        // struct matrix *exp = matrix_m(propP.row, propP.col);

        // //changing data
        // G_serial->elements  = G.data_NoPadding();
        // Q_serial->elements  = Q.data_NoPadding();
        // P_serial->elements  = P.data_NoPadding();
        // F_serial->elements  = F.data_NoPadding();
        // propP_omp->elements = propP.data_NoPadding();

        
        // cov_matrix_serial(F_serial, P_serial, G_serial, Q_serial, exp);
        // if(matrix_isEqual(exp, propP_omp))
        //     printf("Matrices are equal, OK!\n");
        // else 
        //     printf("Matrices are NOT equal, NOT OK!\n");
    // }

    propP.printMat();
    omp_FreeCovMatrix();





    return 0;
}