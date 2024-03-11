#include <stdio.h>
#include <stdlib.h>

#include <sycl/sycl.hpp>

#include "structs.h"
#include "libmat.h"
#include "kalman_filter.h"
#include "Timer.hpp"
#include "libmat.h"
#include "kalman_filter.h"
#include "defines.h"
#include "kernels/mkl_kernels.hpp"

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
    int N = 15, verify = 0;

    sycl::queue main_queue;
    CommandLineOptions options = parseCommandLine(argc, argv);
    N = options.N;
    std::cout << "Testing Time for CovMatrix Prediction using oneAPI\n";
    try 
    {

        if(options.deviceType == DeviceType::GPU)
            main_queue = sycl::queue(sycl::gpu_selector_v);
        else
            main_queue = sycl::queue(sycl::cpu_selector_v);
                                
        std::cout << "Device:                 " << main_queue.get_device().get_info<sycl::info::device::name>() << std::endl
                //   << "Core/EU count:          " << main_queue.get_device().get_info<sycl::info::device::max_compute_units>()  << std::endl
                  << "Maximum clock frequency:" << main_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << " MHz" << std::endl
                  << "Range Defined:          " << options.N << std::endl
                  << "Data Type:              " << (DATA_PRECISION == 1 ? "DOUBLE" : "FLOAT") <<std::endl
#if defined(PAD_LD) 
                  << "Padding                 " << "Yes" << std::endl;
#else
                  << "Padding                 " << "No" << std::endl;
#endif

    } catch (sycl::exception const& e) {
        // Handle exceptions and display error messages
        std::cerr << "Error: " << e.what() << std::endl;
        std::terminate();
    }
    

    auto exception_handler = [] (sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception during Cov Matrix Execution:\n"
                << e.what() << std::endl;
            }
        }
    };
    // create execution queue and buffers of matrix data


    //Matrices decl
    struct mkl_matrix F(main_queue), G(main_queue), Q(main_queue), P(main_queue), propP(main_queue);
    
    G.allocate(N, N);
    Q.allocate(N, N);
    F.allocate(N, N);
    P.allocate(N, N);
    propP.allocate(N,N);
    
    Timer T;

    timer_setconf("testbench_CovMatrix.csv", &T, "size,device_type,padding,type_execution");

    //shuffle matrix
    G.rnd_matrix();
    Q.rnd_matrix();
    F.rnd_matrix();
    P.rnd_matrix();

    //Allocate intermediary matrices
    oneapi_AllocCovMatrix(F,P,G,Q);
    //warm up..
    oneapi_CovMatrix(main_queue,F, P, G, Q, propP);
    std::string label = "CovPred," +         /*function*/
                    std::to_string(N) +","+  /*size*/
        (options.deviceType == DeviceType::GPU ? "GPU,":"CPU,")+ /*dev_type*/
    #if defined(PAD_LD)
                    "Y,"+
    #else
                    "N,"+
    #endif
                    "oneMKL";
    //run sample
    timer_start(&T);
    oneapi_CovMatrix(main_queue,F,P,G, Q, propP);
    timer_stop(&T,label.c_str());


    if(verify)
    {
        std::vector<__data_type> inG(G.row * G.col);
        std::vector<__data_type> inQ(G.row * G.col);
        std::vector<__data_type> inP(G.row * G.col);
        std::vector<__data_type> inF(G.row * G.col);
        
        // std::vector<__data_type> GQ = matmul(G,Q);
        // std::vector<__data_type> FP = matmul(F,P);
    }
    // if(verify){
    //     //create serial matrices
    //     struct matrix *G_serial = matrix_m(G.row,G.col);free(G_serial->elements);
    //     struct matrix *Q_serial = matrix_m(Q.row,Q.col);free(Q_serial->elements);
    //     struct matrix *P_serial = matrix_m(P.row,P.col);free(P_serial->elements);
    //     struct matrix *F_serial = matrix_m(F.row,F.col);free(F_serial->elements);
    //     struct matrix *propP_oneapi = matrix_m(propP.row, propP.col);
    //     struct matrix *exp = matrix_m(propP.row, propP.col);

    //     //changing data
    //     G_serial->elements  = G.data_NoPadding();
    //     Q_serial->elements  = Q.data_NoPadding();
    //     P_serial->elements  = P.data_NoPadding();
    //     F_serial->elements  = F.data_NoPadding();
    //     propP_oneapi->elements = propP.data_NoPadding();

        
    //     cov_matrix_serial(F_serial, P_serial, G_serial, Q_serial, exp);
    //     if(matrix_isEqual(exp, propP_oneapi))
    //         printf("Matrices are equal, OK!\n");
    //     else 
    //         printf("Matrices are NOT equal, NOT OK!\n");
    // }

    propP.printMat();
    oneapi_FreeCovMatrix();





    return 0;
}