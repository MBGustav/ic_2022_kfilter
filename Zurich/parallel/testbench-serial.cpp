#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "structs.h"
#include "Timer.hpp"
#include "libmat.h"
#include "kalman_filter.h"
#include "defines.h"

#define DATA_PRECISION (2)
// Here we create a sample to not be influenced by external compile flags being used in MKL Sections
int main(int argc, char **argv)
{
    int N, verify;
    if(argc == 1)
    {
        printf("Benchmark tester that shows a efficience curve to offload or mkl use\n");
        printf("Here we considere row-major matrices operations\n");
        printf("How to use: ./testbench-omp  <N> [<verifify]\n");
        exit(EXIT_FAILURE);
    }
    if(argc == 2)
        N = std::max(3,atoi(argv[1])); 

    if(argc == 3)
        verify = atoi(argv[2]);
    else 
        verify = 0;
    
    //Matrices decl
    struct matrix *F, *G, *Q, *P, *propP;
    G = rnd_matrix(N,N);
    Q = rnd_matrix(N,N);
    P = rnd_matrix(N,N);
    F = rnd_matrix(N,N);
    propP = matrix_m(N,N);
    
    Timer T;
    timer_setconf("testbench_CovMatrix.csv", &T, "size,device_type,padding,type_execution");

    printf("Testing Time for CovMatrix Prediction\n");
    std::string label = "CovPred," +             /*function*/
                        std::to_string(N) +","+  /*size*/
                        "CPU,"+                  /*dev_type*/
                        "N,"+                    /*padding*/
                        "serial";               /*type_exec*/
    //run sample
    timer_start(&T);
    cov_matrix_serial(F, P, G, Q, propP);

                         
    timer_stop(&T, label.c_str());
    
    printf("Printing top left:\n");
    for (int i = 0; i< std::min(6, propP->n_row); i++)  {        
    for (int j = 0; j< std::min(6, propP->n_col); j++)  {
                printf("  %08.4lf", propP->elements[i*propP->n_col + j]);                                     
            }
        printf(" \n");               
    }
    printf("\n");
    return 0;
}