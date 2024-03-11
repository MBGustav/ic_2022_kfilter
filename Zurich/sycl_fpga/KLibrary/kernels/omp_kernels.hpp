#ifndef OMP_KERNELS_H_
#define OMP_KERNELS_H_

#include <stdlib.h>
#include <stdio.h>
#include "defines.h"
#include "structs.h"

#define NrThread (1<<7)
#define BlockDim (16)
#define BlockRow (512)


int _idx(matrix *A, int x, int y);

void omp_matmul(matrix *A, matrix *B, matrix *C);

void omp_MultRotational(matrix *A, matrix *B, matrix *C);

void omp_MultRotationalv2(matrix *A, matrix *B, matrix *C);

void omp_MultAcc(matrix *A, matrix *B, matrix *C, double scal, matrix *out);


/* Functions To Offload (GPU or CPU)*/
#ifdef _OFFLOAD
void offload_matmul(matrix *A, matrix *B, matrix *C)
#endif /*_OFFLOAD*/




#endif //OMP_KERNELS_H_