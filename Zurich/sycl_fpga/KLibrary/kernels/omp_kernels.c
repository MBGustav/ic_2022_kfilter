#include "omp_kernels.h"
#include "libmat.h"

/*
# DISCLAIMER: in these kernel operations, is considered that
 all these operations are correct in parameters such as n_rows, 
 n_cols and not NULL pointers. Any kind of verification must be done BEFORE.
*/




inline int _idx(matrix *A, int x, int y)
{
    return A->n_col*y + x;
}


inline void omp_matmul(matrix *A, matrix *B, matrix *C)
{
    double *pA, *pB, *pC, aux;
    int nrA, nrB, ncA, ncB, ncC, nrC;
    nrA = A->n_row;
    ncA = A->n_col;
    nrB = B->n_row;
    ncB = B->n_col;
    nrC = C->n_row;
    ncC = C->n_col;
    pA = A->elements;
    pB = B->elements;
    pC = C->elements;
    double* pBt = (double*) malloc(nrB * ncB * sizeof(double));

    int ib = 2;
    int kb = 2;
    double acc00, acc01, acc10, acc11;
    #pragma omp parallel
    {
        #pragma omp for simd reduction(+:acc00, acc01, acc11, acc10)
        for (int ii = 0; ii < nrA; ii += ib){
        for (int kk = 0; kk < ncB; kk += kb)
        {
            for (int j=0; j < nrB; j += 2){
            for (int i = ii; i < ii + ib; i += 2){
                    if (kk == 0)
                        acc00 = acc01 = acc10 = acc11 = 0;
                    else
                    {
                        acc00 = pC[_idx(C, i + 0, j + 0)];
                        acc01 = pC[_idx(C, i + 0, j + 1)];
                        acc10 = pC[_idx(C, i + 1, j + 0)];
                        acc11 = pC[_idx(C, i + 1, j + 1)];
                    }
                    for (int k = kk; k < kk + kb; k++)
                    {
                        acc00 += pB[_idx(B, k, j + 0)] * pA[_idx(A, i + 0, k)];
                        acc01 += pB[_idx(B, k, j + 1)] * pA[_idx(A, i + 0, k)];
                        acc10 += pB[_idx(B, k, j + 0)] * pA[_idx(A, i + 1, k)];
                        acc11 += pB[_idx(B, k, j + 1)] * pA[_idx(A, i + 1, k)];
                    }
                    pC[_idx(C, i + 0,j + 0)] = acc00;
                    pC[_idx(C, i + 0,j + 1)] = acc01;
                    pC[_idx(C, i + 1,j + 0)] = acc10;
                    pC[_idx(C, i + 1,j + 1)] = acc11;
                }
            }
        }
        }
    }
}


inline void omp_MultRotational(matrix *A, matrix *B, matrix *C)
{
/*
#  Takes two matrices and mutiply in the following order:
#  C = A x B x At, where At is transposed
*/

    double *pA, *pB, *pC, aux;
    int nrA, nrB, ncA, ncB, ncC, nrC;
    int i, j, cont;
    nrA = A->n_row;
    ncA = A->n_col;
    nrB = B->n_row;
    ncB = B->n_col;
    nrC = C->n_row;
    ncC = C->n_col;
    pA = A->elements;
    pB = B->elements;
    pC = C->elements;
    double *pAB = (double*) malloc( nrA * ncB * sizeof(double));

    // ** Can be done in Parallel tasks? No!** 
    // First: A x B
    #pragma omp parallel for private(i, j, cont) shared(pAB, pB) reduction(+:aux)
    for (i=0; i<nrA; i++) {
    for (j=0; j<ncB; j++) {
        aux=0;
        for (cont=0; cont<nrB; cont++ )
            aux=aux + pA[_idx(A , i, cont)]*pB[_idx(B, cont, j)];
        pAB[_idx(B, i, j)]= aux;
    }
    }

    //Second: (A x B) x At == C x At
    // #pragma omp parallel for private(i, j) shared(pAB, pA) reduction(+:aux)
    for (i = 0; i<nrA; i++) {
    for (j = 0; j<ncB; j++) {
        aux = 0;
        for (cont=0; cont<nrB; cont++ )
            aux=aux + pAB[ncA*i + cont]*pA[ncB*j + cont];
        pC[ncB*i + j]= aux;
    }
    }
}

/*Similar to omp_MultRotational, but the B matrix is similar to Ident. Matrix*/
inline void omp_MultRotationalv2(matrix *A, matrix *B, matrix *C)
{

    double *pA, *pB, *pC, aux;
    int nrA, nrB, ncA, ncB, ncC, nrC;
    int i, j, cont;
    nrA = A->n_row;
    ncA = A->n_col;
    nrB = B->n_row * B->n_col;
    ncB = B->n_col;
    nrC = C->n_row;
    ncC = C->n_col;
    pA = A->elements;
    pB = B->elements;
    pC = C->elements;
    struct matrix *AB = matrix_m(nrA, ncB);
    double *pAB = AB->elements;


    #pragma omp parallel for
    for (i = 0; i<nrA; i++) {
    for (j = 0; j<ncB; j++) {
        pAB[_idx(AB, i, j)]=  pA[_idx(A, i, j)] * pB[j];
    }
    }
    

     //Second: (A x B) x At == C x At
    #pragma omp parallel for
    for (i = 0; i<nrA; i++) {
    for (j = 0; j<ncB; j++) {
        aux = 0;
        for (cont=0; cont<nrB; cont++ )
            aux += pAB[ncA*i + cont]*pA[ncB*j + cont];
        pC[ncB*i + j]= aux;
    }
    }

}
//Takes three matrices, which: out = A x B + (scal) * C
inline void omp_MultAcc(matrix *A, matrix *B, matrix *C, double scal, matrix *out){
/*

*/
    double *pA, *pB, *pC, *pout, aux;
    int nrA, nrB, ncA, ncB, ncC, nrC;
    int i, j, cont;
    nrA = A->n_row;
    ncA = A->n_col;
    nrB = B->n_row;
    ncB = B->n_col;
    nrC = C->n_row;
    ncC = C->n_col; 
    pA = A->elements;
    pB = B->elements;
    pC = C->elements;
    pout = out->elements;


    #pragma omp parallel for private(i, j, cont) \
                    shared(pA, pB, pC, pout) \
                    reduction(+:aux)
    for (i=0; i<nrA; i++) {
        for (j=0; j<ncB; j++) {
            aux=0;
            for (cont=0; cont<nrB; cont++ )
                aux += pA[ncA*i + cont]*pB[ncB*cont + j] + scal*pC[ncA*i + cont] ;
            pout[ncB*i + j]= aux;
            }
        }
}


#ifdef _OFFLOAD

inline void offload_matmul(matrix *A, matrix *B, matrix *C)
{
    
    double *pA, *pB, *pC, aux;
    int nrA, nrB, ncA, ncB, ncC, nrC, szA, szB, szC;
    int i, j, cont;
    nrA = A->n_row;
    ncA = A->n_col;
    szA = nrA * ncA;
    nrB = B->n_row;
    ncB = B->n_col;
    szB = nrB * ncB;
    nrC = C->n_row;
    ncC = C->n_col; 
    szC = nrC * ncC;
    pA = A->elements;
    pB = B->elements;
    pC = C->elements;

    #pragma omp target data device(0) \
                map(to: pA[0:szA], pB[0:szB]) \
                map(from: pC[0:szC])

    #pragma omp distribute parallel for \
        simd collapse(2) num_threads(NrThread)
    for (i = 0; i<nrA; i++) {/*par*/
    for (j = 0; j<ncB; j++) {/*par*/
        aux = 0;
        for (cont=0; cont<nrB; cont++ ) //sequential
            aux=aux + pA[ncA*i + cont]*pB[ncB*cont + j];
        pC[ncB*i + j]= aux;
    }
    }
}

#endif /*_OFFLOAD*/