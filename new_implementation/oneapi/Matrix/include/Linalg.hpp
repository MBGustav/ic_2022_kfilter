#ifndef LINALG_NAIVE_H
#define LINALG_NAIVE_H

#include "Matrix.hpp"

#include "common_matrix.hpp"


namespace naive{
    
    // Template is at compile time (ijk model)
    // template<data_t alpha, data_t beta> 
    /*TODO: incluir multiplicação de acordo com a transposição*/
    //B = alpha * AB + beta *C
    template<
        bool TransA = false,  
        bool TransB = false,
        bool TransC = false>
    void gemm(Matrix &A,Matrix &B, Matrix &C, data_t alpha, data_t beta)
    {

        int M = TransA ? A.getCol() : A.getRow();
        int K = TransA ? A.getRow() : A.getCol();
        int N = TransB ? B.getRow() : B.getCol();

        #ifndef JUMP_ASSERTION_COMANDS
        assert(K == (TransB ? B.getCol() : B.getRow()));
        assert(M == (TransC ? C.getCol() : C.getRow())); /*op(A) == op(C)*/
        assert(N == (TransC ? C.getRow() : C.getCol())); /*op(A) == op(B)*/

        // assert(M == (TransC ? C.getCol() : C.getRow()));
        #endif /*JUMP_ASSERTION_COMANDS*/

        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                data_t aux = 0.0f;
                for (int k = 0; k < K; k++) {
                    aux += (TransA ? A.at(k, i) : A.at(i, k)) *
                           (TransB ? B.at(j, k) : B.at(k, j));
                }
                
                C.at(i, j) = alpha * aux + beta *  C.at(i, j);
            }
        }
    }
    // Z = alpha(A * B * C) + beta * Z
    template<
    bool TransA = false, 
    bool TransB = false, 
    bool TransC = false, 
    bool TransZ = false
    >
    void trimatrix_mult(Matrix &A, Matrix &B, Matrix &C, Matrix &Z, data_t alpha, 
                        data_t beta)
    {
        #ifndef JUMP_ASSERTION_COMANDS

        #endif /*JUMP_ASSERTION_COMANDS*/

        Matrix aux_AB(A.getRow(), B.getCol());

        // Op: A * B
        naive::gemm<TransA, TransB>(A, B, aux_AB, 1.0, 0.0);

        // Op: alpha AB * C + beta * Z
        naive::gemm<false, TransC, TransZ>(aux_AB, C, Z, alpha, beta);
    }
//     y <- alpha * A * x + beta * y(vectors:x, y)
    void gemv(Matrix &A,Matrix &x, Matrix &y,const data_t alpha,const data_t beta)
    {
        #ifndef JUMP_ASSERTION_COMANDS
        assert(A.getCol() == x.getSize());
        assert(x.getSize() == y.getSize());
        #endif /*JUMP_ASSERTION_COMANDS*/

        #pragma omp parallel for shared(A, x, y) reduction(+:aux)
        for (int i = 0; i<A.getRow(); i++){
            data_t aux = 0.0f;
            for (int k=0; k<x.getSize(); k++) 
                aux += A.at(i,k) * x.idx(k);
            y.idx(i) = alpha * aux + beta * y.idx(i);
        }
    }
    

    // template<data_t alpha> 
    //B = alpha * A + beta * B
    void axpy(Matrix &A, Matrix &B,const data_t alpha, const data_t beta){
    #ifndef JUMP_ASSERTION_COMANDS
    assert(A.getCol() == B.getCol());
    assert(A.getRow() == B.getRow());
    #endif /*JUMP_ASSERTION_COMANDS*/
    
    #pragma omp parallel for
    for(int idx=0; idx<A.getSize(); idx++)
        B.idx(idx) = alpha * A.idx(idx) + beta * B.idx(idx);
    }

    //B = inv(A)
    inline void inverse(Matrix &A, Matrix &B)
    {
    #ifndef JUMP_ASSERTION_COMANDS
    assert(A.getCol() == B.getCol());
    assert(A.getRow() == B.getRow());
    #endif /*JUMP_ASSERTION_COMANDS*/
    data_t sum, aux; 
    int idx2, k, i, j, mem, flag, cont;
    

    Matrix Aux = A.copy_of();

//---------------------------- Partial pivoting --------------------------------
    data_t *b = (data_t*) malloc(A.getRow()*sizeof(data_t));
    data_t *x = (data_t*) malloc(A.getRow()*sizeof(data_t));
    int *idx = (int*)  malloc(A.getRow()*sizeof(int));
    
    for (k = 0; k<A.getRow(); k++)
        idx[k] = k;

    for (i = 0; i<A.getRow(); i++) {
        j = i;
        idx2 = i;
        // if (pa[ncA*i+j] == 0) {
        if (Aux.at(i, j) == 0) {
            flag = 1;
            for (k = i+1; k<A.getRow(); k++ ) {
                // if (fabsf(pa[ncA*k+j]) >= TINY && flag == 1) {
                    if (fabsf(Aux.at(k, j)) >= TINY && flag == 1) {
                    mem  = idx[i];
                    idx[i] = idx[k];
                    idx[k] = mem;
                    idx2 = k;
                    flag = 0;
                }
            }
            if (flag == 1) {
                for (k = 0; k<A.getRow(); k++) {
                    // if (fabsf(pa[ncA*k+j]) > TINY && fabsf(pa[ncA*i+k]) > TINY) {
                        if (fabsf(Aux.at(k, j)) > TINY && fabsf(Aux.at(i, k)) > TINY) {
                        mem = idx[i];
                        idx[i] = idx[k];
                        idx[k] = mem;
                        idx2 = k;
                        flag = 0;
                    }
                }
            }
            if (idx2 == i){
                printf("\n Singular matrix \n \n");
                A.at(i,j) = TINY;
            }
            for (k = 0; k<A.getRow(); k++){
                mem = A.at(i,k);
                A.at(i,k) = Aux.at(idx2, k);
                Aux.at(idx2, k) = mem;
            }
        }

    }
//------------------- Crout's algorithm for LU Decomposition -------------------
    for (j = 0; j<A.getRow(); j++) {
        for (i = 0; i<A.getRow(); i++) {
            if (i<j | i ==j) {
                sum = Aux.at(i,j);
                for (k = 0; k<i; k++) {
                    sum = sum - Aux.at(i,k)*Aux.at(k,j);
                }
                Aux.at(i,j) = sum;
            }
            if (i > j) {
                sum = Aux.at(i,j);
                for (k = 0; k<j; k++) {
                    sum = sum - Aux.at(i,k)*Aux.at(k,j);
                }
                Aux.at(i,j) = sum/Aux.at(j,j);
            }
        }
    }
//---------------------------- Forward substituion -----------------------------
    for (k = 0; k<A.getRow(); k++) {
        #pragma simd vectorlength(AVX_LEN)
        for (cont = 0; cont<A.getRow(); cont ++ ) {
            b[cont] = 0;
        }
        b[k] = 1;
        for (i = 0; i<A.getRow(); i++) {
            sum = b[i];
            #pragma simd vectorlength(AVX_LEN)
            for (j = 0; j<i; j++) {
                sum = sum - Aux.at(i, j)*x[j];
            }
            x[i] = sum;
        }
//---------------------------- Backward substituion ----------------------------
        for (i=(A.getRow()-1); i>=0; i--) {
            sum = x[i];
            for (j = i+1; j<A.getRow(); j++) {
                sum = sum - Aux.at(i, j)*x[j];
            }
            x[i] = sum/Aux.at(i,i);
        }
        for (cont = 0;  cont<A.getRow(); cont++){
            B.at(cont, idx[k]) = x[cont];
        }
    }

    free(b);
    free(x);
    free(idx);
    }       


    void geodetic2ecef(Matrix& llh, Matrix& xyz){
        double lat = angle2r(llh.idx(0));
        double lon = angle2r(llh.idx(1));
        double h = llh.idx(2);

        // World Geodetic System 1984
        // b = 6356752.31424518; SemiminorAxis [m]
        double a = 6378137; // SemimajorAxis [m]
        double e2 = wg84_Eccentricity*wg84_Eccentricity;

        double RN = wg84_SemimajorAxis / (sqrt((1-e2*(sin(lat)*sin(lat)))));

        xyz.idx(0) = (RN + h)*cos(lat)*cos(lon);    // x
        xyz.idx(1) = (RN + h)*cos(lat)*sin(lon);    // y
        xyz.idx(2) = (RN*(1-e2)+h)*sin(lat);    // z
    }


void ecef2ned(Matrix &pe, Matrix &llh0, Matrix &pt){
    double lat0= angle2r(llh0.idx(0));
    double lon0= angle2r(llh0.idx(1));
    double h0  = llh0.idx(2);

    Matrix RTE(3,3);
    Matrix pe0(3,1); 
    Matrix aux(3,1);

    double slo = sin(lon0);
    double sla = sin(lat0);
    double clo = cos(lon0);
    double cla = cos(lat0);

    RTE.idx(0) = -sla*clo; RTE.idx(1) = -sla*slo; RTE.idx(2) = cla;
    RTE.idx(3) = -slo    ; RTE.idx(4) =  clo    ; RTE.idx(5) = 0;
    RTE.idx(6) = -cla*clo; RTE.idx(7) = -cla*slo; RTE.idx(8) = -sla;

    geodetic2ecef(llh0, pe0);
    // less_m(pe, pe0, aux);//aux = pe - pe0
    // times_m(RTE,aux, pt);

    axpy(pe, pe0, 1.0,-1.0);
    gemm(RTE, pe0, pt, 1.0, 0.0); 
    
}

void geodetic2ned(Matrix& llh, Matrix& llh0, Matrix &ned)
{
    Matrix xyz(3,1);
    geodetic2ecef(llh, xyz); 
    ecef2ned(xyz, llh0, ned);
    // delete_m(xyz);
}

}


#endif /*LINALG_H*/