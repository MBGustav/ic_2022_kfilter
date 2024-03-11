#ifndef __QUATERNION_HPP__
#define __QUATERNION_HPP__

// #define oneAPIMatrix
// Smaller cost to generate templates
// #if defined ONEAPI_MATRIX
// #include "oneAPIMatrix.hpp"
// typedef oneAPIMatrix std_Matrix;
// #elif defined CUDA_MATRIX
// #include "cudaMatrix.hpp"
// typedef cudastd_Matrix Matrix;
// #else 
#include "Matrix.hpp"
typedef Matrix std_Matrix;
// #endif


#include "KalmanFilter.hpp"
#include "Linalg.hpp"


data_t quatnormalize(std_Matrix &Q);
void crossM_f(std_Matrix &A, std_Matrix &M);
void quat2dcm(const std_Matrix &Q, std_Matrix &dcm);
void quaternion_inv(std_Matrix &ql, std_Matrix &Qout);
void angle2dcm(data_t x,data_t y, data_t z, std_Matrix &dcm);
void quaternion_mul(std_Matrix &ql, std_Matrix &qr, std_Matrix &Qout);
void deltaAngle2quat_f(std_Matrix &deltaAngle, std_Matrix &LGqpHat, std_Matrix &out_quaternion);


void quaternion_mul(std_Matrix &ql, std_Matrix &qr, std_Matrix &Qout)
{
    data_t ql_w = ql.idx(0);
    data_t ql_x = ql.idx(1);
    data_t ql_y = ql.idx(2);
    data_t ql_z = ql.idx(3);

    data_t qr_w = qr.idx(0);
    data_t qr_x = qr.idx(1);
    data_t qr_y = qr.idx(2);
    data_t qr_z = qr.idx(3);

    Qout.idx(0) = ql_w * qr_w - ql_x * qr_x - ql_y * qr_y - ql_z * qr_z;
    Qout.idx(1) = ql_w * qr_x + ql_x * qr_w + ql_y * qr_z - ql_z * qr_y;
    Qout.idx(2) = ql_w * qr_y + ql_y * qr_w + ql_z * qr_x - ql_x * qr_z;
    Qout.idx(3) = ql_w * qr_z + ql_z * qr_w + ql_x * qr_y - ql_y * qr_x;   
}



void quat2dcm(const std_Matrix &Q, std_Matrix &dcm)
{
    data_t n = Q.norm();

    data_t x = Q.idx(0)/(n+TINY);  /*qn(:,1)*/
    data_t y = Q.idx(1)/(n+TINY);  /*qn(:,2)*/
    data_t z = Q.idx(2)/(n+TINY);  /*qn(:,3)*/
    data_t w = Q.idx(3)/(n+TINY);  /*qn(:,4)*/
    data_t x2 = x*x;
    data_t y2 = y*y;
    data_t z2 = z*z;
    data_t w2 = w*w;

    dcm.idx(0) = x2 + y2 - z2 - w2;
    dcm.idx(1) = 2.0*(y*z + x*w);
    dcm.idx(2) = 2.0*(y*w - x*z);

    dcm.idx(3) = 2.0*(y*z - x*w);
    dcm.idx(4) = x2 - y2 + z2 - w2;
    dcm.idx(5) = 2.0*(z*w + x*y);
    
    dcm.idx(6) = 2.0*(y*w + x*z);
    dcm.idx(7) = 2.0*(z*w - x*y);
    dcm.idx(8) = x2 - y2 - z2 + w2;
}


void quaternion_inv(std_Matrix &ql, std_Matrix &Qout)
{

    data_t q0 = ql.idx(0); //x;
    data_t q1 = ql.idx(1); //y; 
    data_t q2 = ql.idx(2); //z; 
    data_t q3 = ql.idx(3); //w;

    data_t norm_q = quatnormalize(ql); //#TODO: preciso de sqrt ? 

    Qout.idx(0) =  q0 / (norm_q + TINY);
    Qout.idx(1) = -q1 / (norm_q + TINY);
    Qout.idx(2) = -q2 / (norm_q + TINY);
    Qout.idx(3) = -q3 / (norm_q + TINY);
}



data_t quatnormalize(std_Matrix &Q)
{
    data_t x = Q.idx(0);
    data_t y = Q.idx(1);
    data_t z = Q.idx(2);
    data_t w = Q.idx(3);
    return x*x + y*y + z*z + w*w;
}

void angle2dcm(data_t x,data_t y, data_t z, std_Matrix &dcm)
{
    //Tomando como default a sequencia zyx para conversao

    
    // gera matrizes Tx,Ty,Tz
    std_Matrix Tx(3, 3);
    std_Matrix Ty(3, 3);
    std_Matrix Tz(3, 3);
    std_Matrix dcm_aux(3,3);
    // [1  0 0  ; 0  c s  ; 0 -s c ]
    Tx.at(0,0) = 1; Tx.at(1, 0) = 0;       Tx.at(2, 0) = 0;
    Tx.at(0,1) = 0; Tx.at(1, 1) = cos_d(z);Tx.at(2, 1) = sin_d(z);
    Tx.at(0,2) = 0; Tx.at(1, 2) =-sin_d(z);Tx.at(2, 2) = cos_d(z);

    // c 0 -s ; 0 1 0 ; s 0 c ]
    Ty.at(0,0) = cos_d(y);Ty.at(1,0) = 0; Ty.at(2, 0) =-sin_d(y);
    Ty.at(0,1) = 0;       Ty.at(1,1) = 1; Ty.at(2, 1) = 0; 
    Ty.at(0,2) = sin_d(y);Ty.at(1,2) = 0; Ty.at(2, 2) = cos_d(y);

    // c s 0 ; -s c 0 ; 0 0 1 ] 
    Tz.at(0,0) = cos_d(x);Tz.at(1,0) = sin_d(x); Tz.at(2,0) = 0.0f;
    Tz.at(0,1) =-sin_d(x);Tz.at(1,1) = cos_d(x); Tz.at(2,1) = 0.0f;
    Tz.at(0,2) = 0;       Tz.at(1,2) = 0;        Tz.at(2,2) = 1.0f;


    //dcm =  Tx * Ty * Tz 
    naive::trimatrix_mult(Tx,Ty, Tz, dcm, 1.0, 0.0);
}

void deltaAngle2quat_f(std_Matrix &deltaAngle, std_Matrix &LGqpHat, std_Matrix &out_quaternion)
{
   
    std_Matrix deltaqvector(1,3);
    std_Matrix deltaqaux(1, 1); // multiply: At(1x3) * At(3x1) = B(1x1)
    std_Matrix deltaq(4, 1);
    
    for(int i=0; i < deltaqvector.getSize(); i++)
        deltaqvector.idx(i) = deltaAngle.idx(i) * 0.5;
    
    // deltaqaux = (deltaAngle/2).T*(deltaAngle/2);
    naive::gemm<true,false>(deltaqvector, deltaqvector, deltaqaux, 1.0, 0.0);

    if(LessThanOrEq(deltaqaux.idx(0), 1.0)){
        // deltaq = [sqrt(1 - deltaqaux); deltaqvector];
        deltaq.idx(0) = sqrt(1 - deltaqaux.idx(0));
        deltaq.fill_matrix(deltaqvector, 1);
    }else{
        // deltaq = [1/sqrt(1 + deltaqaux); deltaqvector/sqrt(1 + deltaqaux)];
        deltaq.idx(0) = 1.0f / sqrt(1.0f + deltaqaux.idx(0));
        for(int i = 0; i < deltaqvector.getSize(); i++){
            deltaq.idx(1+i) = deltaqvector.idx(1) / sqrt(1.0f + deltaqaux.idx(0));
        }
    }

    quaternion_mul(deltaq, LGqpHat, out_quaternion);
}


void quat2dcm(Matrix &Q, Matrix &dcm){

    double n = Q.norm(); 
    double x ,y, z, w,
           x2,y2,z2,w2;

    x = Q.idx(0)/(n+TINY);  /*qn(:,1)*/
    y = Q.idx(1)/(n+TINY);  /*qn(:,2)*/
    z = Q.idx(2)/(n+TINY);  /*qn(:,3)*/
    w = Q.idx(3)/(n+TINY);  /*qn(:,4)*/
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    w2 = w*w;

    dcm.at(0,0) = x2 + y2 - z2 - w2;
    dcm.at(0,1) = 2.0 * (y*z + x*w);
    dcm.at(0,2) = 2.0 * (y*w - x*z);
    dcm.at(1,0) = 2.0 * (y*z - x*w);
    dcm.at(1,1) = x2 - y2 + z2 - w2;
    dcm.at(1,2) = 2.0 * (z*w + x*y);
    dcm.at(2,0) = 2.0 * (y*w + x*z);
    dcm.at(2,1) = 2.0 * (z*w - x*y);
    dcm.at(2,2) = x2 - y2 - z2 + w2;

}


void crossM_f(std_Matrix &A, std_Matrix &M)
{
/*
% Descrição: Retorna a matriz antissimétrica M de um vetor da velocidade angular
% Entrada: vetor em R3
% Saída: Matriz antissimétrica de produto
*/
        #ifndef JUMP_ASSERTION_COMANDS
        assert(A.getSize() == 3);
        assert(M.getCol() == 3 && M.getRow() == 3);
        
        #endif /*JUMP_ASSERTION_COMANDS*/
        
        double v1,v2,v3;
        v1 = A.idx(0);
        v2 = A.idx(1);
        v3 = A.idx(2);

        M.at(0, 0) =  0; M.at(0, 1) =-v3; M.at(0, 2) = v2;
        M.at(1, 0) = v3; M.at(1, 1) =  0; M.at(1, 2) =-v1;
        M.at(2, 0) =-v2; M.at(2, 1) = v1; M.at(2, 2) =  0;
}




#endif /*__QUATERNION_HPP__*/