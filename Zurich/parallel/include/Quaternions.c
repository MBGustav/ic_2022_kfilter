
#include "Quaternions.h"
#include "libmat.h"
#include "defines.h"
#include "Initialization.h"


void quaternion_mul(struct matrix *ql, struct matrix *qr, struct matrix *Qout){

    if(ql == NULL || qr == NULL || Qout == NULL)
        fprintf(stderr, "Could not multiply!\n");

    
    __data_type ql_w = ql->elements[0];
    __data_type ql_x = ql->elements[1];
    __data_type ql_y = ql->elements[2];
    __data_type ql_z = ql->elements[3];

    __data_type qr_w = qr->elements[0];
    __data_type qr_x = qr->elements[1];
    __data_type qr_y = qr->elements[2];
    __data_type qr_z = qr->elements[3];

    Qout->elements[0] = ql_w * qr_w - ql_x * qr_x - ql_y * qr_y - ql_z * qr_z;
    Qout->elements[1] = ql_w * qr_x + ql_x * qr_w + ql_y * qr_z - ql_z * qr_y;
    Qout->elements[2] = ql_w * qr_y + ql_y * qr_w + ql_z * qr_x - ql_x * qr_z;
    Qout->elements[3] = ql_w * qr_z + ql_z * qr_w + ql_x * qr_y - ql_y * qr_x;   
}

void quat2dcm(struct matrix *Q, struct matrix* dcm){

    __data_type n = norm_m(Q); 
    __data_type x ,y, z, w,
           x2,y2,z2,w2;

    x = Q->elements[0]/(n+TINY);  /*qn(:,1)*/
    y = Q->elements[1]/(n+TINY);  /*qn(:,2)*/
    z = Q->elements[2]/(n+TINY);  /*qn(:,3)*/
    w = Q->elements[3]/(n+TINY);  /*qn(:,4)*/
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    w2 = w*w;

    dcm->elements[0] = x2 + y2 - z2 - w2;
    dcm->elements[1] = 2.0*(y*z + x*w);
    dcm->elements[2] = 2.0*(y*w - x*z);

    dcm->elements[3] = 2.0*(y*z - x*w);
    dcm->elements[4] = x2 - y2 + z2 - w2;
    dcm->elements[5] = 2.0*(z*w + x*y);
    
    dcm->elements[6] = 2.0*(y*w + x*z);
    dcm->elements[7] = 2.0*(z*w - x*y);
    dcm->elements[8] = x2 - y2 - z2 + w2;


    // dcm->elements[3*0 + 0] = x*x + y*y - z*z - w*w;
    // dcm->elements[3*0 + 1] = 2*(y*z + x*w);
    // dcm->elements[3*0 + 2] = 2*(y*w - x*z);
    // dcm->elements[3*1 + 0] = 2*(y*z - x*w);
    // dcm->elements[3*1 + 1] = x*x - y*y + z*z - w*w;
    // dcm->elements[3*1 + 2] = 2*(z*w + x*y);
    // dcm->elements[3*2 + 0] = 2*(y*w + x*z);
    // dcm->elements[3*2 + 1] = 2*(z*w - x*y);
    // dcm->elements[3*2 + 2] = x*x - y*y - z*z + w*w;
}



void quaternion_inv(struct matrix *ql, struct matrix *Qout){

    if(ql == NULL || Qout == NULL){
        fprintf(stderr, "Could not multiply!\n");
        exit(1);
    }
    __data_type q0 = ql->elements[0]; //x;
    __data_type q1 = ql->elements[1]; //y; 
    __data_type q2 = ql->elements[2]; //z; 
    __data_type q3 = ql->elements[3]; //w;

    __data_type norm_q = quatnormalize(ql); //#TODO: preciso de sqrt ? 

    Qout->elements[0] =  q0 / (norm_q + TINY);
    Qout->elements[1] = -q1 / (norm_q + TINY);
    Qout->elements[2] = -q2 / (norm_q + TINY);
    Qout->elements[3] = -q3 / (norm_q + TINY);
}


inline __data_type quatnormalize(struct matrix *Q){
    if(Q == NULL){
        fprintf(stderr, "Could not QuatNormalize!\n");
        exit(1);
    }
    __data_type x = Q->elements[0];
    __data_type y = Q->elements[1];
    __data_type z = Q->elements[2];
    __data_type w = Q->elements[3];
    return x*x + y*y + z*z + w*w;
}

void quatDisplay(struct matrix *Q){
    if(Q == NULL ){
        fprintf(stderr, "Could not Display!\n");
        exit(1);
    }
    __data_type x = Q->elements[0];
    __data_type y = Q->elements[1];
    __data_type z = Q->elements[2];
    __data_type w = Q->elements[3];
    printf("\nQ = (%.3f, %.3f, %.3f, %.3f ) \n", x, y, z, w);
    
}

void angle2dcm(__data_type x,__data_type y, __data_type z, struct matrix* dcm){
    //Tomando como default a sequencia zyx para conversao

    if(dcm == NULL) alloc_error();
    // gera matrizes Tx,Ty,Tz
    struct matrix * Tx = matrix_m(3, 3);
    struct matrix * Ty = matrix_m(3, 3);
    struct matrix * Tz = matrix_m(3, 3);
    struct matrix *dcm_aux = matrix_m(3,3);
    // [1  0 0  ; 0  c s  ; 0 -s c ]
    Tx->elements[0] = 1; Tx->elements[1] = 0;      Tx->elements[2] = 0;
    Tx->elements[3] = 0; Tx->elements[4] = cos_d(z);Tx->elements[5] = sin_d(z);
    Tx->elements[6] = 0; Tx->elements[7] =-sin_d(z);Tx->elements[8] = cos_d(z);

    // c 0 -s ; 0 1 0 ; s 0 c ]
    Ty->elements[0] = cos_d(y);Ty->elements[1] = 0; Ty->elements[2] =-sin_d(y);
    Ty->elements[3] = 0;      Ty->elements[4] = 1; Ty->elements[5] = 0; 
    Ty->elements[6] = sin_d(y);Ty->elements[7] = 0; Ty->elements[8] = cos_d(y);

    // c s 0 ; -s c 0 ; 0 0 1 ] 
    Tz->elements[0] = cos_d(x);Tz->elements[1] = sin_d(x); Tz->elements[2] = 0;
    Tz->elements[3] =-sin_d(x);Tz->elements[4] = cos_d(x); Tz->elements[5] = 0;
    Tz->elements[6] = 0;      Tz->elements[7] = 0;       Tz->elements[8] = 1;


    //dcm =  Tx * Ty * Tz 
    times_m(Tx, Ty, dcm_aux);
    times_m(dcm_aux, Tz, dcm);

    delete_m(Tx);
    delete_m(Ty);
    delete_m(Tz);
    delete_m(dcm_aux);

}
void angle2quat(__data_type z,__data_type y, __data_type x, struct matrix* quat){
    
    //theta = z, y, x
    __data_type c1 = cos(z/2); 
    __data_type s1 = sin(z/2);
    __data_type c2 = cos(y/2);
    __data_type s2 = sin(y/2);
    __data_type s3 = sin(x/2);
    __data_type c3 = cos(x/2);

    quat->elements[0] = c1 * c2 * c3 + s1 * s2 * s3;
    quat->elements[1] = c1 * c2 * s3 - s1 * s2 * c3;
    quat->elements[2] = c1 * s2 * c3 + s1 * c2 * s3;
    quat->elements[3] = s1 * c2 * c3 - c1 * s2 * s3;    
}


void deltaAngle2quat_f(struct matrix* deltaAngle, struct matrix* LGqpHat, struct matrix* out_quaternion){

// deltaAngle (3x1) - angle(rpy)
// LGqpHat(4x1) - quaternion
// out_quaternion(4x1)

    if(!deltaAngle || !LGqpHat || !out_quaternion) alloc_error();
    
    struct matrix* deltaq_vector = matrix_m(1,3);
    struct matrix* deltaq_vector_trp = matrix_m(3,1);
    struct matrix* deltaq_aux = matrix_m(1, 1); // multiply: At(1x3) * At(3x1) = B(1x1)

    struct matrix* delta_q = matrix_m(1, 4);
    // % parte imaginário quatérnio de erro de atitude
    timesc_m(0.5, deltaAngle, deltaq_vector);
    transp_m(deltaq_vector, deltaq_vector_trp);

    // % Parâmetro usado para definir o cálculo do quatérnio do erro de atitude
    times_m(deltaq_vector, deltaq_vector_trp, deltaq_aux);



    __data_type deltaqaux = deltaq_aux->elements[0]; 
    if(LessThanOrEq(deltaqaux, (__data_type) 1.0)){
        delta_q->elements[0] = sqrt(1-deltaqaux);
        getpart_mat(deltaq_vector, 0, 1, delta_q);
        
    } else {
        __data_type rsqrt = 1.0 / (sqrt(1+deltaqaux));
        delta_q->elements[0] = rsqrt;

        timesc_m(rsqrt, deltaq_vector, deltaq_vector);//delta_q[1:4]
        getpart_mat(deltaq_vector, 0, 1, delta_q);
    }
    quaternion_mul(delta_q, LGqpHat, out_quaternion);

    delete_m( deltaq_vector) ;
    delete_m( deltaq_vector_trp) ;
    delete_m( deltaq_aux) ;
    delete_m( delta_q);
}


void quat2angle(struct matrix* Q, struct matrix *xyz){
    
    if(Q ==NULL || xyz ==NULL) alloc_error();
    __data_type *pQ    =   Q->elements;
    __data_type *p_xyz = xyz->elements;
    __data_type norm_q = norm_m(Q);
    __data_type qn1, qn2,qn3,qn4;
    __data_type x,y,z;

    qn1 = pQ[0]/(norm_q);
    qn2 = pQ[1]/(norm_q);
    qn3 = pQ[2]/(norm_q);
    qn4 = pQ[3]/(norm_q);

    z = atan2( 2*(qn2*qn3 + qn1*qn4), qn1*qn1 + qn2*qn2 - qn3*qn3 - qn4*qn4);
    y =  asin(-2*(qn2*qn4 - qn1*qn3));
    x = atan2( 2*(qn3*qn4 + qn1*qn2), qn1*qn1 - qn2*qn2 - qn3*qn3 + qn4*qn4);

    p_xyz[0]= x; 
    p_xyz[1]= y;
    p_xyz[2]= z;
}
