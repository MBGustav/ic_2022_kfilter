#ifndef _KALMAN_FILTER_HPP_
#define _KALMAN_FILTER_HPP_

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

#include "Quaternion.hpp"
#include "Linalg.hpp"



class KalmanFilter 
{
private:
    //Kalman Filter Parameters
    std_Matrix state_transition_matrix;         // F
    std_Matrix measurement_function;            // H
    std_Matrix measurement_noise;               // R
    std_Matrix process_noise;                   // Q
    std_Matrix priori_state_covariance;         // P
    std_Matrix transition_matrix;               // ? 
    std_Matrix kalman_gain;                     // K
    std_Matrix measurement_matrix;              // ? 
    std_Matrix control_transition_matrix;       // B  ? 


    //System Parameters
    std_Matrix lambda_gyr;
    std_Matrix lambda_acc;
    std_Matrix velocity_gps;
    std_Matrix position_gps;
    std_Matrix acceleration_bias;
    std_Matrix magnetometer; 
    std_Matrix gyroscope;
    std_Matrix delta_x; 
    std_Matrix delta_z; 
    
    //Constant Parameters
    std_Matrix constant_Rz;
    std_Matrix constant_Rx;

    //Space State
    std_Matrix state_bias_gyroscope;
    std_Matrix state_bias_acceleration;
    std_Matrix state_quaternion;
    std_Matrix state_position;
    std_Matrix state_velocity;
    std_Matrix state_acceleration;
    std_Matrix state_omegag;
    std_Matrix state_P; // state covariance

    //Propagated States
    std_Matrix propagated_bias_gyroscope;
    std_Matrix propagated_bias_acceleration;
    std_Matrix propagated_quaternion;
    std_Matrix propagated_position;
    std_Matrix propagated_velocity;
    std_Matrix propagated_acceleration;
    std_Matrix propagated_omegag;
    std_Matrix propagated_P; // covariance

    // System Measure constants
    std_Matrix acce;
    std_Matrix magn;
    std_Matrix pos_gps;
    std_Matrix vel_gps;
    // Variables to deal with mode size
    int mode;
    int dim_xm, dim_z;

    // Auxiliar Matrices (reduce OS time alloc/dealloc)
    std_Matrix me;
    std_Matrix ge;
    std_Matrix aux_prioriAndmeasure; /*step_kalman_gain*/
    std_Matrix HPHt;

    void meas_mode00();
    void meas_mode01();
    void meas_mode02();
    void meas_mode03();
    void meas_mode04();
    void meas_mode05();
    void meas_mode06();
    void meas_mode07();
    void meas_mode08();

public:
    void deallocate();
    KalmanFilter(int dim_x, int dim_z);

    void set_initial_state(Matrix &x0);
    void set_state_transition_matrix(Matrix &F); 
    void set_measurement_function(Matrix &H);
    void set_covariance(Matrix &P);
    void set_measurement_noise(Matrix &R);
    void set_process_noise(Matrix &Q);
    
    void predict();
    void update();/*preciso usar o delta_z ?*/
    void step_kalman_gain();
    
    void set_mode(int i);
    int get_mode();

    // void set_constant_parameters(
    //             std_Matrix constant_Rz,
    //             std_Matrix constant_Rx
    //             );
    
    ~KalmanFilter();

};

void deltaAngle2quat_f(Matrix &deltaAngle, Matrix &LGqpHat, Matrix &out_quaternion);


//#TODO: divisao de casos - elaborar testes com entradas fixas para os casos abaixo


// static inline struct matrix* 
//         meas_mode00(sys_param *sys, state_space *prop){
    
//     struct matrix *delta_z = matrix_m(3, 1);
//     struct matrix* H = matrix_m(1,15); 
//     struct matrix* R = matrix_m(3,15);
//     zeros(H);
//     zeros(R);

//     sys->H = H;
//     sys->R = R;
    
//     return NULL;
// }

// static inline struct matrix* 
//         meas_mode01(sys_param *sys, state_space *prop, struct matrix *me, struct matrix *magn, struct matrix *mag){
    
//     int Hr = 3, Hc = 15;
//     int Rr = 3, Rc = 3;
//     //saida matriz
//     struct matrix *delta_z = matrix_m(3, 1);
//     struct matrix *rot = matrix_m(3, 3);
//     struct matrix *subH1= matrix_m(3, 3);
//     struct matrix *skew_me= matrix_m(3, 3);
//     struct matrix *z_prop = matrix_m(3, 1);
//     // struct matrix *z_prop = matrix_m(3, 1);
//     struct matrix *me_tr = matrix_m(3,1);
//     struct matrix *mag_t = matrix_m(3,1);
//     struct matrix *ge_tr = matrix_m(3,1);

//     struct matrix *H = matrix_m(Hr, Hc);
//     struct matrix *R = matrix_m(Rr, Rc);
    
//     // rot = quat2dcm(prop.qtn');
//     quat2dcm(prop->qtn,rot);
//     //subH1 = rot * skew_me
//     crossM_f(me, skew_me);
//     times_m(rot,skew_me, subH1);
    
//     //Preenchimento H e R
//     zeros(H);
//     getpart_mat(subH1, 0, 0, H);

//     zeros(R);
//     mat_cpy(magn, R);

    
//     //z = mag; 
//     //z_prop = rot*me;
//     //delta_z = z - z_prop;
//     transp_m(me, me_tr);
//     times_m (rot, me_tr, z_prop);//z_prop = rot*me == (3,1)
//     transp_m(mag, mag_t);
//     less_m(mag_t, z_prop, delta_z);
//     // debug(z, 1000);
//     // debug(z_prop, 1000);

//     sys->H = H;
//     sys->R = R;
//     // debug(me,1);
//     // debug(z_prop,1);
//     // free aux memory
//     delete_m(rot);
//     delete_m(subH1);
//     delete_m(skew_me);
//     delete_m(z_prop);
//     delete_m(me_tr);
//     delete_m(mag_t);
//     delete_m(ge_tr);
//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode02(sys_param *sys, state_space *prop,struct matrix* acce ,struct matrix *ge)
// {
//     int Hr = 3, Hc = 15;
//     int Rr = 3, Rc = 3;

//     //matriz de saida - delta_z
//     struct matrix *delta_z = matrix_m(3, 1);
//     //Alloc aux matrices
//     struct matrix* rot    = matrix_m(3, 3);
//     struct matrix* ge_t   = matrix_m(3, 1); 
//     struct matrix* skew_ge= matrix_m(3, 3);
//     struct matrix* auxM1  = matrix_m(3, 3);/*rot * skew_ge = 3x3*/
//     struct matrix* I_3x3  = matrix_m(3, 3);eye_m(1, I_3x3);
//     struct matrix* z_prop = matrix_m(3, 1);
//     struct matrix* acc_tr = matrix_m(3, 1); // transposta de prop->acc
    
//     struct matrix* H = matrix_m(Hr,Hc);
//     struct matrix* R = matrix_m(Rr,Rc); 
//     zeros(H);
//     zeros(R);

//     quat2dcm(prop->qtn, rot);
//     crossM_f(ge, skew_ge);
//     times_m(rot, skew_ge, auxM1);
//     timesc_m(-1.0, auxM1, auxM1);


//     //preenchimento de H
//     getpart_mat(auxM1, 0, 0, H); //H(1:3,   1:3) = -rot * skew_ge
//     getpart_mat(I_3x3, 0, 9, H); //H(1:3, 10:12) = I
    

//     //preench. R
//     mat_cpy(acce, R);
//     transp_m(ge, ge_t);
//     times_m(rot, ge_t, z_prop);
    

//     //calculo de dz
//     transp_m(prop->acc, acc_tr);
//     // printf("time_m: ");
    
//     sum_m(acc_tr, z_prop, delta_z); //...dz = z - (-rot*ge) == z + rot*ge
//     // debug(z, 1000);
    
//     sys->R = R;
//     sys->H = H;

//     //Free memory
//     delete_m(rot); 
//     delete_m(ge_t); 
//     delete_m(skew_ge);
//     delete_m(auxM1); 
//     delete_m(I_3x3); 
//     delete_m(z_prop);
//     delete_m(acc_tr); 
//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode03(sys_param *sys, state_space *prop,struct matrix *magn,
//                             struct matrix *acce, struct matrix *mag, struct matrix *ge, struct matrix *me){
//     //matriz de saida - delta_z
//     int Hr = 6, Hc = 15;
//     int Rr = 6, Rc = 6;

//     struct matrix *delta_z = matrix_m(6, 1);
//     struct matrix *rot    = matrix_m(3, 3);
//     struct matrix *skew_me= matrix_m(3, 3);
//     struct matrix *skew_ge= matrix_m(3, 3);
//     struct matrix *skew_me_aux= matrix_m(3, 3);
//     struct matrix *skew_ge_aux= matrix_m(3, 3);
//     struct matrix *I_3x3  = matrix_m(3, 3);eye_m(1, I_3x3);
//     struct matrix *z_prop = matrix_m(6, 1); 
//     struct matrix *z = matrix_m(6, 1); 
//     struct matrix *auxV1  = matrix_m(3, 1); /*==  rot * me*/
//     struct matrix *auxV2  = matrix_m(3, 1); /*== -rot * ge*/
//     struct matrix* mag_t = matrix_m(3, 1);
//     struct matrix* acc_t = matrix_m(3, 1);
//     struct matrix* ge_t  = matrix_m(3, 1);
//     struct matrix* me_t  = matrix_m(3, 1);
    
//     struct matrix* H = matrix_m(Hr, Hc);
//     struct matrix* R = matrix_m(Rr, Rc);
    
//     // skew_me = -rot x crossM_f(me) 
//     quat2dcm(prop->qtn, rot);
//     crossM_f(me,skew_me_aux); 
//     times_m(rot, skew_me_aux, skew_me);
    
//     // skew_ge = -rot x crossM_f(ge)
//     crossM_f(ge,skew_ge_aux); 
//     times_m(rot, skew_ge_aux, skew_ge);
//     timesc_m(-1, skew_ge, skew_ge);
        
//     //Inclui matrizes em H
//     zeros(H);
//     getpart_mat(skew_me, 0, 0, H);
//     getpart_mat(skew_ge, 3, 0, H);
//     getpart_mat(I_3x3 ,  3, 9, H);

//     //Escrita na matriz R
//     zeros(R);
//     getpart_mat(magn, 0, 0, R);
//     getpart_mat(acce, 3, 3, R);//escrito na diagonal

//     // Calculo Z
//     transp_m(mag, mag_t);
//     transp_m(prop->acc, acc_t);
    
//     getpart_mat(mag_t, 0, 0, z);
//     getpart_mat(acc_t, 3, 0, z);//copia a coluna em z = [mag ; prop.acc]

//     // rot*me
//     transp_m(me, me_t);
//     times_m(rot, me_t, auxV1);

//     // -rot*ge
//     transp_m(ge, ge_t);
//     times_m(rot, ge_t , auxV2);
//     timesc_m(-1, auxV2, auxV2);

//     getpart_mat(auxV1, 0, 0, z_prop);
//     getpart_mat(auxV2, 3, 0, z_prop);
//     less_m(z, z_prop, delta_z);
    
//     sys->H = H; 
//     sys->R = R;
    
//     //Free memory
//     delete_m( rot ); 
//     delete_m( skew_me ); 
//     delete_m( skew_ge ); 
//     delete_m( I_3x3 ); 
//     delete_m( z_prop ); 
//     delete_m( auxV1 ); 
//     delete_m( auxV2 ); 
//     delete_m(mag_t);
//     delete_m(acc_t);
//     delete_m(ge_t);
//     delete_m(me_t);
    
//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode04(sys_param *sys, state_space *prop,
//     struct matrix* vel_gps, struct matrix* pos_gps, struct matrix* gps_position){

//     int Hr = 6, Hc = 15;
//     int Rr = 6, Rc = 6;

//     //Alloc aux matrices
//     struct matrix *I_3x3  = matrix_m(3, 3);eye_m(1.0, I_3x3);
//     struct matrix *z      = matrix_m(6, 1);
//     struct matrix *z_prop = matrix_m(6, 1);
//     struct matrix *gps_position_t= matrix_m(3,1);
//     struct matrix *prop_vel_t = matrix_m(3,1);
//     struct matrix *prop_pos_t= matrix_m(3,1);
//     //matriz de saida - delta_z
//     struct matrix *delta_z = matrix_m(6, 1);
//     struct matrix *H = matrix_m(Hr, Hc);
//     struct matrix *R = matrix_m(Rr, Rc);

//     //montagem H
//     zeros(H);
//     getpart_mat(I_3x3, 0,  6, H); /*H(1:3,7:9)   = eye(3);*/
//     getpart_mat(I_3x3, 3, 12, H); /*H(4:6,13:15) = eye(3);*/

//     //montagem R
//     zeros(R);
//     getpart_mat(vel_gps, 0, 0, R);
//     getpart_mat(pos_gps, 3, 3, R);

//     //Montagem z
//     zeros(z);
//     transp_m(gps_position, gps_position_t);
//     getpart_mat(gps_position_t, 3, 0, z);
    
//     //Montagem z_prop
//     transp_m(prop->vel, prop_vel_t);
//     transp_m(prop->pos, prop_pos_t);

//     getpart_mat(prop_vel_t , 0, 0, z_prop);
//     getpart_mat(prop_pos_t , 3, 0, z_prop);
    
//     // debug(prop_pos_t,1);
//     // debug(prop_vel_t,1);
//     // debug(z_prop, 1);
//     // debug(z,1);
    
//     less_m(z, z_prop, delta_z);
//     // debug(delta_z,1);
//     sys->H = H;
//     sys->R = R;


//     //Free memory
//     delete_m(I_3x3);
//     delete_m(z);
//     delete_m(z_prop);
//     delete_m(gps_position_t);
//     delete_m(prop_vel_t);
//     delete_m(prop_pos_t);

//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode05(sys_param *sys, state_space *prop, 
//                             struct matrix* magn, struct matrix* vel_gps, struct matrix* pos_gps,
//                             struct matrix* mag, struct matrix* gps_pos,struct matrix* ge, struct matrix *me){
//     //matriz de saida - delta_z
//     struct matrix *delta_z = matrix_m(9, 1);
//     //Alloc aux matrices
//     struct matrix *rot    = matrix_m(3, 3);
//     struct matrix *skew_me= matrix_m(3, 3);
//     struct matrix *subH1 = matrix_m(3, 3);
//     struct matrix *I_3x3  = matrix_m(3, 3);
//     struct matrix *z      = matrix_m(9, 1);
//     struct matrix *z_prop = matrix_m(9, 1);
//     struct matrix *auxV1  = matrix_m(3, 1);//auxV1 = rot(3x3) x me(3x1);
//     struct matrix *mag_t  = matrix_m(3, 1);
//     struct matrix *gps_pos_t= matrix_m(3, 1);
//     struct matrix *prop_vel_t=matrix_m(3,1);
//     struct matrix *prop_pos_t=matrix_m(3,1);
//     struct matrix *me_t=matrix_m(3,1);

//     struct matrix* H = matrix_m(9, 15);
//     struct matrix* R = matrix_m(9,  9);
//     eye_m(1, I_3x3);

//     quat2dcm(prop->qtn, rot);
//     crossM_f(me, skew_me); 

//     //montagem matriz H
//     zeros(H);
//     times_m(rot, skew_me, subH1);
//     getpart_mat(subH1,  0, 0, H); //H(1:3,1:3) = rot*skew_me;
//     getpart_mat(I_3x3,  3, 6, H); //H(4:6,7:9) = eye(3);
//     getpart_mat(I_3x3,  6,12, H); //H(7:9,13:15) = eye(3);
    
//     //montagem matriz R
//     zeros(R);
//     getpart_mat(magn   , 0, 0, R);
//     getpart_mat(vel_gps, 3, 3, R); 
//     getpart_mat(pos_gps, 6, 6, R); 

//     //montagem de z
//     zeros(z);
//     transp_m(mag, mag_t);
//     transp_m(gps_pos, gps_pos_t);
//     getpart_mat(mag_t    , 0, 0, z);
//     getpart_mat(gps_pos_t, 6, 0, z);
    
//     //montagem z_prop
//     transp_m(me, me_t);
//     transp_m(prop->vel, prop_vel_t); //prop.vel
//     transp_m(prop->pos, prop_pos_t); //prop.pos
//     times_m(rot, me_t, auxV1);
//     getpart_mat(auxV1     , 0, 0, z_prop);
//     getpart_mat(prop_vel_t, 3, 0, z_prop);
//     getpart_mat(prop_pos_t, 6, 0, z_prop);
//     less_m(z, z_prop, delta_z);
//     sys->H = H;
//     sys->R = R;
//     //Free memory
//     delete_m(rot);
//     delete_m(subH1);
//     delete_m(skew_me);
//     delete_m(I_3x3);
//     delete_m(z);
//     delete_m(z_prop);
//     delete_m(auxV1);
//     // exit(0);
//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode06(sys_param *sys, state_space *prop, struct matrix *acce, struct matrix* gps_pos,
//                                struct matrix *vel_gps,  struct matrix *pos_gps, struct matrix *ge){
//     //matriz de saida - delta_z
//     struct matrix *delta_z = matrix_m(9, 1);
//     //Alloc aux matrices
//     struct matrix *rot      = matrix_m(3, 3);
//     struct matrix *rot_ge   = matrix_m(3, 1);
//     struct matrix *skew_ge  = matrix_m(3, 3);
//     struct matrix *I_3x3    = matrix_m(3, 3);eye_m(1.0f, I_3x3);
//     struct matrix *z        = matrix_m(9, 1);
//     struct matrix *z_prop   = matrix_m(9, 1);
//     struct matrix *auxM1    = matrix_m(3, 3);/*-rot x skew_ge */
//     struct matrix *ge_t     = matrix_m(3, 1);
//     struct matrix *acc_t    = matrix_m(3, 1);
//     struct matrix *gps_pos_t= matrix_m(3, 1);
//     struct matrix *prop_vt  = matrix_m(3, 1);
//     struct matrix *prop_pt  = matrix_m(3, 1);
//     struct matrix *H = matrix_m(9, 15);
//     struct matrix *R = matrix_m(9 ,9);


//     quat2dcm(prop->qtn, rot);
//     crossM_f(ge, skew_ge);
//     transp_m(ge, ge_t);
//     times_m(rot, skew_ge, auxM1);
//     timesc_m(-1, auxM1, auxM1);

//     // Montagem H
//     getpart_mat(auxM1, 0,  0, H); /*H(1:3,1:3) = -rot*skew_ge;*/
//     getpart_mat(I_3x3, 0,  9, H); /*H(1:3,10:12) = eye(3);*/
//     getpart_mat(I_3x3, 3,  6, H); /*H(4:6,  7:9) = eye(3);*/
//     getpart_mat(I_3x3, 6, 12, H); /*H(7:9,13:15) = eye(3);*/

//     //Montagem R
//     getpart_mat(acce   , 0, 0, R); // montagem de R
//     getpart_mat(vel_gps, 3, 3, R);
//     getpart_mat(pos_gps, 6, 6, R);

//     //Montagem Z
//     transp_m(prop->acc, acc_t);  
//     transp_m(gps_pos, gps_pos_t);  
//     getpart_mat(acc_t     , 0, 0, z); // montagem de z
//     getpart_mat(gps_pos_t , 6, 0, z);
    
//     //Montagem z_prop
//     times_m(rot, ge_t, rot_ge);
//     timesc_m(-1.0, rot_ge, rot_ge);
//     transp_m(prop->vel, prop_vt);
//     transp_m(prop->pos, prop_pt);
//     getpart_mat(rot_ge , 0, 0, z_prop);
//     getpart_mat(prop_vt, 3, 0, z_prop);
//     getpart_mat(prop_pt, 6, 0, z_prop);
//     less_m(z, z_prop, delta_z);
    
//     sys->H = H;
//     sys->R = R;

//     // debug(sys->H, 1000);
//     // debug(H, 1000);

//     // debug(sys->R, 1000);
//     // debug(R, 1000);
    
//     //Free memory
//     delete_m(rot_ge);
//     delete_m(rot);
//     delete_m(skew_ge);
//     delete_m(I_3x3);
//     delete_m(z);
//     delete_m(z_prop);
//     delete_m(auxM1);
//     delete_m(ge_t);
//     delete_m(acc_t);
//     delete_m(gps_pos_t);
//     delete_m(prop_vt);
//     delete_m(prop_pt);


//     return delta_z;
// }

// static inline struct matrix* 
//         meas_mode07(sys_param *sys, state_space *prop, struct matrix* gps_pos, struct matrix *mag,
//                     struct matrix* magn, struct matrix* acce, struct matrix* vel_gps, struct matrix* pos_gps, struct matrix* ge, struct matrix* me ){
//     //matriz de saida - delta_z
//     struct matrix *delta_z = matrix_m(12, 1);
//     struct matrix *delta_z_t= matrix_m(1, 12);
    
//     //Alloc aux matrices
//     struct matrix*rot = matrix_m(3, 3); 
//     struct matrix*skew_me = matrix_m(3, 3); 
//     struct matrix*skew_ge = matrix_m(3, 3); 
//     struct matrix*auxM1  = matrix_m(3, 3); //  rot * skew_me
//     struct matrix*auxM2  = matrix_m(3, 3); // -rot * skew_ge
//     struct matrix*auxM3  = matrix_m(3, 1); //  rot * me
//     struct matrix*auxM4  = matrix_m(3, 1); // -rot * ge
//     struct matrix*auxM3_t= matrix_m(1, 3); // ( rot * me)t
//     struct matrix*auxM4_t= matrix_m(1, 3); // (-rot * ge)t
//     struct matrix*I_3x3  = matrix_m(3, 3); eye_m(1.0, I_3x3);
//     struct matrix*z_prop = matrix_m(1, 12); 
//     struct matrix*z      = matrix_m(1, 12); 
//     struct matrix*me_t  = matrix_m(3, 1); 
//     struct matrix*ge_t  = matrix_m(3, 1); 
//     struct matrix*prop_acc_t= matrix_m(3, 1); 
    
//     struct matrix* H = matrix_m(12, 15); 
//     struct matrix* R = matrix_m(12, 12);

//     quat2dcm(prop->qtn, rot);
//     crossM_f(me, skew_me);
//     crossM_f(ge, skew_ge);

//     //rot * skew_me
//     times_m(rot, skew_me, auxM1);

//     // -rot * skew_ge
//     times_m(rot, skew_ge, auxM2);
//     timesc_m(-1.0, auxM2, auxM2);

//     //rot * me
//     transp_m(me, me_t);
//     times_m(rot, me_t, auxM3);
    
//     // - rot*ge
//     transp_m(ge, ge_t);
//     times_m(rot, ge_t, auxM4);
//     timesc_m(-1.0, auxM4, auxM4);

//     //montagem H
//     zeros(H);
//     getpart_mat(auxM1 , 0, 0, H); //H(1:3,1:3) = rot*skew_me;
//     getpart_mat(auxM2 , 3, 0, H); //H(4:6,1:3) = -rot*skew_ge;
//     getpart_mat(I_3x3 , 3, 9, H); //H(4:6,10:12) = eye(3);
//     getpart_mat(I_3x3 , 6, 6, H); //H(7:9,7:9) = eye(3);
//     getpart_mat(I_3x3 , 9, 12,H); //H(10:12,13:15) = eye(3);

//     // montagem R
//     zeros(R);
//     getpart_mat(magn   , 0, 0, R); 
//     getpart_mat(acce   , 3, 3, R);
//     getpart_mat(vel_gps, 6, 6, R);
//     getpart_mat(pos_gps, 9, 9, R);
    
//     //montagem z 
//     zeros(z);
//     getpart_mat(mag      , 0, 0, z); 
//     getpart_mat(prop->acc, 0, 3, z);
//     getpart_mat(gps_pos  , 0, 9, z);

//     //montagem z_prop
//     transp_m(auxM3, auxM3_t);
//     transp_m(auxM4, auxM4_t);
//     getpart_mat(auxM3_t   , 0, 0, z_prop); // rot*me
//     getpart_mat(auxM4_t   , 0, 3, z_prop); //-rot*ge
//     getpart_mat(prop->vel , 0, 6, z_prop);
//     getpart_mat(prop->pos , 0, 9, z_prop);

//     less_m(z,z_prop, delta_z_t);
//     transp_m(delta_z_t, delta_z);

//     sys->H = H;
//     sys->R = R;


//     //Free memory
//     delete_m( rot );
//     delete_m( skew_me );
//     delete_m( skew_ge );
//     delete_m( auxM1 );
//     delete_m( auxM2 );
//     delete_m( auxM3 );
//     delete_m( auxM4 );
//     delete_m( I_3x3 );
//     delete_m( z_prop );
//     delete_m( z );
//     delete_m(me_t);
//     delete_m(ge_t);
//     delete_m( prop_acc_t);
//     delete_m( auxM3_t);
//     delete_m( auxM4_t);
//     delete_m( delta_z_t);


//     return delta_z;
// }

// void KalmanFilter::measure()
// {
//     std_Matrix acce(3,3);
//     std_Matrix magn(3,3);
//     std_Matrix pos_gps(3,3);
//     std_Matrix vel_gps(3,3);



//     switch(mode){
//         case 0:{
//             this->meas_mode00();
//             break;
//         } 
//         case 1:{
//             this->meas_mode01();
//             break;
//         }
//         case 2:{
//             this->meas_mode02();
//             break;
//         }
//         case 3:{
//             this->meas_mode03();
//             break;
//         }
//         case 4:{
//             this->meas_mode04();
//             break;
//         }
//         case 5:{
//             this->meas_mode05();
//             break;
//         }
//         case 6:{
//             this->meas_mode06();
//             break;
//         }
//         case 7:{
//             this->meas_mode07();
//             break;
//         }
//         default:

//     }


// }

// //Calculo do vetor(coluna) de medidas, sendo seu tamanho variavel. Faz chamada de acordo com o valor "mode"
// struct matrix*  measurement(sys_param *sys, int mode,state_space *prop,struct matrix* mag,
//                  struct matrix* gps_position,struct matrix* ge,
//                  struct matrix* me){ 
    
    
//     //  Variância da aceleração linear via acelerômetro (calculada de 0 a 10 segundos)
//     struct matrix* acce = matrix_m(3,3);
//     // % Variância da intensidade do campo magnético via magnetômetro
//     struct matrix* magn = matrix_m(3,3);
//     // % Variância da posição linear XY via GPS
//     struct matrix* pos_gps = matrix_m(3,3);
//     // % Variância da velocidade linear XY via GPS
//     struct matrix* vel_gps = matrix_m(3,3);
    
    
//     struct matrix* delta_z=NULL;//#TODO : como ldar com H (ou delta_z abstrato como 0 ?)
//     // #todo: remover esses zeros ?
//     zeros(acce);
//     zeros(magn);
//     zeros(pos_gps);
//     zeros(vel_gps);
//     // Valores inseridos na diagonal
//     acce->elements[0] = acce_x;acce->elements[4] = acce_y;acce->elements[8] = acce_z;
//     magn->elements[0] = magn_x;magn->elements[4] = magn_y;magn->elements[8] = magn_z;
//     pos_gps->elements[0] = pos_gps_x;pos_gps->elements[4] = pos_gps_y;pos_gps->elements[8] = pos_gps_z;
//     vel_gps->elements[0] = vel_gps_x;vel_gps->elements[4] = vel_gps_y;vel_gps->elements[8] = vel_gps_z;
//     //todos os casos modificam o tamanho de matriz

//     delete_m(sys->H);
//     delete_m(sys->R);


//     switch (mode){
//     case 0:{ //% ---/---/---/---   
//         delta_z = meas_mode00(sys, prop);// OK
//         break;
//     }
//     case 1:{//% ---/---/---/mag
//         delta_z = meas_mode01(sys, prop, me, magn, mag);
//         break;
//     }
//     case 2:{//% ---/---/acc/---
//         delta_z = meas_mode02(sys, prop, acce, ge); //OK
//         break;
//     }
//     case 3:{//% ---/---/acc/mag
//         delta_z = meas_mode03(sys, prop, magn, acce, mag, ge, me); // OK
//         break;
//     }
//     case 4:{ // ---/gps/---/---
//         delta_z = meas_mode04(sys,prop, vel_gps, pos_gps, gps_position); // #CHECK
//         break; 
//     }
//     case 5:{//% ---/gps/---/mag   
//         delta_z = meas_mode05(sys, prop, magn, vel_gps, pos_gps, mag, gps_position, ge, me);
//         break; 
//     }
//     case 6:{ //% ---/gps/acc/---
//         delta_z = meas_mode06(sys,prop, acce, gps_position, vel_gps, pos_gps, ge); 
//         break; 
//     }
//     case 7:{ //% ---/gps/acc/mag
//         delta_z = meas_mode07(sys, prop, gps_position, mag, magn, acce, vel_gps, pos_gps, ge,me);
//         break; 
//     }
//     default: bound_error();
//         break;
//     }

//     //considerando que nao houve falha de borda (bound)
//     // cont_mode->elements[mode]+=1;
    
//     // debug(sys->H, 1000);
//     // debug(sys->R, 1000);
//     // debug(delta_z,1000);


//     delete_m(acce);
//     delete_m(magn);
//     delete_m(pos_gps);
//     delete_m(vel_gps);


//     return delta_z;
// }




#endif /*_KALMAN_FILTER_HPP_*/
