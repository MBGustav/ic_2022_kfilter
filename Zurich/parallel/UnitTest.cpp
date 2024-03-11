//  MKL kernels
#include <sycl/sycl.hpp>
#include "oneapi/mkl.hpp"
#include "kernels/mkl_kernels.hpp"


#include <iostream>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
//usr libs
#include "Timer.h"
#include "structs.h"
#include "defines.h"
#include "libmat.h"
#include "kalman_filter.h"


bool unit_kalman_filter()
{
    bool test = true;
    int mode = 1;
    int nx   = 8;
    state_space prop;
    sys_param sys;
    
    sys.R = read_matrix((char*) "output/UnitTests/unit1_R.csv");
    sys.H = read_matrix((char*) "output/UnitTests/unit1_H.csv");
    struct matrix* exp_dx = read_matrix((char*) "output/UnitTests/unit1_dx.csv");
    struct matrix* exp_P = read_matrix((char*) "output/UnitTests/unit1_P.csv");
    struct matrix* dz = read_matrix((char*) "output/UnitTests/unit1_dz.csv");
    prop.P = read_matrix((char*) "output/UnitTests/unit1_Pp.csv");    
    struct matrix *P = matrix_m(nx, nx);
    struct matrix *dx = kalman_filter(mode, &prop, dz, &sys, nx,P);
    
    test = test && matrix_isEqual(P, exp_P);
    test = test && matrix_isEqual(dx, exp_dx); 

    delete_m(P);
    delete_m(dx);
    delete_m(dz);
    delete_m(exp_P);
    delete_m(exp_dx);
    delete_m(sys.H);
    delete_m(sys.R);

    return test;
}
//Damn, thats too much cases... Dx
bool unit_measurement()
{
    bool res = true;
    int mode = 1;

    
    struct matrix *exp_H1 = read_matrix((char*) "output/UnitTests/unit2_H1.csv");
    struct matrix *exp_R1 = read_matrix((char*) "output/UnitTests/unit2_R1.csv");
    struct matrix *ex_dz1 = read_matrix((char*) "output/UnitTests/unit2_dz1.csv");

    struct matrix *exp_H2 = read_matrix((char*) "output/UnitTests/unit2_H2.csv");
    struct matrix *exp_R2 = read_matrix((char*) "output/UnitTests/unit2_R2.csv");
    struct matrix *ex_dz2 = read_matrix((char*) "output/UnitTests/unit2_dz2.csv");
    
    struct matrix *exp_H3 = read_matrix((char*) "output/UnitTests/unit2_H3.csv");
    struct matrix *exp_R3 = read_matrix((char*) "output/UnitTests/unit2_R3.csv");
    struct matrix *ex_dz3 = read_matrix((char*) "output/UnitTests/unit2_dz3.csv");
    
    struct matrix *exp_H4 = read_matrix((char*) "output/UnitTests/unit2_H4.csv");
    struct matrix *exp_R4 = read_matrix((char*) "output/UnitTests/unit2_R4.csv");
    struct matrix *ex_dz4 = read_matrix((char*) "output/UnitTests/unit2_dz4.csv");
    
    struct matrix *exp_H5 = read_matrix((char*) "output/UnitTests/unit2_H5.csv");
    struct matrix *exp_R5 = read_matrix((char*) "output/UnitTests/unit2_R5.csv");
    struct matrix *ex_dz5 = read_matrix((char*) "output/UnitTests/unit2_dz5.csv");
    
    struct matrix *exp_H6 = read_matrix((char*) "output/UnitTests/unit2_H6.csv");
    struct matrix *exp_R6 = read_matrix((char*) "output/UnitTests/unit2_R6.csv");
    struct matrix *ex_dz6 = read_matrix((char*) "output/UnitTests/unit2_dz6.csv");
    
    struct matrix *exp_H7 = read_matrix((char*) "output/UnitTests/unit2_H7.csv");
    struct matrix *exp_R7 = read_matrix((char*) "output/UnitTests/unit2_R7.csv");
    struct matrix *ex_dz7 = read_matrix((char*) "output/UnitTests/unit2_dz7.csv");
    

    state_space prop; 
    sys_param sys;
    sys.H = NULL;
    sys.R = NULL;

    struct matrix*ge  = matrix_m(1,3);
    struct matrix*me  = matrix_m(1,3);
    struct matrix*mag = matrix_m(1,3);
    struct matrix*gps = matrix_m(1,3);

    prop.qtn = matrix_m(1,4);
    prop.acc = matrix_m(1,3);
    prop.vel = matrix_m(1,3);
    prop.pos = matrix_m(1,3);


    ge->elements[0] = 0.0;
    ge->elements[1] = 0.0;
    ge->elements[2] = 9.81;

    mag->elements[0] = 5.7;
    mag->elements[1] = 6.5;
    mag->elements[2] = 4.5;

    me->elements[0] = norm_m(mag);
    me->elements[1] = 0.0;
    me->elements[2] = 0.0;

    mag->elements[0] = 5.7;
    mag->elements[1] = 6.5;
    mag->elements[2] = 4.5;

    prop.qtn->elements[0] =0.8147;
    prop.qtn->elements[1] =0.9058;
    prop.qtn->elements[2] =0.1270;
    prop.qtn->elements[3] =0.9134;

    prop.acc->elements[0] =0.5469;
    prop.acc->elements[1] =0.9575;
    prop.acc->elements[2] =0.9649;
    
    prop.vel->elements[0]= 0.1576;        
    prop.vel->elements[1]= 0.9706;
    prop.vel->elements[2]= 0.9572;

    prop.pos->elements[0]= 0.4854;
    prop.pos->elements[1]= 0.8003;
    prop.pos->elements[2]= 0.1419;

    gps->elements[0] =0.6324;
    gps->elements[1] =0.0975;
    gps->elements[2] =0.2785;
    
    struct matrix* dz1 = measurement(&sys, 1, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz1, ex_dz1);
    res = res && matrix_isEqual(sys.H, exp_H1);
    res = res && matrix_isEqual(sys.R, exp_R1);

    struct matrix* dz2 = measurement(&sys, 2, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz2, ex_dz2);
    res = res && matrix_isEqual(sys.H, exp_H2);
    res = res && matrix_isEqual(sys.R, exp_R2);

    struct matrix* dz3 = measurement(&sys, 3, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz3, ex_dz3);
    res = res && matrix_isEqual(sys.H, exp_H3);
    res = res && matrix_isEqual(sys.R, exp_R3);
    
    struct matrix* dz4 = measurement(&sys, 4, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(sys.H, exp_H4);
    res = res && matrix_isEqual(sys.R, exp_R4);
    res = res && matrix_isEqual(dz4, ex_dz4);

    struct matrix* dz5 = measurement(&sys, 5, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz5, ex_dz5);
    res = res && matrix_isEqual(sys.H, exp_H5);
    res = res && matrix_isEqual(sys.R, exp_R5);

    struct matrix* dz6 = measurement(&sys, 6, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz6, ex_dz6);
    res = res && matrix_isEqual(sys.H, exp_H6);
    res = res && matrix_isEqual(sys.R, exp_R6);
    
    struct matrix* dz7 = measurement(&sys, 7, &prop, mag, gps, ge, me);
    res = res && matrix_isEqual(dz7, ex_dz7);
    res = res && matrix_isEqual(sys.H, exp_H7);
    res = res && matrix_isEqual(sys.R, exp_R7);
    
    return res;
}

struct matrix* cov_matrix_serial(struct matrix* F, struct matrix* P, struct matrix* G, struct matrix* Q)
{
    int Gr = G->n_row;
    int Gc = G->n_col;
    int Qr = Q->n_row;
    int Qc = Q->n_col;
    int Fr = F->n_row;
    int Fc = F->n_col;
    int Pr = P->n_row;
    int Pc = P->n_col;
    
    struct matrix* Gt = matrix_m(Gc, Gr);
    struct matrix* Ft = matrix_m(Fc, Fr);
    struct matrix* FP = matrix_m(Fr, Pc);
    struct matrix* GQ = matrix_m(Gr, Qc);
    struct matrix* FPF = matrix_m(Fr, Pc);
    struct matrix* GQG = matrix_m(Fr, Pc);

    struct matrix* propP = matrix_m(Fr, Pc);
    
    // FPF  = F * P * F'
    times_m(F, P, FP);      
    transp_m(F, Ft);
    times_m(FP, Ft, FPF);
    
    // GQG = G x Q x G'
    times_m(G, Q, GQ);
    transp_m(G, Gt);
    times_m(GQ, Gt, GQG);

    sum_m(FPF, GQG, propP); 

    delete_m (Gt);
    delete_m (Ft);
    delete_m (FP);
    delete_m (GQ);
    delete_m (FPF);
    delete_m (GQG);

    return propP;
}
void cov_matrix_test(int m, Timer *T)
{
    sycl::device dev = device_runner();
    
    sycl::queue q;

    struct matrix *F = rnd_matrix_usm(q, m, m);
    struct matrix *P = rnd_matrix_usm(q, m, m);
    struct matrix *G = rnd_matrix_usm(q, m, m);
    struct matrix *Q = rnd_matrix_usm(q, m, m);

    struct matrix *prpP = matrix_m(m, m);
    

    timer_start(T);
    mkl_cov_pred(q, F, P, G, Q, prpP);
    timer_stop(T, "par");

    timer_start(T);
    struct matrix *expP = cov_matrix_serial(F, P, G, Q);
    timer_stop(T, "ser");

    
    delete_m_usm(F,q);
    delete_m_usm(P,q);
    delete_m_usm(G,q);
    delete_m_usm(Q,q);
    delete_m(prpP);
    delete_m(expP);
}


bool unit_cov_matrix_offload()
{

    struct matrix *F = read_matrix((char*) "output/UnitTests/Unit3_F.csv");
    struct matrix *P = read_matrix((char*) "output/UnitTests/Unit3_P.csv");
    struct matrix *G = read_matrix((char*) "output/UnitTests/Unit3_G.csv");
    struct matrix *Q = read_matrix((char*) "output/UnitTests/Unit3_Q.csv");
    // Declare generalized square matrix 
    int m = F->n_row;
    struct matrix *propP = matrix_m(m, m);
    struct matrix *expP = cov_matrix_serial(F, P, G, Q);
    // print_m(expP,1);

    sycl::device dev = device_runner();
    sycl::queue q;
    // data offload and select device
    mkl_cov_pred(q, F, P, G, Q, propP);


    //check results
    return matrix_isEqual(expP, propP );
}


bool unit_gemm()
{
    int m = 8; 
    int n = 3;

    struct matrix*A = rnd_matrix(m, n);
    struct matrix*At  = matrix_m(n, m);
    struct matrix*C   = matrix_m(m, n);
    struct matrix*Ct  = matrix_m(n, m);

    struct matrix* I_m =matrix_m(m, m);
    struct matrix* I_n =matrix_m(n, n);

    eye_m(1.0f, I_m); 
    eye_m(1.0f, I_n); 
    transp_m(A, At);


    sycl::device dev = device_runner();
    //without transpose
    oneapi_matmul(dev, A, I_n, C, false, false, 1.0, 0.0);

    //with transpose
    oneapi_matmul(dev, At, I_m, Ct, false, false, 1.0, 0.0);

    delete_m(A);
    delete_m(At);
    delete_m(C);
    delete_m(Ct);
    delete_m(I_m);
    delete_m(I_n);

    return matrix_isEqual(A,C) && matrix_isEqual(At,Ct);
}

#define check_res(cmd) if(!(cmd)) {printf(" Incorrect result from %s\n", #cmd); exit(EXIT_FAILURE);}

int main()
{
    bool result = true;

    printf("Testing\n");
    // check_res(unit_kalman_filter());
    // check_res(unit_measurement()); 
    // check_res(unit_gemm());
    // check_res(unit_cov_matrix_offload());

    sycl::queue q{sycl::gpu_selector_v};
    Timer T;
    timer_setconf("compare_cvm.csv", &T, "name,size");
    char cov_name[] = "cov_matrix,200";
    // for(int i = 0; i < 30; i+=15)    
    cov_matrix_test(200, &T);



    //if everything is okay,
    printf("All tests passed\n");

    return 0;


}