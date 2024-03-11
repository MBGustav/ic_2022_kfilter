#include <stdio.h>
#include <stdlib.h>
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
    
    sys.R = read_matrix("output/UnitTests/unit1_R.csv");
    sys.H = read_matrix("output/UnitTests/unit1_H.csv");
    struct matrix* exp_dx = read_matrix("output/UnitTests/unit1_dx.csv");
    struct matrix* exp_P = read_matrix("output/UnitTests/unit1_P.csv");
    struct matrix* dz = read_matrix("output/UnitTests/unit1_dz.csv");
    prop.P = read_matrix("output/UnitTests/unit1_Pp.csv");    
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

    
    struct matrix *exp_H1 = read_matrix("output/UnitTests/unit2_H1.csv");
    struct matrix *exp_R1 = read_matrix("output/UnitTests/unit2_R1.csv");
    struct matrix *ex_dz1 = read_matrix("output/UnitTests/unit2_dz1.csv");

    struct matrix *exp_H2 = read_matrix("output/UnitTests/unit2_H2.csv");
    struct matrix *exp_R2 = read_matrix("output/UnitTests/unit2_R2.csv");
    struct matrix *ex_dz2 = read_matrix("output/UnitTests/unit2_dz2.csv");
    
    struct matrix *exp_H3 = read_matrix("output/UnitTests/unit2_H3.csv");
    struct matrix *exp_R3 = read_matrix("output/UnitTests/unit2_R3.csv");
    struct matrix *ex_dz3 = read_matrix("output/UnitTests/unit2_dz3.csv");
    
    struct matrix *exp_H4 = read_matrix("output/UnitTests/unit2_H4.csv");
    struct matrix *exp_R4 = read_matrix("output/UnitTests/unit2_R4.csv");
    struct matrix *ex_dz4 = read_matrix("output/UnitTests/unit2_dz4.csv");
    
    struct matrix *exp_H5 = read_matrix("output/UnitTests/unit2_H5.csv");
    struct matrix *exp_R5 = read_matrix("output/UnitTests/unit2_R5.csv");
    struct matrix *ex_dz5 = read_matrix("output/UnitTests/unit2_dz5.csv");
    
    struct matrix *exp_H6 = read_matrix("output/UnitTests/unit2_H6.csv");
    struct matrix *exp_R6 = read_matrix("output/UnitTests/unit2_R6.csv");
    struct matrix *ex_dz6 = read_matrix("output/UnitTests/unit2_dz6.csv");
    
    struct matrix *exp_H7 = read_matrix("output/UnitTests/unit2_H7.csv");
    struct matrix *exp_R7 = read_matrix("output/UnitTests/unit2_R7.csv");
    struct matrix *ex_dz7 = read_matrix("output/UnitTests/unit2_dz7.csv");
    

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

#define check_res(cmd) if(!(cmd)) {printf(" Incorrect result from %s\n", #cmd); exit(EXIT_FAILURE);}

int main()
{
    bool result = true;

    printf("Testing\n");
    check_res(unit_kalman_filter());
    check_res(unit_measurement()); 
    


}