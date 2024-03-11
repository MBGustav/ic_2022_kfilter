
#include "Euler.h"


struct matrix *NewEuler(int tam){
    struct matrix *Euler = matrix_m(tam, 3);
    return Euler;
}

struct matrix* static_attitude(IMU_data *imu_data, MAG_data *mag_data){
    //#TODO: Como checar esta função ? 
    
    //Declaring initial variables
    int tam_angle, tam;
    int j = 1; 
    double num, den;
    double roll, pitch, yaw;
    tam_angle = mag_data->tim->n_row;
    tam = imu_data->tim->n_row;

    struct matrix* acc = matrix_m(1, 3);
    struct matrix* euler_m = matrix_m(tam_angle, 3); // roll, pitch, yaw
    struct matrix* acc_imu = matrix_m(1, 3); //imu_data.acc(i);
    if(!euler_m) alloc_error();

    double *pAcc = acc->elements;
    double *pMag = mag_data->mag->elements;
    double *pMag_tim = mag_data->tim->elements;
    double *pImu_tim = imu_data->tim->elements; 
    
    for(int i = 0; i < tam; i++){
        if(pMag_tim[j] < pImu_tim[i] || FLT_EQ(pMag_tim[j], pImu_tim[i])){ //fix: usar "<" em vez de "=<"
            //calculo acc
            mat_cpyRow(imu_data->acc, i, acc_imu, 0);
            ned2enu(acc_imu, acc);//OK
            
            //calculo -> roll
            roll  = atan2d(pAcc[1], pAcc[2])*(180/PI); /*roll = euler[j,1]*/
            
            //calculo -> pitch
            pitch = atan2d(-pAcc[0], sqrt(pAcc[1]*pAcc[1] + pAcc[2]*pAcc[2]))*(180/PI);
            
            //calculo -> yaw
            num = pMag[3*(j)+2] * sind(roll)  - pMag[3*j+1] * cosd(roll);
            den = pMag[3*(j)] * cosd(pitch) + pMag[3*j+1] * sind(pitch) * sind(roll) + pMag[3*j+2] * sind(pitch) * cosd(roll);
            yaw = atan2d(num, den) * (180/PI) - 2*constant_beta;
            //printf("nessted loop:(%d, %d) %.4lf,  %.4lf,  %.4lf\n",i+1, j , roll, pitch, yaw);
            //if(j== 20) exit(0);
            euler_m->elements[(j)*3+0] = roll;
            euler_m->elements[(j)*3+1] = pitch;
            euler_m->elements[(j)*3+2] = yaw;
            if(j != tam_angle) j++;
            //debug(acc);
        }


    }

    
    
    delete_m( acc);
    delete_m( acc_imu);
    return euler_m;
}


void ned2enu(struct matrix *v, struct matrix *w){

    // Converte vetores da coordenada NED para ENU.
    w->elements[0] =  v->elements[1];
    w->elements[1] =  v->elements[0];
    w->elements[2] = -v->elements[2];
}