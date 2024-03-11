
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "Quaternions.h"
#include "libmat.h"
#include "structs.h"
#include "Euler.h"
#include "Initialization.h"
#include "err_handler.h"

bool LessThanOrEq(__data_type x, __data_type y){
    const __data_type epsilon = 0.000000000000001f;
    return ((x-y) < epsilon || fabs(x-y) < epsilon);
}


bool check_available(struct matrix *imu_time, struct matrix *data_tim,
                         int k, int *last_val, int max_idx)
{
    //check memory
    if(!imu_time || !data_tim) alloc_error();
    
    //parameters
    __data_type *dt_tim = data_tim->elements;
    __data_type *imu_tim = imu_time->elements;
    int sz_imu = imu_time->n_row;
    int sz_par = data_tim->n_row;
    bool eq = (LessThanOrEq( dt_tim[*last_val+1], imu_tim[k]));

    if(eq && (*last_val+1 != max_idx)){
        *last_val+=1;
        return true;
    }
        
    return false;
}


void initial_data(__data_type initial_time,
                IMU_data *imu_data,MAG_data *mag_data,
                GPS_data *gps_data,RPY_data *rpy_data,
                data0 *data_0) /*saida de dados*/
{
    

    if(!imu_data||!gps_data||!mag_data||!rpy_data||!data_0)
        alloc_error();

    //alocacoes em data_0
    data_0->rpy = matrix_m(1,3);
    data_0->acc = matrix_m(1,3);
    data_0->gyr = matrix_m(1,3);
    data_0->qtn = matrix_m(1,4);
    data_0->dcm = matrix_m(3,3);
    data_0->mag = matrix_m(1,3);
    data_0->bias_gyr = matrix_m(1,3);zeros(data_0->bias_gyr);
    data_0->bias_acc = matrix_m(1,3);zeros(data_0->bias_acc);
    
    // data_0->mode = matrix_m(1,15); //substitui "cont.mode"
    // campo gravitacional local com o drone parado
    int idx=0;
    int max_bound = imu_data->tim->n_row;
    
    while(imu_data->tim->elements[idx] <= initial_time){
        data_0->index = idx;
        idx = idx + 1;
    }
    

    mat_cpyRow( imu_data->acc, data_0->index, data_0->acc, 0);
    mat_cpyRow( imu_data->gyr, data_0->index, data_0->gyr, 0);
    
    debug(data_0->acc, 10);
    debug(data_0->gyr, 10);


    //No Matlab: Euler.tim == mag_data.tim
    int mgtim_sz = mag_data->tim->n_row;
    struct matrix *euler_tim = matrix_m(mgtim_sz, 1);
    mat_cpy(mag_data->tim, euler_tim);

    // data_0->euler_m = static_attitude(imu_data, mag_data);
    data_0->euler_m = read_matrix("DadosZurich/euler.csv");


    mat_cpyRow(data_0->euler_m, data_0->index, data_0->rpy, 0);
        
    //exit(0);
    // printf("data_0->euler_m\n");print_m(data_0->euler_m);
    
    //Dados iniciais do magnetometro
    idx=0;
    max_bound = mag_data->tim->n_row;
    while(LessThanOrEq(mag_data->tim->elements[idx], (__data_type) initial_time )){
        //if(idx < max_bound)
            data_0->last_mag = idx++;
        //else break;
    }
    
    // % Aceleração linear inicial
    mat_cpyRow(mag_data->mag, data_0->last_mag-1, data_0->mag, 0);
    //printf("Acel Lin. inicial\n");
    
    // printf("\ntest - %d\n", idx);
    // printf("mag %d x %d\n",mag_data->mag->n_row, mag_data->mag->n_col);
    // printf("dat %d x %d\n",data_0->mag->n_row, data_0->mag->n_col);

    // Dado inicial do GPS
    idx =0;
    max_bound = gps_data->tim->n_row;
    while(gps_data->tim->elements[idx] <= initial_time ){
        if(idx < max_bound)
            data_0->last_gps = idx++;
        else break;
    }
    // % Dado inicial do RPY
    idx=0;
    while(LessThanOrEq(rpy_data->tim->elements[idx], initial_time)){
            data_0->last_rpy = idx;
            idx = idx+1;
    }

    __data_type roll = angle2r(data_0->rpy->elements[0]);
    __data_type pitch= angle2r(data_0->rpy->elements[1]);
    __data_type yaw  = angle2r(data_0->rpy->elements[2]);
    
    angle2quat(yaw, pitch, roll, data_0->qtn);
    angle2dcm (yaw, pitch, roll, data_0->dcm); 
    
    // data_0->bias_gyr = (__data_type*) malloc(sizeof(__data_type)*15); //conforme definido em MATLAB

}

//#TODO: como corrigir esta função ? sol. provisoria !

//#TODO: como corrigir esta função ? sol. provisoria -> dbg_data_available
void data_available(int k, IMU_data *imu_data,MAG_data *mag_data,
                    GPS_data *gps_data,RPY_data *rpy_data,
                    int *last_gps,int *last_mag,int *last_rpy,
                    bool *ava_gps ,bool *ava_mag ,bool *ava_rpy){

    int max_mag =(int) fmax(mag_data->tim->n_row,mag_data->tim->n_col);//maxIdx_v(mag_data->tim->elements,mag_data->tim->n_row);//retorna o indice do ultimo elemento
    int max_gps =(int) fmax(gps_data->tim->n_row,gps_data->tim->n_col);//maxIdx_v(gps_data->tim->elements,gps_data->tim->n_row);//#TODO:retorna o valor: max(row,col)
    int max_rpy =(int) fmax(rpy_data->tim->n_row,rpy_data->tim->n_col);//maxIdx_v(rpy_data->tim->elements,rpy_data->tim->n_row);

    //last - ult indice disponivel do magnetometro/gps/angle
    
    bool mag = check_available(imu_data->tim, mag_data->tim, k, last_mag, max_mag);
    bool gps = check_available(imu_data->tim, gps_data->tim, k, last_gps, max_gps);
    bool rpy = check_available(imu_data->tim, rpy_data->tim, k, last_rpy, max_rpy);
    //printf("return: - %d\n",rpy);
    

    *ava_mag = mag;
    *ava_gps = gps;
    *ava_rpy = rpy;
    // printf("in - %d%d%d\n", mag,gps,rpy);

}