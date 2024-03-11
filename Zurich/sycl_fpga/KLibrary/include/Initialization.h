#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "libmat.h"
#include "structs.h"
#include "defines.h"



//Set predefined matrices (see:defines.h)
void setCovMatrices(struct matrix *Q, struct matrix*R);

sys_param* Parameters();

//set new state-space
state_space *NewStateSpace(); 

//kill state-space
void freeStateSpace(state_space *state);

//initi
void initial_data(double initial_time,
                IMU_data *imu_data,MAG_data *mag_data,
                GPS_data *gps_data,RPY_data *rpy_data,
                data0 *data_0); 
                //#TODO: preciso do vetor cont ?

void data_available(int k, IMU_data *imu_data,MAG_data *mag_data,
                    GPS_data *gps_data,RPY_data *rpy_data,
                    int *last_gps,int *last_mag,int *last_rpy,
                    bool *ava_gps ,bool *ava_mag ,bool *ava_rpy);

//function to check if data is available: rpy_data, gps_data or mag_data
bool check_available(struct matrix *imu_time, struct matrix *par_time,
                         int k, int *p_lastVal, int max_idx);

bool LessThanOrEq(double x, double y);


#endif//_INITIALIZATION_H