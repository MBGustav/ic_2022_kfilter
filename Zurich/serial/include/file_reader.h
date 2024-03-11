#ifndef _FILE_READER_H
#define _FILE_READER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "libmat.h"
#include "structs.h"
#include "err_handler.h"

struct matrix* read_file_nD(char name[], int nD);

IMU_data *import_IMU();
GPS_data *import_GPS();
GPS_data *import_GPS2(); //for debug
RPY_data *import_RPY();
MAG_data *import_MAG();

dataset *load_dataset();

void freeMAG( MAG_data *data );
void freeIMU( IMU_data *data );
void freeGPS( GPS_data *data );
void freeRPY( RPY_data *data );


double local_magnetic(MAG_data *mag_data, double initial_time);
double local_gravity (IMU_data *imu_data, double initial_time);



#endif //_FILE_READER_H