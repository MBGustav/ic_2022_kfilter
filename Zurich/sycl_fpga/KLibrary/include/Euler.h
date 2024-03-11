#ifndef EULER_H
#define EULER_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "libmat.h"
#include "structs.h"
#include "defines.h"
#include "err_handler.h"

// No matlab - > Euler.tim == mag_data.tim
void ned2enu(struct matrix *v, struct matrix *w);
struct matrix * static_attitude(IMU_data *imu_data, MAG_data *mag_data);


void quaternion_mul(struct matrix *ql, struct matrix *qr, struct matrix *Qout);
void quaternion_inv(struct matrix *ql, struct matrix *Qout);
void quatDisplay(struct matrix *Q);


#endif /*EULER_H*/