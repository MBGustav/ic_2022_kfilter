#ifndef _QUATERNIONS_H_
#define _QUATERNIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libmat.h"
#include "structs.h"
#include "Euler.h"

//Declaração de Estrutura de Quaternion -> usando matriz (4 x 1)

void quaternion_mul(struct matrix *ql, struct matrix *qr, struct matrix *Qout);
void quaternion_inv(struct matrix *ql, struct matrix *Qout);
__data_type quatnormalize(struct matrix *Q);
void quatDisplay(struct matrix *Q);

void quat2dcm(struct matrix *Q, struct matrix* dcm);
void angle2dcm(__data_type z,__data_type y, __data_type x, struct matrix* dcm_out);
void angle2quat(__data_type z,__data_type y, __data_type x, struct matrix* quat);

void quat2angle(struct matrix* Q, struct matrix *angle);
void deltaAngle2quat_f(struct matrix* deltaAngle, struct matrix* LGqpHat, struct matrix* out_quaternion);






#endif /*_QUATERNIONS_H_*/