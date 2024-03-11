#ifndef _KALMAN_FILTER_H_
#define _KALMAN_FILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <stdbool.h>

#include "libmat.h"
#include "structs.h"
#include "defines.h"
#include "Quaternions.h"

//set new state-space
state_space *NewStateSpace(); 

//kill state-space
void freeStateSpace(state_space *state);

//Defining Parameters
sys_param* Parameters();
//Set predefined matrices (see:parameters.h)
void setCovMatrices(struct matrix *Q, struct matrix*R);

void freeParameters(sys_param* sys);

//Propagation
void propagation(state_space *State, struct matrix* omega_fil,struct matrix* acc_fil, 
                 double time, struct matrix* ge, state_space *Prop);

//Cov Matrix Pred
void cov_matrix_serial(struct matrix* F, struct matrix* P, struct matrix* G, struct matrix* Q, struct matrix* propP);
//state_equations
void state_equations(sys_param *sys,state_space *prop,double time);


void data_reliable(int last_mag,bool ava_mag,state_space* prop,MAG_data *mag_data,double constant_gravity,double constant_magnetic,bool *rel_acc,bool *rel_mag);

int observation_modes(bool ava_mag, bool ava_gps, bool ava_bar,bool rel_acc, bool rel_mag);

struct matrix *measurement(sys_param *sys, int mode,state_space *prop,struct matrix* mag,
                 struct matrix* gps,struct matrix* ge,
                 struct matrix* me);


struct matrix* kalman_filter(int mode, state_space *prop, struct matrix* delta_z, sys_param *sys, int nx, struct matrix *P);

void update(state_space *prop,struct matrix* delta_x,int mode,state_space *state);/*return in State pointer*/





#endif//_KALMAN_FILTER_H_