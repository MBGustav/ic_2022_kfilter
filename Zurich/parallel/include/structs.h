// Structs

#ifndef STRUCTS_H
#define STRUCTS_H

#include "defines.h"
struct matrix {
        __data_type *elements;
        int n_row;
        int n_col;
#ifdef MKL_ONEAPI
        int lda;
        
#endif
};
//################################## IMU's structs #############################
typedef struct axes {
        __data_type x;
        __data_type y;
        __data_type z;
}axes;
// Data from IMU
struct IMU{
        axes acel;
        axes rate;
        axes mag;
};

struct data_IMU {
        struct IMU bin;
        struct IMU scaled;
        struct IMU adj;
};

typedef struct euler {
        __data_type yaw;
        __data_type pitch;
        __data_type roll;
}euler;




// ############################### State Space's structs ###############################
// Weight matrices
struct wm {
        struct matrix *Q;
        struct matrix *Qeta;
        struct matrix *Qeta_sr;
        struct matrix *R;
        struct matrix *R_sr;
        struct matrix *P;
        struct matrix *P_sr;
        struct matrix *Pp;
};
// Discrete space state
struct ss_d {
        struct matrix *x;
        struct matrix *xp;
        struct matrix *z;
        struct matrix *Phi;
        struct matrix *G;
        struct matrix *H;
        struct wm wm;
};


// ############################### Zurich Structs ##########################


//#TODO: viavel? prop / state
typedef struct state_space
{
        //used in state and prop
        struct matrix *bias_gyr;
        struct matrix *bias_acc;
        
        struct matrix *qtn;
        struct matrix *pos;
        struct matrix *vel;
        struct matrix *acc;

        struct matrix *omegag;

        //used in prop ?
        struct matrix *P;

}state_space;

typedef struct sys_param{
        struct matrix *F;
        struct matrix *G;
        struct matrix *H;
        struct matrix *lambda_gyr;
        struct matrix *lambda_acc;
        struct matrix *R;
        struct matrix *Q;
}sys_param;


typedef struct data0{
        int index;

        int last_mag, last_gps, last_rpy;
        // axes mag;
        // axes acc;
        // axes gyr;
        struct matrix *mag;
        struct matrix *gyr;
        struct matrix *acc;

        struct matrix *dcm;
        struct matrix *qtn;
        struct matrix *rpy;
        struct matrix *bias_gyr;
        struct matrix *bias_acc;
        struct matrix *euler_m;
}data0;



//Definicao de Estrutura de Dados - dataset Zurich

typedef struct GPS_data{
        struct matrix* tim;
        struct matrix* pos;

        /*considering:
          x = lat; 
          y = lon;
          z = alt; 
        */
}GPS_data;
typedef struct pos_vel {
        struct matrix* tim;
        struct matrix *coord;//#TODO:converter para matrizes?
        
        /*considering:
          x = east; 
          y = north;
          z = up; 
        */
}pos_vel;
typedef struct bar_data{
        struct matrix* tim; 
        struct matrix* alt;
}BAR_data;
typedef struct IMU_data{
        struct matrix* tim; 
        struct matrix *gyr;
        struct matrix *acc;
}IMU_data;

typedef struct MAG_data{ 
        struct matrix* tim;
        struct matrix *mag;
}MAG_data;

typedef struct RPY_data{
        struct matrix* tim;
        struct matrix *rpy;

}RPY_data;

typedef struct dataset{
        BAR_data *bar;
        IMU_data *imu;
        MAG_data *mag;
        RPY_data *rpy;
        GPS_data *gps;
}dataset;


// #############################################################################



#endif /*STRUCTS_H*/