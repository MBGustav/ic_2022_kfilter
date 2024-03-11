// Defines
#ifndef _DEFINES_H_
#define _DEFINES_H_

// Inclusao de constantes de constant.m
#define PI 3.14159265358979
#define TINY (1.0e-200)
#define f  25     // Hz
#define g -9.81
#define T  1/f   // s
#define tau  100  // s
#define T2G (1) /*conversao de unidade: fluxo magnético de Tesla para Gauss */
#define constant_beta (-4.0 * (angle2rad))


#define sind(x) (sin((x)*angle2rad))
#define cosd(x) (cos((x)*angle2rad))
#define atan2d(x,y) (atan2(angle2r(x),angle2r(y)))

#define zeros(Matrix_in) (timesc_m(0.0, Matrix_in, Matrix_in))
#define angle2rad (PI/180)
#define rad2angle (180/PI)
#define angle2r(x) ((x) * (PI/180))
#define r2angle(x) ((x) * (180/PI))
#define EPSILON (0.000000005f)
#define IsEqual(x,y) (( ((x)-(y))*((x)-(y))  <= EPSILON * 10000)? 1 : 0)


#define LEFT_INTERVAL  0.00000000001
#define RIGHT_INTERVAL 0.000000001

#define FLOAT_COMPARE_WITH_INTERVAL(left, right) \
    ((fabs((right) - (left)) <= LEFT_INTERVAL) || (fabs((left) - (right)) <= RIGHT_INTERVAL))
    
#define FLT_EQ(x,y) (( ((x)-(y))*((x)-(y))  <= EPSILON)? 1 : 0)

// ############################### Parameters  Definition ##########################
//define frequencias
#define freq_imu (250)
#define freq_mag (83.33)
#define freq_gps (2.5)
#define T_imu (0.004)

//Definition gyro Variances
#define  gyro_x 1e-7
#define  gyro_y 1e-7
#define  gyro_z 1e-7
#define  gyro_bias_x 0.0000000001
#define  gyro_bias_y 0.0000000001
#define  gyro_bias_z 0.0000000001

//Definition Accelerometer Variances
#define  acce_x 0.01
#define  acce_y 0.01
#define  acce_z 0.01
#define  acce_bias_x 0.000000001
#define  acce_bias_y 0.000000001
#define  acce_bias_z 0.000000001

//Define Magnetometer Variances
#define  magn_x 0.1
#define  magn_y 0.1
#define  magn_z 0.1

//Define gps Variances
#define  pos_gps_x 0.0001
#define  pos_gps_y 0.0001
#define  pos_gps_z 0.0001
#define  vel_gps_x 1
#define  vel_gps_y 1
#define  vel_gps_z 1 

// ############################### Ellipsoid Definition ##########################

// % Informações do elipsóide de referência para o Sistema Geodésico Mundial de 
// % 1984 com unidade em metros
#define wg84_SemimajorAxis (6378137) //(meters)
#define wg84_Eccentricity  (0.0818191908426215) //(meters)

#ifdef _DEBUG
    
    #define debug(XX, a)  {\
        if(XX==NULL) printf("%s is null!\n", #XX); \
        else{                       \
            printf("(%s) \n", #XX); \
            mat_size(XX);           \
            print_m(XX, a);         \
        }}            

    #define dbg_cmd(cmd) do { \
        cmd; \
    } while (0)

    #define dbg_eq(X, Y) printf("(%s == %s) ? ", #X,#Y); \
                     printf(matrix_isEqual((X), (Y))? "  equal\n": " not equal\n")
#else
    #define debug(XX, a) do{}while(0)
    #define dbg_cmd(cmd) do{}while(0)
    #define dbg_eq(X, Y) do{}while(0)

#endif

#ifdef _PROFILING
    #define dbg_profiling(cmd) do{ \
        cmd; \
    }while(0)
#else
    #define dbg_profiling(cmd) do{}while(0)

#endif 



#endif /*_DEFINES_H_*/