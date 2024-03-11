
#ifndef __SYSTEM_DEFINE_H__
#define __SYSTEM_DEFINE_H__


#define PI 3.14159265358979
#define TINY (1.0e-200)
#define GRAVITY -9.81
#define Time  1/freq   // s
#define Tau  100  // s
#define T2G (1) /*conversao de unidade: fluxo magnético de Tesla para Gauss */
#define constant_beta (-4.0 * (angle2rad))

#define sin_d(x) (sin((x)*angle2rad))
#define cos_d(x) (cos((x)*angle2rad))
#define atan_2d(x,y) (atan2(angle2r(x),angle2r(y)))
#define angle2r(x) ((x) * (PI/180))
#define r2angle(x) ((x) * (180/PI))

#define FREQUENCY  25 /*Hz*/
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


#endif /*__SYSTEM_DEFINE_H__*/