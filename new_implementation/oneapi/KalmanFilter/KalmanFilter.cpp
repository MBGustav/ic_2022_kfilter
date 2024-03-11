
#include "Matrix.hpp"
#include "system_define.h"
#include "Quaternion.hpp"
#include "Linalg.hpp"
#include "KalmanFilter.hpp"

KalmanFilter::~KalmanFilter(){this->deallocate();}

KalmanFilter::KalmanFilter(int dim_x, int dim_z)
{
    priori_state_covariance.resize(dim_x, dim_z);
    transition_matrix.resize(dim_x, dim_z);
    kalman_gain.resize(dim_x, dim_z);
    process_noise.resize(dim_x, dim_z);
    measurement_matrix.resize(dim_x, dim_z);
    measurement_noise.resize(dim_x, dim_z);
    measurement_function.resize(dim_x, dim_z);
    state_transition_matrix.resize(dim_x, dim_z);
    control_transition_matrix.resize(dim_x, dim_z);

    //Constant Values
    acce.resize(3,3);
    magn.resize(3,3);
    pos_gps.resize(3,3);
    vel_gps.resize(3,3);
    this->acce.fill_diag({acce_x,acce_y,acce_z});
    this->magn.fill_diag({magn_x,magn_y,magn_z});
    this->pos_gps.fill_diag({pos_gps_x,pos_gps_y,pos_gps_z});
    this->vel_gps.fill_diag({vel_gps_x,vel_gps_y,vel_gps_z});


}
void KalmanFilter::deallocate(){
    std::cout <<"dealloc\n";
}


// void KalmanFilter::set_constant_parameters(
//                     Matrix constant_Rz,
//                     Matrix constant_Rx){



void KalmanFilter::update()
{
    constexpr int vecSize_3d = 3;
    if(this->mode==0)
        propagated_quaternion.mem_copy_to(state_quaternion);
    else{
        Matrix d_x1(vecSize_3d,1);
        this->delta_x.mem_copy_to(d_x1, 0, 2);
        deltaAngle2quat_f(d_x1, propagated_quaternion, state_quaternion);
    }

    for(int i =0; i<  vecSize_3d; i++){
        //delta_x2
        state_bias_gyroscope.idx(i)    = propagated_bias_gyroscope.idx(i) + delta_x.idx(i+3);
        //delta_x3
        state_velocity.idx(i)          = propagated_velocity.idx(i) + delta_x.idx(i+6);
        //delta_x4
        state_bias_acceleration.idx(i) = propagated_bias_acceleration.idx(i) + delta_x.idx(i+9);
        //delta_x5
        state_position.idx(i)          = propagated_position.idx(i) + delta_x.idx(i+12);
    }
}


void KalmanFilter::step_kalman_gain()
{
    std_Matrix aux_prioriAndmeasure(1,1);/*PHt*/
    std_Matrix HPHt(1,1);
    std_Matrix aux_inverse(1,1);
    std_Matrix aux_kalmanAndmeasurement(1,1);
    if(mode == 0)
        kalman_gain.zeros();
    else{
        naive::gemm<false, true>(priori_state_covariance, measurement_function, aux_prioriAndmeasure, 1.0, 0.0);
        naive::gemm<false, true>(aux_prioriAndmeasure, measurement_function, HPHt, 1.0, 0.0);
        
        for(int i = 0; i < measurement_noise.getSize(); i++)
            aux_inverse.idx(i) = HPHt.idx(i) + measurement_noise.idx(i);

        naive::inverse(aux_inverse, aux_inverse);
        naive::gemm(aux_prioriAndmeasure, aux_inverse, kalman_gain, 1.0, 0.0);
        naive::gemm(kalman_gain, delta_z, delta_x, 1.0, 0.0);

    }
    naive::gemm(kalman_gain, measurement_function, aux_kalmanAndmeasurement, 1.0, 0.0);

    for(int i =0; i < kalman_gain.getRow(); i++)
        aux_kalmanAndmeasurement.at(i,i) = 1.0f - aux_kalmanAndmeasurement.at(i,i);

    naive::gemm(aux_kalmanAndmeasurement, propagated_P, priori_state_covariance, 1.0, 0.0);

}

void KalmanFilter::set_mode(int i){
    #ifndef JUMP_ASSERTION_COMANDS
    assert(i >= 0 && i < 7);
    #endif /*JUMP_ASSERTION_COMANDS*/
    this->mode  = i;
}

int KalmanFilter::get_mode(){return this->mode;}

void KalmanFilter::meas_mode00()
{
    delta_z.resize(3,1);
    measurement_function.zeros();
    measurement_noise.zeros();
}

void KalmanFilter::meas_mode01()
{
    std_Matrix rot(3, 3);
    std_Matrix subH1(3, 3);
    std_Matrix me(3, 3);
    std_Matrix skew_me(3, 3);

    delta_z.resize(3,1);

    quat2dcm(this->propagated_quaternion, rot);
    crossM_f(me, skew_me);

    naive::gemm(rot, skew_me, me, 1.0, 0.0); 

    
    

}
void KalmanFilter::meas_mode02()
{



}
void KalmanFilter::meas_mode03()
{

}
void KalmanFilter::meas_mode04()
{

}
void KalmanFilter::meas_mode05()
{

}
void KalmanFilter::meas_mode06()
{

}
void KalmanFilter::meas_mode07()
{

}
void KalmanFilter::meas_mode08()
{

}

// }

// KalmanFilter::KalmanFilter(
//                            Matrix measurement_noise, 
//                            Matrix measurement_matrix, 
//                            Matrix transition_matrix,
//                            Matrix initial_state_mean,
//                            Matrix process_noise
//                           )
// {
//     this->priori_state_covariance = priori_state_covariance;
//     this->transition_matrix = transition_matrix;
//     this->kalman_gain = kalman_gain;
//     this->process_noise = process_noise;
//     this->measurement_matrix = measurement_matrix;
//     this->measurement_noise = measurement_noise;
//     this->measurement_function = measurement_function;
//     this->state_transition_matrix = state_transition_matrix;
//     this->control_transition_matrix = control_transition_matrix;
// }

// void KalmanFilter::set_system_parameters(
//                                 Matrix lambda_gyr,
//                                 Matrix lambda_acc,
//                                 Matrix velocity_gps,
//                                 Matrix position_gps,
//                                 Matrix acceleration_bias,
//                                 Matrix magnetometer,
//                                 Matrix gyroscope)
// {
//     this->lambda_gyr = lambda_gyr;
//     this->lambda_acc = lambda_acc;
//     this->velocity_gps = velocity_gps;
//     this->position_gps = position_gps;
//     this->acceleration_bias = acceleration_bias;
//     this->magnetometer = magnetometer;;
//     this->gyroscope = gyroscope;
// }




// KalmanFilter::~KalmanFilter()
// {

// }






