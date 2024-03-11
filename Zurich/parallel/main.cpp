#include <sycl/sycl.hpp>
#include "oneapi/mkl.hpp"
#include "kernels/mkl_kernels.hpp"

#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <stdbool.h>

#include "structs.h"
#include "defines.h"
#include "libmat.h"
#include "Quaternions.h"
#include "file_reader.h"
#include "Initialization.h"
#include "kalman_filter.h"



//Ferramentas de Perfilemanto / debug
#include "Timer.hpp"




int idx=0; //global (debugging)

int main(){
    printf("Starting Execution\n");
    Timer timer, general_timer;

    timer_setconf("fun_time.out", &timer, "function,time_us\n");
    timer_setconf("time.out", &general_timer,"function,time_us\n");


    //setting initial parameters
    double initial_time = 5;
    double constant_gravity, constant_magnetic;
    // double roll, pitch, yaw;

    
    int nx = 15;
    // int *cont;//estaria relacionado com nx ? 
    int last_mag=0, last_gps=0, last_rpy=0;
    bool ava_mag = false;
    bool ava_rpy = false;
    bool ava_gps = false; 
    bool ava_bar = false; // #TODO: onde obter este dado ? 
    bool rel_acc = false; 
    bool rel_mag = false; 

   
    //Inicialização de estruturas
    /*% dados do giroscópio e acelerômetro*/
    dataset *Zurich = load_dataset();
    IMU_data *imu_data = Zurich->imu;
    MAG_data *mag_data = Zurich->mag;
    GPS_data *gps_data = Zurich->gps;
    RPY_data *rpy_data = Zurich->rpy;

    int N = imu_data->tim->n_row -1207;// Qtd. dados?
    // Inicializando Matrizes
    data0 *data_0 = (data0*) malloc(sizeof(data0));
    state_space *State = NewStateSpace(); 
    state_space *prop  = NewStateSpace();
    sys_param   *sys   = Parameters();

    prop->P = matrix_m(nx, nx);
    struct matrix* P = matrix_m(nx, nx);eye_m(1e-2, P);
    
    struct matrix* constant_rot = matrix_m(3,3);
    struct matrix* omegag = matrix_m(1, 3);
    struct matrix* mag = matrix_m(1,3);
    struct matrix* acc = matrix_m(1,3);
    struct matrix* me  = matrix_m(1,3);
    struct matrix* ge  = matrix_m(1,3);
    struct matrix* vel_gps = matrix_m(1,3);
    struct matrix* gps_pos = matrix_m(1, 3);//substitui: gps.pos
    struct matrix* prev_gps_pos = matrix_m(1,3);//substitui: prev.gps_pos
    struct matrix* delta_z; //= matrix_m(1,3);
    struct matrix* delta_x = NULL;

    struct matrix *time = matrix_m((N-200), 1);//atribuicao em loop
    struct matrix *input_vel_gps   = matrix_m(N-200, 3);//atribuicao em loop
    struct matrix *output_qtn      = matrix_m(N-200, 4);
    struct matrix *output_bias_gyr = matrix_m(N-200, 3);
    struct matrix *output_vel      = matrix_m(N-200, 3);
    struct matrix *output_bias_acc = matrix_m(N-200, 3);
    struct matrix *output_pos      = matrix_m(N-200, 3);
    struct matrix *output_rpy      = matrix_m(N-200, 3);
    struct matrix *output_acc      = matrix_m(N-200, 3);
    
    //matrizes auxiliares
    struct matrix *cont_mode= matrix_m(1 , 15); //#TODO: uso de array[15] simples?
    struct matrix *aux_rpy  = matrix_m(1 ,  3);

    angle2dcm((double) constant_beta, 0.0, 0.0, constant_rot);

    //definicoes locais
    constant_gravity  = local_gravity(imu_data, initial_time);
    constant_magnetic = local_magnetic(mag_data, initial_time);
    // setCovMatrices(sys_R, sys_Q);    ==> definido em Parameters

    dbg_cmd(printf("Setting initial_data\n"));
    //sycl::device dev = device_runner();

    initial_data(initial_time, imu_data, mag_data, gps_data, rpy_data, data_0);
    last_mag = (data_0->last_mag);
    last_gps = (data_0->last_gps);
    last_rpy = (data_0->last_rpy);
    
    mat_cpy(data_0->qtn, State->qtn);
    
    // Giroscopio
    mat_cpy(data_0->gyr, omegag);
    mat_cpy(data_0->bias_gyr, State->bias_gyr);

    //Acelerometro
    mat_cpy(data_0->bias_acc, State->bias_acc);
    mat_cpy(data_0->acc, acc);
    
    //Magnetometro
    timesc_m(1000, data_0->mag, mag);
        
    // % Posição e Velocidade via GPS
    mat_cpyRow(gps_data->pos, data_0->last_gps,State->pos, 0);
    mat_cpy(State->pos, gps_pos);
    zeros(vel_gps);zeros(State->vel);
    
    // % Intensidade do campo gravidade local e campo magnético local alinhado com o frame (NED) N: Norte Magnético
    dbg_cmd(printf("Setting me, ge and others for Loop...\n"));
    zeros(ge);ge->elements[2] = constant_gravity;
    zeros(me);me->elements[0] = norm_m(mag);

    int max_len = N;
    // int max_len = data_0->index + 500;
    
    for(int k= data_0->index+1; k < max_len; k++, idx++){
        timer_start(&general_timer);
        dbg_cmd(printf("Disponibilidade dos Dados Sensoriais\n"));

        dbg_profiling(timer_start(&timer));
        data_available(k,  imu_data, mag_data,
                       gps_data, rpy_data,
                       &last_gps, &last_mag, &last_rpy,/*return values*/
                       &ava_gps,  &ava_mag,  &ava_rpy);/*return values*/
        
        //% Velocidade angular no instante atual via giroscópio
        mat_cpyRow(imu_data->gyr, k, omegag, 0);
        //% Aceleração linear no instante atual via acelerômetro
        mat_cpyRow(imu_data->acc, k, acc, 0);
        if(ava_mag){
            mag->elements[0] = T2G * mag_data->mag->elements[last_mag*3 + 0];
            mag->elements[1] = T2G * mag_data->mag->elements[last_mag*3 + 1];
            mag->elements[2] = T2G * mag_data->mag->elements[last_mag*3 + 2] - mag_data->mag->elements[0*3 + 2];
        }
        dbg_profiling(timer_stop(&timer, "data_available"));



        //  %% Propagação dos Estados e Medidas Sensoriais
        dbg_profiling(timer_start(&timer));
        
        
        
        propagation(State, omegag, acc, T_imu, ge, prop);


        if(ava_gps){    
            // % Valor da posição no instante anterior via GPS
            mat_cpy(gps_pos, prev_gps_pos);
            // % Posição no instante atual via GPS
            mat_cpyRow(gps_data->pos, last_gps, gps_pos ,0);
            // % Cálculo da velocidade a partir dos dados de posição via GPS com dt = 0.1s
            less_m(gps_pos, prev_gps_pos, vel_gps);
            timesc_m((freq_gps), vel_gps, vel_gps);
        }
        // % Salva valores de velocidade com filtro via GPS no vetor
        mat_cpyRow(vel_gps, 0, input_vel_gps, idx);
        dbg_profiling(timer_stop(&timer, "Propagation"));

        // %% Representação em espaço de estados em tempo discreto
        dbg_cmd(printf("State equations step\n"));
        dbg_profiling(timer_start(&timer));
        state_equations(sys, prop, T_imu);
        dbg_profiling(timer_stop(&timer, "State Equations step"));


        // %% Predição da matriz de covariância
        dbg_cmd(printf("Covariance Matrix Prediction\n"));
        dbg_profiling(timer_start(&timer));
        // oneapi_cov_matrix_pred(dev, sys->F, P, sys->G, sys->Q, prop->P);
        cov_matrix_serial(sys->F, P,sys->G, sys->Q, prop->P);
        dbg_profiling(timer_stop(&timer, "Cov. Matrix Prediction"));


        dbg_cmd(printf("Data Reliable\n"));
        dbg_profiling(timer_start(&timer));
        data_reliable(last_mag, ava_mag, prop, mag_data, constant_gravity, constant_magnetic, &rel_acc, &rel_mag);
        dbg_profiling(timer_stop(&timer, "data_reliable"));

        dbg_cmd(printf("Observation Nodes\n"));
        dbg_profiling(timer_start(&timer));
        int mode = observation_modes(ava_mag, ava_gps, ava_bar/*?*/, rel_acc, rel_mag);
        dbg_profiling(timer_stop(&timer, "Observation Nodes "));

        dbg_cmd(printf("Measurement\n"));
        dbg_profiling(timer_start(&timer));
        delta_z = measurement(sys, mode, prop, mag, gps_pos, ge, me);
        dbg_profiling(timer_stop(&timer, "Measurement"));
        
        dbg_cmd(printf("Kalman Filter\n"));
        dbg_profiling(timer_start(&timer));
        delta_x = kalman_filter(mode, prop, delta_z, sys, nx, P);
        dbg_profiling(timer_stop(&timer, "Kalman Filter"));
        
        // %% Atualização dos Estados
        
        dbg_profiling(timer_start(&timer));
        update(prop, delta_x, mode, State);/*return in State pointer*/

        dbg_profiling(timer_stop(&timer, "Update"));
        mat_cpyRow(State->bias_gyr, 0, output_bias_gyr  , idx);
        mat_cpyRow(State->vel     , 0, output_vel       , idx);
        mat_cpyRow(State->bias_acc, 0, output_bias_acc  , idx);
        mat_cpyRow(State->pos     , 0, output_pos       , idx);
        mat_cpyRow(State->qtn     , 0, output_qtn       , idx);
        quat2angle(State->qtn, aux_rpy);
        mat_cpyRow(aux_rpy        , 0, output_rpy       , idx);

        if(delta_z != NULL) delete_m(delta_z);
        // create a new matriz on each loop(measurement)

        timer_stop(&general_timer, "Tempo por Loop");
    }

    dbg_cmd(printf("Saving Outputs\n"));
    save_matrix(output_bias_gyr, "output/output_bias_gyr.csv");
    save_matrix(output_vel,      "output/output_vel.csv");
    save_matrix(output_bias_acc, "output/output_bias_acc.csv");
    save_matrix(output_pos,      "output/output_pos.csv" );
    save_matrix(output_qtn,      "output/output_qtn.csv" ) ;
    save_matrix(output_rpy,      "output/output_rpy.csv");





    dbg_cmd(printf("Freeing memory\n"));
    delete_m(me);
    delete_m(ge);
    delete_m(constant_rot);
    delete_m(P);
    delete_m(omegag);
    delete_m(acc);
    delete_m(mag);
    // delete_m(vel_gps);

    //Output matrices  
    delete_m(time);
    delete_m(input_vel_gps);
    delete_m(output_qtn);
    delete_m(output_bias_gyr);
    delete_m(output_vel);
    delete_m(output_bias_acc);
    delete_m(output_pos);
    delete_m(output_rpy);
    // delete_m(mode);

    //freeing aux matrices
    delete_m(cont_mode);
    delete_m(aux_rpy);
    //AND more Frees ....
    freeIMU(imu_data);
    freeMAG(mag_data);
    freeGPS(gps_data);
    freeRPY(rpy_data);
    freeStateSpace(State);
    freeStateSpace(prop);

    printf("End of Test\n");

    return 0;
}
