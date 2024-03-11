
#include "kalman_filter.h"


//Usado para vetor state e prop
state_space *NewStateSpace(){

    state_space *out = (state_space*)malloc(sizeof(state_space));
    
    if (out == NULL) alloc_error();

    out->bias_gyr = matrix_m(1, 3); zeros(out->bias_gyr); //roll, pitch, yaw
    out->bias_acc = matrix_m(1, 3); zeros(out->bias_acc); //x, y, z
    out->qtn = matrix_m(1, 4); zeros(out->qtn); // w, x, y, z
    out->acc = matrix_m(1, 3); zeros(out->acc); //x,y,z 
    out->vel = matrix_m(1, 3); zeros(out->vel); //x,y,z 
    out->pos = matrix_m(1, 3); zeros(out->pos); //x,y,z

    out->omegag = matrix_m(1, 3);  zeros(out->omegag);
    out->P = NULL;//#TODO: to be defined .. ?

    return out;
}


void freeStateSpace(state_space *state) {
    
    if(state == NULL) return;

    free(state->bias_gyr);
    free(state->bias_acc);
    free(state->qtn);
    free(state->pos);
    free(state->vel);
    free(state->acc);
    if ((state->P != NULL)) free(state->P);
    
    free(state);
}


void propagation(state_space *State, struct matrix* omega_fil,struct matrix* acc_fil, 
                 double time, struct matrix* ge, state_space *Prop){
    
    if(omega_fil==NULL ||acc_fil==NULL ||Prop==NULL) alloc_error();

    struct matrix *qtn_transp= matrix_m(4,1);
    struct matrix *M = matrix_m(3, 3); 
    struct matrix *M_t = matrix_m(3, 3); 
    struct matrix *w = matrix_m(1, 3);
    struct matrix *Omega = matrix_m(4, 4);
    struct matrix *Omega_aux = matrix_m(4, 4);

    struct matrix *acc_t= matrix_m(3, 1); 
    
    struct matrix *auxM1 = matrix_m(4,4); 
    struct matrix *f_m = matrix_m(1,3); 
    struct matrix *f_mt = matrix_m(3, 1); 
    struct matrix *auxM3 = matrix_m(1,3);
    double *pw, *p_Omega;
    double nw = 0.0f; 
    double w1, w2,w3;

    //  propaga state

    mat_cpy(State->bias_gyr, Prop->bias_gyr);
    less_m(omega_fil, Prop->bias_gyr, Prop->omegag);

    // % Propagação do Bias da aceleração linear
    mat_cpy(State->bias_acc, Prop->bias_acc);
    // % Aceleração linear corrigida
    less_m(acc_fil, Prop->bias_acc, Prop->acc);
    
    double scal = (0.5*time);
    timesc_m(scal, Prop->omegag, w);

    p_Omega = Omega->elements;
    pw = w->elements;
    w1 = pw[0];
    w2 = pw[1];
    w3 = pw[2];   

    p_Omega[0] = 0 ; p_Omega[1] = -w1; p_Omega[2] = -w2; p_Omega[3] = -w3;
    p_Omega[4] = w1; p_Omega[5] =  0 ; p_Omega[6] =  w3; p_Omega[7] = -w2;
    p_Omega[8] = w2; p_Omega[9] = -w3; p_Omega[10]=  0 ; p_Omega[11]=  w1;
    p_Omega[12]= w3; p_Omega[13]=  w2; p_Omega[14]= -w1; p_Omega[15]=  0;

    nw = norm_m(w);
    
    // qtn = (cos(norm_m(w)) * eye(4) + (sin(norm_m(w))/norm_m(w)) * OMEGA) * state.qtn;
    //TODO: eu consigo reduzir o nro de operações alterando auxM1 (transposicao) ?
    
    eye_m(cos(nw), auxM1); 
    
    timesc_m(sin(nw)/(nw), Omega, Omega_aux); // OMEGA = (sin(norm_m(w))/norm_m(w)) * OMEGA
    sum_m(auxM1, Omega_aux, auxM1); // auxM1 = cos(norm_m(w)) * eye(4) + (sin(norm_m(w))/norm_m(w)) * OMEGA

    transp_m(State->qtn, qtn_transp);
    times_m(auxM1, qtn_transp, qtn_transp);
    transp_m(qtn_transp, Prop->qtn);

    // % Conversão para matriz de rotação
    quat2dcm(State->qtn, M);
    transp_m(M, M_t);


    // % Cálculo da aceleração cinemática    
    transp_m(Prop->acc, acc_t);
    times_m(M_t, acc_t, f_mt);
    transp_m(f_mt, f_m);
    sum_m(ge, f_m, f_m);
    
    // % Propagação da velocidade
    timesc_m(time, f_m, Prop->vel);
    sum_m(Prop->vel, State->vel, Prop->vel);
    

    // % Propagação da posição
    timesc_m(time, Prop->vel, auxM3);
    sum_m(auxM3, State->pos, Prop->pos);

    

    //free memory
    delete_m(Omega);
    delete_m(Omega_aux);
    delete_m(qtn_transp);
    delete_m(M);
    delete_m(M_t);
    delete_m(acc_t);
    delete_m(w);
    delete_m(auxM1);
    delete_m(f_m);
    delete_m(auxM3);
}


inline void cov_matrix_serial(struct matrix* F, struct matrix* P, struct matrix* G, struct matrix* Q, struct matrix* propP)
{
    int Gr = G->n_row;
    int Gc = G->n_col;
    int Qr = Q->n_row;
    int Qc = Q->n_col;
    int Fr = F->n_row;
    int Fc = F->n_col;
    int Pr = P->n_row;
    int Pc = P->n_col;
    
    struct matrix* Gt = matrix_m(Gc, Gr);
    struct matrix* Ft = matrix_m(Fc, Fr);
    struct matrix* FP = matrix_m(Fr, Pc);
    struct matrix* GQ = matrix_m(Gr, Qc);
    struct matrix* FPF = matrix_m(Fr, Pc);
    struct matrix* GQG = matrix_m(Fr, Pc);
    // struct matrix* propP = matrix_m(Fr, Pc);
    
    // FPF  = F * P * F'
    times_m(F, P, FP);      
    transp_m(F, Ft);
    times_m(FP, Ft, FPF);
    
    // GQG = G x Q x G'
    times_m(G, Q, GQ);
    transp_m(G, Gt);
    times_m(GQ, Gt, GQG);

    sum_m(FPF, GQG, propP); 

    delete_m (Gt);
    delete_m (Ft);
    delete_m (FP);
    delete_m (GQ);
    delete_m (FPF);
    delete_m (GQG);
    // return propP;
}

sys_param* Parameters(){

    int nx = 15;
    int nw = 12;
    int nv =  3;
    sys_param* sys_out = (sys_param*) malloc(sizeof(sys_param));
    
    if (!sys_out) alloc_error();
    
    sys_out->R = matrix_m(nw, nw);// Q = matrix_m(12, 12);
    sys_out->Q = matrix_m(nw, nw);// R = matrix_m(12, 12);
    sys_out->lambda_gyr = matrix_m(3, 3);
    sys_out->lambda_acc = matrix_m(3, 3);


    //#TODO:qual o tamanho destas matrizes? Generalizado ===========
    sys_out->F = matrix_m(nx, nx);
    sys_out->G = matrix_m(nx, nw);
    sys_out->H = NULL; // matrix_m( 1, nx); // confirmar tamanho
// ====================================================
    setCovMatrices(sys_out->Q, sys_out->R);
    
    double *p_lg = sys_out->lambda_gyr->elements;
    double *p_la = sys_out->lambda_acc->elements;

    // diag([200; 300; 100])
    zeros(sys_out->lambda_gyr);
    p_lg[0] = 200;p_lg[4] = 300;p_lg[8] = 100;

    // diag([1/100; 1/100; 1/100]);
    eye_m(0.01, sys_out->lambda_acc);
    

    return sys_out;
}

void freeParameters(sys_param* sys){
    if( sys == NULL )
        return;
    delete_m( sys->R) ;
    delete_m( sys->Q) ;
    delete_m( sys->lambda_gyr) ;
    delete_m( sys->lambda_acc) ;
    delete_m( sys->F) ;
    delete_m( sys->G) ;
    if(sys->H !=NULL)  delete_m( sys->H) ;

    free(sys);
    
}

//modifcado - altera por meio de ponteiros
void state_equations(sys_param *sys, state_space *prop, double time){
    int nv = 3,
        nx = 15,
        nw = 12;
    //#TODO: alterar modo de cópia da matriz A -> nao precisaria transpor
    struct matrix *minusRT = matrix_m(nv, nv);
    struct matrix *minusRT_aux = matrix_m(nv, nv);
    struct matrix *OmegaCross = matrix_m(nv, nv);
    struct matrix *AccelCross = matrix_m(nv, nv);
    struct matrix *I = matrix_m(nv, nv);
    struct matrix *I_15 = matrix_m(nx, nx);//ser usado no passo sist. dinamico 
    struct matrix *n_I = matrix_m(nv, nv); //negative I
    struct matrix *A = matrix_m(nx, nx);zeros(A);
    struct matrix *A_t = matrix_m(nx, nx);
    struct matrix *B = matrix_m(nx, nw);zeros(B);
    struct matrix *n_lambda_gyr = matrix_m(nv, nv);
    struct matrix *n_lambda_acc = matrix_m(nv, nv);
    struct matrix *auxM1 = matrix_m(nv, nv); // aux: minusRT * ACCELcross
    eye_m(1.0 , I_15);
    eye_m(1.0 , I);
    eye_m(-1.0, n_I);
    // % Matriz de rotação do frame Global para o frame Local
    quat2dcm(prop->qtn, minusRT);
    transp_m(minusRT, minusRT_aux);
    timesc_m(-1.0, minusRT_aux, minusRT);


    // % Matriz antissimétrica usada na atualização do quatérnio de atitude
    crossM_f(prop->omegag, OmegaCross);
    timesc_m(-1, OmegaCross, OmegaCross);//sera usado "A", negativo.

    // % Matriz antissimétrica usada na atualização da velocidade linear
    crossM_f(prop->acc, AccelCross);
    
    //obtendo aux: minusRT * ACCELcross
    times_m(minusRT, AccelCross, auxM1);
    //negative val from lambda_acc
    timesc_m(-1.0, sys->lambda_acc, n_lambda_acc);
    //negative val from lambda_gyr
    timesc_m(-1.0, sys->lambda_gyr, n_lambda_gyr);

    // % Preenchimento da matriz de estados A
    getpart_mat(OmegaCross  ,  0, 0, A);//A(0:2,0:2) = -OMEGAcross;
    getpart_mat(n_I         ,  0, 3, A);//A(0:2,3:5) = -I;                   
    getpart_mat(n_lambda_gyr,  3, 3, A);//A(3:5,3:5) = -system.lambda.gyr; 
    getpart_mat(auxM1       ,  6, 0, A);//A(6:8,0:2) = minusRT * ACCELcross;
    getpart_mat(minusRT     ,  6, 9, A);//A(6:8,9:11) = minusRT;
    getpart_mat(n_lambda_acc,  9, 9, A);//A(9:11,9:11) = -system.lambda.acc;
    getpart_mat(I           , 12, 6, A);//A(13:15,6:8) = I;

    // % Preenchimento da matriz de ruídos do processo
    getpart_mat(n_I, 0, 0, B);
    getpart_mat(I, 3, 3, B);
    getpart_mat(minusRT, 6, 6, B);
    getpart_mat(I, 9, 9, B);
    // % Discretização do sistema dinâmico: sys->F
    eye_m(1, I_15);
    timesc_m(time, A, A);
    
    sum_m(A, I_15, sys->F);


    // % Discretização do sistema dinâmico: sys->G
    double scal = sqrt(time);
    timesc_m(scal, B, sys->G);
    

    //free memory
    delete_m(minusRT);
    delete_m(OmegaCross);
    delete_m(AccelCross);
    delete_m(I);
    delete_m(I_15);
    delete_m(n_I);
    delete_m(A);
    delete_m(B);
    delete_m(n_lambda_gyr);
    delete_m(n_lambda_acc);
    delete_m(auxM1);
}

int observation_modes(bool ava_mag, bool ava_gps, bool ava_bar,bool rel_acc, bool rel_mag){
    // % Descrição: Determina o modo de markov para a observação dos sensores.
    // % 0000) ---/---/---/---
    // % 0001) ---/---/---/mag
    // % 0010) ---/---/acc/---
    // % 0011) ---/---/acc/mag
    // % 0100) ---/gps/---/---
    // % 0101) ---/gps/---/mag
    // % 0110) ---/gps/acc/---
    // % 0111) ---/gps/acc/mag
    
    //#TODO: nao esta sendo usado o quarto digito ??
    // % 1000) slm/---/---/---
    // % 1001) slm/---/---/mag
    // % 1010) slm/---/acc/---
    // % 1011) slm/---/acc/mag
    // % 1100) slm/gps/---/---
    // % 1101) slm/gps/---/mag
    // % 1110) slm/gps/acc/---
    // % 1111) slm/gps/acc/mag
    int mode,
        mode_mag=0, mode_gps=0, mode_acc=0;
    
    if (ava_mag && rel_mag) mode_mag=1;
    if(rel_acc) mode_acc =1;
    if(ava_gps) mode_gps =1;
    mode = 4*mode_gps + 2*mode_acc + mode_mag;

    return mode;
}

void data_reliable(int last_mag,bool ava_mag,state_space* prop,MAG_data *mag_data,double constant_gravity,double constant_magnetic,bool *rel_acc,bool *rel_mag){

    double rho_acc = 0.05;
    double rho_mag = 1000000.0f;

    // % Teste de influência de aceleração cinemática na medida do acelerômetro    
    double lim_acc = norm_m(prop->acc) / constant_gravity;
    double lim_mag;
    struct matrix *mag = matrix_m(1,3);
    
    
    if(fabs(lim_acc-1) < rho_acc) *rel_acc = true;
    else *rel_acc = false; 

    if(ava_mag){
        // % Utiliza-se (last.mag-1) pois o valor do índice de mag foi atualizado anteriormente pela função data_available
        mat_cpyRow(mag_data->mag, (last_mag-1),mag, 0);
        lim_mag =  norm_m(mag) / constant_magnetic;
        if(fabs(lim_mag/constant_magnetic -1.0f) < rho_mag) *rel_mag = true;
        else *rel_mag = false;
    } else {
        *rel_mag = false;
        }

}

//Importa matrix R e Q, predefinida no arquivo parameter.h
void setCovMatrices(struct matrix *Q, struct matrix *R){

    // Q = matrix_m(12, 12);
    // R = matrix_m(12, 12);
    
    //Ruido do Processo(diagonal)
    Q->elements[12*0 +0] = gyro_x;
    Q->elements[12*1 +1] = gyro_y;
    Q->elements[12*2 +2] = gyro_z;
    
    Q->elements[12*3 +3] = gyro_bias_x;
    Q->elements[12*4 +4] = gyro_bias_y;
    Q->elements[12*5 +5] = gyro_bias_z;

    Q->elements[12*6 +6] = acce_x;
    Q->elements[12*7 +7] = acce_y;
    Q->elements[12*8 +8] = acce_z;

    Q->elements[12*9 +9 ] = acce_bias_x;
    Q->elements[12*10+10] = acce_bias_y;
    Q->elements[12*11+11] = acce_bias_z;

    //Ruido das medidas(diagonal)
    R->elements[12*0 +0] = magn_x;
    R->elements[12*1 +1] = magn_y;
    R->elements[12*2 +2] = magn_z;

    R->elements[12*3 +3] = acce_x;
    R->elements[12*4 +4] = acce_y;
    R->elements[12*5 +5] = acce_z;

    R->elements[12*6 +6] = vel_gps_x;
    R->elements[12*7 +7] = vel_gps_y;
    R->elements[12*8 +8] = vel_gps_z;

    R->elements[12*9 +9 ] = pos_gps_x;
    R->elements[12*10+10] = pos_gps_y;
    R->elements[12*11+11] = pos_gps_z;
}

//#TODO: divisao de casos - elaborar testes com entradas fixas para os casos abaixo
static inline struct matrix* 
        meas_mode00(sys_param *sys, state_space *prop){
    
    struct matrix *delta_z = matrix_m(3, 1);
    struct matrix* H = matrix_m(1,15);
    struct matrix* R = matrix_m(3,15);
    zeros(H);
    zeros(R);

    sys->H = H;
    sys->R = R;
    
    return NULL;
}

static inline struct matrix* 
        meas_mode01(sys_param *sys, state_space *prop, struct matrix *me, struct matrix *magn, struct matrix *mag){
    
    int Hr = 3, Hc = 15;
    int Rr = 3, Rc = 3;
    //saida matriz
    struct matrix *delta_z = matrix_m(3, 1);
    struct matrix *rot = matrix_m(3, 3);
    struct matrix *subH1= matrix_m(3, 3);
    struct matrix *skew_me= matrix_m(3, 3);
    struct matrix *z_prop = matrix_m(3, 1);
    // struct matrix *z_prop = matrix_m(3, 1);
    struct matrix *me_tr = matrix_m(3,1);
    struct matrix *mag_t = matrix_m(3,1);
    struct matrix *ge_tr = matrix_m(3,1);

    struct matrix *H = matrix_m(Hr, Hc);
    struct matrix *R = matrix_m(Rr, Rc);
    
    // rot = quat2dcm(prop.qtn');
    quat2dcm(prop->qtn,rot);
    //subH1 = rot * skew_me
    crossM_f(me, skew_me);
    times_m(rot,skew_me, subH1);
    
    //Preenchimento H e R
    zeros(H);
    getpart_mat(subH1, 0, 0, H);

    zeros(R);
    mat_cpy(magn, R);

    
    //z = mag; 
    //z_prop = rot*me;
    //delta_z = z - z_prop;
    transp_m(me, me_tr);
    times_m (rot, me_tr, z_prop);//z_prop = rot*me == (3,1)
    transp_m(mag, mag_t);
    less_m(mag_t, z_prop, delta_z);
    // debug(z, 1000);
    // debug(z_prop, 1000);

    sys->H = H;
    sys->R = R;
    // debug(me,1);
    // debug(z_prop,1);
    // free aux memory
    delete_m(rot);
    delete_m(subH1);
    delete_m(skew_me);
    delete_m(z_prop);
    delete_m(me_tr);
    delete_m(mag_t);
    delete_m(ge_tr);
    return delta_z;
}

static inline struct matrix* 
        meas_mode02(sys_param *sys, state_space *prop,struct matrix* acce ,struct matrix *ge)
{
    int Hr = 3, Hc = 15;
    int Rr = 3, Rc = 3;

    //matriz de saida - delta_z
    struct matrix *delta_z = matrix_m(3, 1);
    //Alloc aux matrices
    struct matrix* rot    = matrix_m(3, 3);
    struct matrix* ge_t   = matrix_m(3, 1); 
    struct matrix* skew_ge= matrix_m(3, 3);
    struct matrix* auxM1  = matrix_m(3, 3);/*rot * skew_ge = 3x3*/
    struct matrix* I_3x3  = matrix_m(3, 3);eye_m(1, I_3x3);
    struct matrix* z_prop = matrix_m(3, 1);
    struct matrix* acc_tr = matrix_m(3, 1); // transposta de prop->acc
    
    struct matrix* H = matrix_m(Hr,Hc);
    struct matrix* R = matrix_m(Rr,Rc); 
    zeros(H);
    zeros(R);

    quat2dcm(prop->qtn, rot);
    crossM_f(ge, skew_ge);
    times_m(rot, skew_ge, auxM1);
    timesc_m(-1.0, auxM1, auxM1);


    //preenchimento de H
    getpart_mat(auxM1, 0, 0, H); //H(1:3,   1:3) = -rot * skew_ge
    getpart_mat(I_3x3, 0, 9, H); //H(1:3, 10:12) = I
    

    //preench. R
    mat_cpy(acce, R);
    transp_m(ge, ge_t);
    times_m(rot, ge_t, z_prop);
    

    //calculo de dz
    transp_m(prop->acc, acc_tr);
    // printf("time_m: ");
    
    sum_m(acc_tr, z_prop, delta_z); //...dz = z - (-rot*ge) == z + rot*ge
    // debug(z, 1000);
    
    sys->R = R;
    sys->H = H;

    //Free memory
    delete_m(rot); 
    delete_m(ge_t); 
    delete_m(skew_ge);
    delete_m(auxM1); 
    delete_m(I_3x3); 
    delete_m(z_prop);
    delete_m(acc_tr); 
    return delta_z;
}

static inline struct matrix* 
        meas_mode03(sys_param *sys, state_space *prop,struct matrix *magn,
                            struct matrix *acce, struct matrix *mag, struct matrix *ge, struct matrix *me){
    //matriz de saida - delta_z
    int Hr = 6, Hc = 15;
    int Rr = 6, Rc = 6;

    struct matrix *delta_z = matrix_m(6, 1);
    struct matrix *rot    = matrix_m(3, 3);
    struct matrix *skew_me= matrix_m(3, 3);
    struct matrix *skew_ge= matrix_m(3, 3);
    struct matrix *skew_me_aux= matrix_m(3, 3);
    struct matrix *skew_ge_aux= matrix_m(3, 3);
    struct matrix *I_3x3  = matrix_m(3, 3);eye_m(1, I_3x3);
    struct matrix *z_prop = matrix_m(6, 1); 
    struct matrix *z = matrix_m(6, 1); 
    struct matrix *auxV1  = matrix_m(3, 1); /*==  rot * me*/
    struct matrix *auxV2  = matrix_m(3, 1); /*== -rot * ge*/
    struct matrix* mag_t = matrix_m(3, 1);
    struct matrix* acc_t = matrix_m(3, 1);
    struct matrix* ge_t  = matrix_m(3, 1);
    struct matrix* me_t  = matrix_m(3, 1);
    
    struct matrix* H = matrix_m(Hr, Hc);
    struct matrix* R = matrix_m(Rr, Rc);
    
    // skew_me = -rot x crossM_f(me) 
    quat2dcm(prop->qtn, rot);
    crossM_f(me,skew_me_aux); 
    times_m(rot, skew_me_aux, skew_me);
    
    // skew_ge = -rot x crossM_f(ge)
    crossM_f(ge,skew_ge_aux); 
    times_m(rot, skew_ge_aux, skew_ge);
    timesc_m(-1, skew_ge, skew_ge);
        
    //Inclui matrizes em H
    zeros(H);
    getpart_mat(skew_me, 0, 0, H);
    getpart_mat(skew_ge, 3, 0, H);
    getpart_mat(I_3x3 ,  3, 9, H);

    //Escrita na matriz R
    zeros(R);
    getpart_mat(magn, 0, 0, R);
    getpart_mat(acce, 3, 3, R);//escrito na diagonal

    // Calculo Z
    transp_m(mag, mag_t);
    transp_m(prop->acc, acc_t);
    
    getpart_mat(mag_t, 0, 0, z);
    getpart_mat(acc_t, 3, 0, z);//copia a coluna em z = [mag ; prop.acc]

    // rot*me
    transp_m(me, me_t);
    times_m(rot, me_t, auxV1);

    // -rot*ge
    transp_m(ge, ge_t);
    times_m(rot, ge_t , auxV2);
    timesc_m(-1, auxV2, auxV2);

    getpart_mat(auxV1, 0, 0, z_prop);
    getpart_mat(auxV2, 3, 0, z_prop);
    less_m(z, z_prop, delta_z);
    
    sys->H = H; 
    sys->R = R;
    
    //Free memory
    delete_m( rot ); 
    delete_m( skew_me ); 
    delete_m( skew_ge ); 
    delete_m( I_3x3 ); 
    delete_m( z_prop ); 
    delete_m( auxV1 ); 
    delete_m( auxV2 ); 
    delete_m(mag_t);
    delete_m(acc_t);
    delete_m(ge_t);
    delete_m(me_t);
    
    return delta_z;
}

static inline struct matrix* 
        meas_mode04(sys_param *sys, state_space *prop,
    struct matrix* vel_gps, struct matrix* pos_gps, struct matrix* gps_position){

    int Hr = 6, Hc = 15;
    int Rr = 6, Rc = 6;

    //Alloc aux matrices
    struct matrix *I_3x3  = matrix_m(3, 3);eye_m(1.0, I_3x3);
    struct matrix *z      = matrix_m(6, 1);
    struct matrix *z_prop = matrix_m(6, 1);
    struct matrix *gps_position_t= matrix_m(3,1);
    struct matrix *prop_vel_t = matrix_m(3,1);
    struct matrix *prop_pos_t= matrix_m(3,1);
    //matriz de saida - delta_z
    struct matrix *delta_z = matrix_m(6, 1);
    struct matrix *H = matrix_m(Hr, Hc);
    struct matrix *R = matrix_m(Rr, Rc);

    //montagem H
    zeros(H);
    getpart_mat(I_3x3, 0,  6, H); /*H(1:3,7:9)   = eye(3);*/
    getpart_mat(I_3x3, 3, 12, H); /*H(4:6,13:15) = eye(3);*/

    //montagem R
    zeros(R);
    getpart_mat(vel_gps, 0, 0, R);
    getpart_mat(pos_gps, 3, 3, R);

    //Montagem z
    zeros(z);
    transp_m(gps_position, gps_position_t);
    getpart_mat(gps_position_t, 3, 0, z);
    
    //Montagem z_prop
    transp_m(prop->vel, prop_vel_t);
    transp_m(prop->pos, prop_pos_t);

    getpart_mat(prop_vel_t , 0, 0, z_prop);
    getpart_mat(prop_pos_t , 3, 0, z_prop);
    
    // debug(prop_pos_t,1);
    // debug(prop_vel_t,1);
    // debug(z_prop, 1);
    // debug(z,1);
    
    less_m(z, z_prop, delta_z);
    // debug(delta_z,1);
    sys->H = H;
    sys->R = R;


    //Free memory
    delete_m(I_3x3);
    delete_m(z);
    delete_m(z_prop);
    delete_m(gps_position_t);
    delete_m(prop_vel_t);
    delete_m(prop_pos_t);

    return delta_z;
}

static inline struct matrix* 
        meas_mode05(sys_param *sys, state_space *prop, 
                            struct matrix* magn, struct matrix* vel_gps, struct matrix* pos_gps,
                            struct matrix* mag, struct matrix* gps_pos,struct matrix* ge, struct matrix *me){
    //matriz de saida - delta_z
    struct matrix *delta_z = matrix_m(9, 1);
    //Alloc aux matrices
    struct matrix *rot    = matrix_m(3, 3);
    struct matrix *skew_me= matrix_m(3, 3);
    struct matrix *subH1 = matrix_m(3, 3);
    struct matrix *I_3x3  = matrix_m(3, 3);
    struct matrix *z      = matrix_m(9, 1);
    struct matrix *z_prop = matrix_m(9, 1);
    struct matrix *auxV1  = matrix_m(3, 1);//auxV1 = rot(3x3) x me(3x1);
    struct matrix *mag_t  = matrix_m(3, 1);
    struct matrix *gps_pos_t= matrix_m(3, 1);
    struct matrix *prop_vel_t=matrix_m(3,1);
    struct matrix *prop_pos_t=matrix_m(3,1);
    struct matrix *me_t=matrix_m(3,1);

    struct matrix* H = matrix_m(9, 15);
    struct matrix* R = matrix_m(9,  9);
    eye_m(1, I_3x3);

    quat2dcm(prop->qtn, rot);
    crossM_f(me, skew_me); 

    //montagem matriz H
    zeros(H);
    times_m(rot, skew_me, subH1);
    getpart_mat(subH1,  0, 0, H); //H(1:3,1:3) = rot*skew_me;
    getpart_mat(I_3x3,  3, 6, H); //H(4:6,7:9) = eye(3);
    getpart_mat(I_3x3,  6,12, H); //H(7:9,13:15) = eye(3);
    
    //montagem matriz R
    zeros(R);
    getpart_mat(magn   , 0, 0, R);
    getpart_mat(vel_gps, 3, 3, R); 
    getpart_mat(pos_gps, 6, 6, R); 

    //montagem de z
    zeros(z);
    transp_m(mag, mag_t);
    transp_m(gps_pos, gps_pos_t);
    getpart_mat(mag_t    , 0, 0, z);
    getpart_mat(gps_pos_t, 6, 0, z);
    
    //montagem z_prop
    transp_m(me, me_t);
    transp_m(prop->vel, prop_vel_t); //prop.vel
    transp_m(prop->pos, prop_pos_t); //prop.pos
    times_m(rot, me_t, auxV1);
    getpart_mat(auxV1     , 0, 0, z_prop);
    getpart_mat(prop_vel_t, 3, 0, z_prop);
    getpart_mat(prop_pos_t, 6, 0, z_prop);
    less_m(z, z_prop, delta_z);
    sys->H = H;
    sys->R = R;
    //Free memory
    delete_m(rot);
    delete_m(subH1);
    delete_m(skew_me);
    delete_m(I_3x3);
    delete_m(z);
    delete_m(z_prop);
    delete_m(auxV1);
    // exit(0);
    return delta_z;
}

static inline struct matrix* 
        meas_mode06(sys_param *sys, state_space *prop, struct matrix *acce, struct matrix* gps_pos,
                               struct matrix *vel_gps,  struct matrix *pos_gps, struct matrix *ge){
    //matriz de saida - delta_z
    struct matrix *delta_z = matrix_m(9, 1);
    //Alloc aux matrices
    struct matrix *rot      = matrix_m(3, 3);
    struct matrix *rot_ge   = matrix_m(3, 1);
    struct matrix *skew_ge  = matrix_m(3, 3);
    struct matrix *I_3x3    = matrix_m(3, 3);eye_m(1.0f, I_3x3);
    struct matrix *z        = matrix_m(9, 1);
    struct matrix *z_prop   = matrix_m(9, 1);
    struct matrix *auxM1    = matrix_m(3, 3);/*-rot x skew_ge */
    struct matrix *ge_t     = matrix_m(3, 1);
    struct matrix *acc_t    = matrix_m(3, 1);
    struct matrix *gps_pos_t= matrix_m(3, 1);
    struct matrix *prop_vt  = matrix_m(3, 1);
    struct matrix *prop_pt  = matrix_m(3, 1);
    struct matrix *H = matrix_m(9, 15);
    struct matrix *R = matrix_m(9 ,9);


    quat2dcm(prop->qtn, rot);
    crossM_f(ge, skew_ge);
    transp_m(ge, ge_t);
    times_m(rot, skew_ge, auxM1);
    timesc_m(-1, auxM1, auxM1);

    // Montagem H
    getpart_mat(auxM1, 0,  0, H); /*H(1:3,1:3) = -rot*skew_ge;*/
    getpart_mat(I_3x3, 0,  9, H); /*H(1:3,10:12) = eye(3);*/
    getpart_mat(I_3x3, 3,  6, H); /*H(4:6,  7:9) = eye(3);*/
    getpart_mat(I_3x3, 6, 12, H); /*H(7:9,13:15) = eye(3);*/

    //Montagem R
    getpart_mat(acce   , 0, 0, R); // montagem de R
    getpart_mat(vel_gps, 3, 3, R);
    getpart_mat(pos_gps, 6, 6, R);

    //Montagem Z
    transp_m(prop->acc, acc_t);  
    transp_m(gps_pos, gps_pos_t);  
    getpart_mat(acc_t     , 0, 0, z); // montagem de z
    getpart_mat(gps_pos_t , 6, 0, z);
    
    //Montagem z_prop
    times_m(rot, ge_t, rot_ge);
    timesc_m(-1.0, rot_ge, rot_ge);
    transp_m(prop->vel, prop_vt);
    transp_m(prop->pos, prop_pt);
    getpart_mat(rot_ge , 0, 0, z_prop);
    getpart_mat(prop_vt, 3, 0, z_prop);
    getpart_mat(prop_pt, 6, 0, z_prop);
    less_m(z, z_prop, delta_z);
    
    sys->H = H;
    sys->R = R;

    // debug(sys->H, 1000);
    // debug(H, 1000);

    // debug(sys->R, 1000);
    // debug(R, 1000);
    
    //Free memory
    delete_m(rot_ge);
    delete_m(rot);
    delete_m(skew_ge);
    delete_m(I_3x3);
    delete_m(z);
    delete_m(z_prop);
    delete_m(auxM1);
    delete_m(ge_t);
    delete_m(acc_t);
    delete_m(gps_pos_t);
    delete_m(prop_vt);
    delete_m(prop_pt);


    return delta_z;
}

static inline struct matrix* 
        meas_mode07(sys_param *sys, state_space *prop, struct matrix* gps_pos, struct matrix *mag,
                    struct matrix* magn, struct matrix* acce, struct matrix* vel_gps, struct matrix* pos_gps, struct matrix* ge, struct matrix* me ){
    //matriz de saida - delta_z
    struct matrix *delta_z = matrix_m(12, 1);
    struct matrix *delta_z_t= matrix_m(1, 12);
    
    //Alloc aux matrices
    struct matrix*rot = matrix_m(3, 3); 
    struct matrix*skew_me = matrix_m(3, 3); 
    struct matrix*skew_ge = matrix_m(3, 3); 
    struct matrix*auxM1  = matrix_m(3, 3); //  rot * skew_me
    struct matrix*auxM2  = matrix_m(3, 3); // -rot * skew_ge
    struct matrix*auxM3  = matrix_m(3, 1); //  rot * me
    struct matrix*auxM4  = matrix_m(3, 1); // -rot * ge
    struct matrix*auxM3_t= matrix_m(1, 3); // ( rot * me)t
    struct matrix*auxM4_t= matrix_m(1, 3); // (-rot * ge)t
    struct matrix*I_3x3  = matrix_m(3, 3); eye_m(1.0, I_3x3);
    struct matrix*z_prop = matrix_m(1, 12); 
    struct matrix*z      = matrix_m(1, 12); 
    struct matrix*me_t  = matrix_m(3, 1); 
    struct matrix*ge_t  = matrix_m(3, 1); 
    struct matrix*prop_acc_t= matrix_m(3, 1); 
    
    struct matrix* H = matrix_m(12, 15); 
    struct matrix* R = matrix_m(12, 12);

    quat2dcm(prop->qtn, rot);
    crossM_f(me, skew_me);
    crossM_f(ge, skew_ge);

    //rot * skew_me
    times_m(rot, skew_me, auxM1);

    // -rot * skew_ge
    times_m(rot, skew_ge, auxM2);
    timesc_m(-1.0, auxM2, auxM2);

    //rot * me
    transp_m(me, me_t);
    times_m(rot, me_t, auxM3);
    
    // - rot*ge
    transp_m(ge, ge_t);
    times_m(rot, ge_t, auxM4);
    timesc_m(-1.0, auxM4, auxM4);

    //montagem H
    zeros(H);
    getpart_mat(auxM1 , 0, 0, H); //H(1:3,1:3) = rot*skew_me;
    getpart_mat(auxM2 , 3, 0, H); //H(4:6,1:3) = -rot*skew_ge;
    getpart_mat(I_3x3 , 3, 9, H); //H(4:6,10:12) = eye(3);
    getpart_mat(I_3x3 , 6, 6, H); //H(7:9,7:9) = eye(3);
    getpart_mat(I_3x3 , 9, 12,H); //H(10:12,13:15) = eye(3);

    // montagem R
    zeros(R);
    getpart_mat(magn   , 0, 0, R); 
    getpart_mat(acce   , 3, 3, R);
    getpart_mat(vel_gps, 6, 6, R);
    getpart_mat(pos_gps, 9, 9, R);
    
    //montagem z 
    zeros(z);
    getpart_mat(mag      , 0, 0, z); 
    getpart_mat(prop->acc, 0, 3, z);
    getpart_mat(gps_pos  , 0, 9, z);

    //montagem z_prop
    transp_m(auxM3, auxM3_t);
    transp_m(auxM4, auxM4_t);
    getpart_mat(auxM3_t   , 0, 0, z_prop); // rot*me
    getpart_mat(auxM4_t   , 0, 3, z_prop); //-rot*ge
    getpart_mat(prop->vel , 0, 6, z_prop);
    getpart_mat(prop->pos , 0, 9, z_prop);

    less_m(z,z_prop, delta_z_t);
    transp_m(delta_z_t, delta_z);

    sys->H = H;
    sys->R = R;


    //Free memory
    delete_m( rot );
    delete_m( skew_me );
    delete_m( skew_ge );
    delete_m( auxM1 );
    delete_m( auxM2 );
    delete_m( auxM3 );
    delete_m( auxM4 );
    delete_m( I_3x3 );
    delete_m( z_prop );
    delete_m( z );
    delete_m(me_t);
    delete_m(ge_t);
    delete_m( prop_acc_t);
    delete_m( auxM3_t);
    delete_m( auxM4_t);
    delete_m( delta_z_t);


    return delta_z;
}


//Calculo do vetor(coluna) de medidas, sendo seu tamanho variavel. Faz chamada de acordo com o valor "mode"
struct matrix*  measurement(sys_param *sys, int mode,state_space *prop,struct matrix* mag,
                 struct matrix* gps_position,struct matrix* ge,
                 struct matrix* me){ 
    
    
    //  Variância da aceleração linear via acelerômetro (calculada de 0 a 10 segundos)
    struct matrix* acce = matrix_m(3,3);
    // % Variância da intensidade do campo magnético via magnetômetro
    struct matrix* magn = matrix_m(3,3);
    // % Variância da posição linear XY via GPS
    struct matrix* pos_gps = matrix_m(3,3);
    // % Variância da velocidade linear XY via GPS
    struct matrix* vel_gps = matrix_m(3,3);
    
    
    struct matrix* delta_z=NULL;//#TODO : como ldar com H (ou delta_z abstrato como 0 ?)
    // #todo: remover esses zeros ?
    zeros(acce);
    zeros(magn);
    zeros(pos_gps);
    zeros(vel_gps);
    // Valores inseridos na diagonal
    acce->elements[0] = acce_x;acce->elements[4] = acce_y;acce->elements[8] = acce_z;
    magn->elements[0] = magn_x;magn->elements[4] = magn_y;magn->elements[8] = magn_z;
    pos_gps->elements[0] = pos_gps_x;pos_gps->elements[4] = pos_gps_y;pos_gps->elements[8] = pos_gps_z;
    vel_gps->elements[0] = vel_gps_x;vel_gps->elements[4] = vel_gps_y;vel_gps->elements[8] = vel_gps_z;
    //todos os casos modificam o tamanho de matriz

    delete_m(sys->H);
    delete_m(sys->R);


    switch (mode){
    case 0:{ //% ---/---/---/---   
        delta_z = meas_mode00(sys, prop);// OK
        break;
    }
    case 1:{//% ---/---/---/mag
        delta_z = meas_mode01(sys, prop, me, magn, mag);
        break;
    }
    case 2:{//% ---/---/acc/---
        delta_z = meas_mode02(sys, prop, acce, ge); //OK
        break;
    }
    case 3:{//% ---/---/acc/mag
        delta_z = meas_mode03(sys, prop, magn, acce, mag, ge, me); // OK
        break;
    }
    case 4:{ // ---/gps/---/---
        delta_z = meas_mode04(sys,prop, vel_gps, pos_gps, gps_position); // #CHECK
        break; 
    }
    case 5:{//% ---/gps/---/mag   
        delta_z = meas_mode05(sys, prop, magn, vel_gps, pos_gps, mag, gps_position, ge, me);
        break; 
    }
    case 6:{ //% ---/gps/acc/---
        delta_z = meas_mode06(sys,prop, acce, gps_position, vel_gps, pos_gps, ge); 
        break; 
    }
    case 7:{ //% ---/gps/acc/mag
        delta_z = meas_mode07(sys, prop, gps_position, mag, magn, acce, vel_gps, pos_gps, ge,me);
        break; 
    }
    default: bound_error();
        break;
    }

    //considerando que nao houve falha de borda (bound)
    // cont_mode->elements[mode]+=1;
    
    // debug(sys->H, 1000);
    // debug(sys->R, 1000);
    // debug(delta_z,1000);


    delete_m(acce);
    delete_m(magn);
    delete_m(pos_gps);
    delete_m(vel_gps);


    return delta_z;
}


struct matrix* kalman_filter(int mode, state_space *prop, struct matrix* delta_z, sys_param *sys, int nx, struct matrix *P)
{

    //verifica estruturas necessarias
    if(!prop || !sys || !P) alloc_error();
    int nrP = prop->P->n_row;
    int ncP = prop->P->n_row;
    struct matrix *delta_x = matrix_m(nx, 1);

    //dimensoes de matriz
    int nrH = sys->H->n_row;
    int ncH = sys->H->n_col;
    struct matrix *Ht     = matrix_m(ncH, nrH); /*PHt = P * H'*/
    struct matrix *PHt    = matrix_m(nrP, nrH); /*PHt = P * H'*/
    struct matrix *HPHt   = matrix_m(nrH, nrH); /*HpHt= H * (P * H') */
    struct matrix *K      = matrix_m(nrP, nrH); 
    struct matrix *KH     = matrix_m(nrP, nrP);
    struct matrix *KH_aux = matrix_m(nrP, nrP);
    struct matrix *I_nx   = matrix_m(nx,  nx);
    struct matrix *aux_inv= matrix_m(nrH, nrH);  //matriz quadrada (A , At) => [A(row) x A(row)]
    
    eye_m(1.0, I_nx);
    // struct matrix *aux_P = matrix_m(nrP, ncP); 
    zeros(delta_x);
    

    // Ganho de Kalman
    if(mode == 0){
        zeros(K);
    } else {
        // K = (P x H') x inv(H x (P x H') + R);
        transp_m(sys->H, Ht);
        times_m(prop->P, Ht, PHt);
        times_m(sys->H, PHt, HPHt);     //    H x (P x H')
        sum_m(HPHt, sys->R, HPHt);
        inv_m(HPHt,  aux_inv);          //    H x (P x H') + R
        times_m(PHt,aux_inv, K);        //inv(H x (P x H') + R)
        times_m(K, delta_z, delta_x); 
    }

    //% Correção da matriz de covariância
    // P = (eye(nx) - K * system.H) * prop.P;
    times_m(K, sys->H, KH_aux); // KH = K x H
    less_m(I_nx, KH_aux, KH);// KH = I - KH
    times_m(KH, prop->P, P);


    // debug(prop->P, 1e3);
    // debug(sys->H, 1e3);
    // debug(PHt, 100);
    // debug(HPHt, 1000);
    // debug(aux_inv, 100);
    // debug(P, 100);


    //free memory
    delete_m(Ht);
    delete_m(PHt);
    delete_m(HPHt);
    delete_m(KH_aux);
    delete_m(K);
    delete_m(aux_inv);
    delete_m(KH);
    delete_m(I_nx);

    return delta_x;   
}

void update(state_space *prop,struct matrix* delta_x,int mode,state_space *state){

    // if(delta_x == NULL)
    //     delta_x = matrix_m(15,1);

    if(mode == 0)
        mat_cpy(prop->qtn, state->qtn);
    else{
        struct matrix* delta_x1 = matrix_m(1, 3);
        getpart_vec(delta_x, 0, 2, delta_x1);
        deltaAngle2quat_f(delta_x1, prop->qtn, state->qtn);
        delete_m(delta_x1);
    }
    struct matrix *delta_x2 = matrix_m(1, 3);//delta(4:6)
    struct matrix *delta_x3 = matrix_m(1, 3);//delta(7:9)
    struct matrix *delta_x4 = matrix_m(1, 3);//delta(10:12)
    struct matrix *delta_x5 = matrix_m(1, 3);//delta(13:15)

    //retiro parte do vetor
    
    // printf("mat sz: "); mat_size(delta_x);
    getpart_vec(delta_x, 3,  5, delta_x2);
    getpart_vec(delta_x, 6,  8, delta_x3);
    getpart_vec(delta_x, 9, 11, delta_x4);
    getpart_vec(delta_x, 12,14, delta_x5);

    // % Atualiza o bias da velocidade angular
    sum_m(delta_x2 ,prop->bias_gyr, state->bias_gyr); 

    // % Atualiza a valocidade linear
    sum_m(delta_x3 ,prop->vel, state->vel);
    
    // % Atualiza o bias da aceleração linear
    sum_m(delta_x4 ,prop->bias_acc, state->bias_acc);
    // debug(delta_x, 1000);
    // debug(delta_x4, 1000); 
    // debug(prop->bias_acc, 1000);
    
    sum_m(delta_x5 ,prop->pos, state->pos);

    delete_m( delta_x2 );
    delete_m( delta_x3 );
    delete_m( delta_x4 );
    delete_m( delta_x5 );
}