#include "file_reader.h"

#define MAX_LEN (30000)


dataset* load_dataset()
{
    dataset *out;
    out = (dataset*)malloc(sizeof(dataset));
    if(!out) alloc_error();
    
    
    out->gps = import_GPS2();///verificar
    out->imu = import_IMU();
    out->mag = import_MAG(); 
    out->rpy = import_RPY();
    

        // Ajuste de t0
    double t0 = out->imu->tim->elements[0];

    int size1 = out->mag->tim->n_row;
    int size2 = out->gps->tim->n_row;
    
    for(int i = 0; i < size1; i++)
         out->mag->tim->elements[i] =out->mag->tim->elements[i] -  t0;

    for(int i = 0; i < size2; i++)
         out->gps->tim->elements[i] =out->gps->tim->elements[i] -  t0;
     


    return out;
}


//stores a csv file in matrix structure
struct matrix* read_file_nD(char name[], int nD){
    FILE *file;
    struct matrix* M;
    int ctr = 0; // limiter counter
    
    if ((file = fopen(name, "r")) == NULL) file_error();
    
    double data[MAX_LEN][nD];
    while (ctr < MAX_LEN && fscanf(file, "%lf", &data[ctr][0]) == 1) {
        for (int i = 1; i < nD; i++) {
            int a = fscanf(file, ",%lf", &data[ctr][i]);
        }
        ctr++;
    }
    fclose(file);// Close file
    //alloc in matrix
    printf("allocated: %d x %d\n", ctr, nD);
    M = matrix_m(ctr, nD);
    if(M == NULL) alloc_error();
    
    double *pM = M->elements;
    int idx = 0;
    for (int i = 0; i < ctr; i++) {
        for (int j = 0; j < nD; j++) {
            pM[idx] = data[i][j];
            idx++;
        }
    }

    return M;
}



IMU_data *import_IMU()
{
    IMU_data* data;
    data = (IMU_data*) malloc(sizeof(IMU_data));
    
    if(data == NULL) alloc_error();
    
    char fileTim[] = "DadosZurich/imu_tim.csv";
    char fileGyr[] = "DadosZurich/imu_gyr.csv";
    char fileAcc[] = "DadosZurich/imu_acc.csv";

    data->tim = read_matrix(fileTim);
    
    data->acc = read_matrix(fileAcc);
    
    data->gyr = read_matrix(fileGyr);

    return data;
}
GPS_data *import_GPS2()
{
    GPS_data* data;
    data = (GPS_data*) malloc(sizeof(GPS_data));

    if(data == NULL) alloc_error();
    
    char fileTim[] = "DadosZurich/gps_tim.csv";
    char filePos[] = "DadosZurich/gps_pos.csv";

    data->tim = read_matrix(fileTim);
    data->pos = read_matrix(filePos);

    return data;
}

GPS_data *import_GPS()
{
    GPS_data* data;
    data = (GPS_data*) malloc(sizeof(GPS_data));

    if(data == NULL) alloc_error();
    
    char fileLat[] = "DadosZurich/gps_lat.csv";
    char fileTim[] = "DadosZurich/gps_tim.csv";
    char fileLon[] = "DadosZurich/gps_lon.csv";
    char fileHgt[] = "DadosZurich/gps_hgt.csv";

    
    struct matrix *aux_lat, *aux_lon, *aux_hgt;
    struct matrix *aux_llh0= matrix_m(3,1);
    struct matrix *aux_llh = matrix_m(3,1);
    struct matrix *aux_ned = matrix_m(3,1);

    //Alocacao e atribuição de dados
    data->tim = read_matrix(fileTim);
    aux_lat   = read_matrix(fileLat);
    aux_lon   = read_matrix(fileLon);    
    aux_hgt   = read_matrix(fileHgt);
    int data_size = aux_lat->n_row; 

    data->pos = matrix_m(data_size, 3);

        
    aux_llh0->elements[0] = aux_lat->elements[0];
    aux_llh0->elements[1] = aux_lon->elements[0];
    aux_llh0->elements[2] = aux_hgt->elements[0];
    double t0 = data->tim->elements[0];

    //#TODO: possivel melhoria no desempenho(memoria), sendo uso de lat,lon.. ?
    for(int i=0; i < data_size;i++){
        aux_llh->elements[0] = aux_lat->elements[i];
        aux_llh->elements[1] = aux_lon->elements[i];
        aux_llh->elements[2] = aux_hgt->elements[i];
        geodetic2ned(aux_llh, aux_llh0, aux_ned);
    
        data->pos->elements[3*i]   = aux_ned->elements[0];        
        data->pos->elements[3*i+1] = aux_ned->elements[1];
        //mudança referencial
        data->pos->elements[3*i+2] = -aux_ned->elements[2];
         data->tim->elements[i] -= t0;
    }
    //free memory
    delete_m(aux_llh0);
    delete_m(aux_llh);
    delete_m(aux_ned);
    delete_m(aux_lat);
    delete_m(aux_lon);
    delete_m(aux_hgt);
    return data;
}
RPY_data *import_RPY()
{
    RPY_data* data;
    data = (RPY_data*) malloc(sizeof(RPY_data));
    if(data ==NULL) alloc_error();


    char fileTim[] = "DadosZurich/rpy_tim.csv";
    char fileRPY[] = "DadosZurich/rpy_rpy.csv";

    data->tim = read_matrix(fileTim);
    data->rpy = read_matrix(fileRPY);
    return data;

}

MAG_data *import_MAG()
{
    MAG_data* data;
    data = (MAG_data*) malloc(sizeof(MAG_data));
    
    if(data ==NULL) alloc_error();
    
    
    char fileMag[] ="DadosZurich/mag_mag.csv";
    char fileTim[] ="DadosZurich/mag_tim.csv";    
    
    data->mag = read_matrix(fileMag);
    data->tim = read_matrix(fileTim);

    return data;
}
// ================================
void freeIMU(IMU_data *data)
{   
    if(data ==  NULL) return;
    delete_m(data->tim);
    delete_m(data->acc);
    delete_m(data->gyr);
    free(data);
}

void freeGPS(GPS_data *data)
{
    if(data ==  NULL) return;
    delete_m(data->tim);
    delete_m(data->pos);
    free(data);
}

void freeRPY(RPY_data *data)
{   
    if(data ==  NULL) return;
    delete_m(data->tim);
    delete_m(data->rpy);
    free(data);
}

void freeMAG(MAG_data *data)
{
    if(data ==  NULL) return;
    delete_m(data->tim);
    delete_m(data->mag);
    free(data);
}

double local_gravity(IMU_data *imu_data, double initial_time){
    
    double x2,y2,z2;
    double gravity = 0.0f;
    int count=0, j=1, nD=3;
    double *p_imu = imu_data->acc->elements; 
    double *p_time= imu_data->tim->elements;

    while(p_time[j]<= initial_time){
        j+=1;
        x2 = p_imu[j*nD+0] * p_imu[j*nD+0];
        y2 = p_imu[j*nD+1] * p_imu[j*nD+1];
        z2 = p_imu[j*nD+2] * p_imu[j*nD+2];
        gravity += sqrt((x2+y2+z2));
        count+=1;
    }
    
    return gravity / count;
}

double local_magnetic(MAG_data *mag_data,double initial_time){
    double x2,y2,z2;
    double magnetic = 0.0f;
    int count=0, j=1, nD=3;
    double *p_mag = mag_data->mag->elements; 
    double *p_time= mag_data->tim->elements;

    while(p_time[j]<= initial_time){
        j+=1;
        x2 = p_mag[j*nD+0] * p_mag[j*nD+0];
        y2 = p_mag[j*nD+1] * p_mag[j*nD+1];
        z2 = p_mag[j*nD+2] * p_mag[j*nD+2];
        magnetic += sqrt((x2+y2+z2));
        count+=1;
    }
    
    return magnetic / count;
}


