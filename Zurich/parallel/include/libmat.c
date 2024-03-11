// Libray of matrices


#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "libmat.h"
#include "defines.h"
#include "err_handler.h"
#include <stdbool.h>

// Create a matrix
struct matrix *matrix_m(int n_row, int n_col) {
        struct matrix *M = (struct matrix *) malloc(sizeof(struct matrix));
        if(!M) alloc_error();
        
        M->elements = (__data_type*) (__data_type*) calloc(n_row * n_col, sizeof(__data_type));
        if(M->elements == NULL) alloc_error();
        
        M->n_row = n_row;
        M->n_col = n_col;
        return M;
}
// Delete a matrix
void delete_m(struct matrix *M) {
        if(M == NULL) return;
        
        free(M->elements);
        free(M);
}

struct matrix* rnd_matrix(int row, int col)
{
        struct matrix *A= matrix_m(row, col);
        
        for(int i=0; i < row; i++)
        for(int j=0; j < col; j++)
                A->elements[i* col + j] = (__data_type) (rand()%100);

        return A;
}
// Print a matrix
void print_m(struct matrix *M , __data_type a) {
        __data_type *p;
        int nr, nc;
        int i,j;

        p =  M->elements;
        nr = M->n_row;
        nc = M->n_col;
        printf("[\n");
        for (i=0; i<nr; i++)  {
            // if(i>10) printf("%d", i);
            for (j=0; j<nc; j++)  {
                IsEqual(a * p[i*nc + j], 0.0) ? printf("  %8i", 0):
                			printf("  %08.4lf", a * p[i*nc + j]); 
                                                    
            }
             printf(" ;\n");               
        }
        printf("\n]");
}
        


// Eye of a matrix
void eye_m(__data_type value, struct matrix *M) {

        __data_type *p;
        int nr, nc;
        int i,j;

        p = M->elements;
        nr = M->n_row;
        nc = M->n_col;

        for (i=0; i<nr; i++) {
                for (j=0; j<nc; j++) {
                        if (i == j) {
                                p[i*nc + j]=value;
                        }
                        else {
                                p[i*nc + j]=0;
                        }
                }
        }
}

// Equality of matrix
short int equal_m(struct matrix *A, struct matrix *B) {
        int nrA, ncA, nrB, ncB;
        int i,j;
        __data_type *pA, *pB;

        pA = A->elements;
        pB = B->elements;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;

        if (nrA != nrB | ncA != ncB) {
                printf("\n Wrong dimension \n") ;
                return 0;
        }
        else {
                for (i = 0; i< nrA; i++) {
                        for (j = 0; j< ncA; j++) {
                                pA[i*ncA + j] = pB [i*ncB + j];
                        }
                }
                return 1;
        }
}

void equalpart_m(struct matrix *A, int i_0, int j_0, struct matrix *B) {
        int nrA, ncA, nrB, ncB;
        int nr_f, nc_f;
        int i_A,j_A, i_B, j_B;
        __data_type *pA, *pB;



        pA = A->elements;
        pB = B->elements;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;

        nr_f = (i_0 + nrB);  // final row
        nc_f = (j_0 + ncB);  // final column

        for (i_A = i_0, i_B =0; i_A < nr_f; i_A++, i_B++) {
                for (j_A = j_0, j_B = 0; j_A < nc_f; j_A++, j_B++) {
                        pA[i_A*ncA + j_A] = pB [i_B*ncB + j_B];
                }
        }

}

void getpart_m(struct matrix *A, int i_0, int j_0, struct matrix *B) {

        int nrA, ncA, nrB, ncB;
        int nr_f, nc_f;
        int i_A,j_A, i_B, j_B;
        __data_type *pA, *pB;

        pA = A->elements;
        pB = B->elements;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;

        nr_f = (i_0 + nrB);  // final row
        nc_f = (j_0 + ncB);  // final column

        for (i_A = i_0, i_B =0; i_A < nr_f; i_A++, i_B++) {
                for (j_A = j_0, j_B = 0; j_A < nc_f; j_A++, j_B++) {
                         pB [i_B*ncB + j_B] = pA[i_A*ncA + j_A];
                }
        }

}


// Sum of matrices
short int sum_m(struct matrix *A, struct matrix *B, struct matrix *result) {

        __data_type *pA, *pB, *pR;
        int nrA, nrB, nrresult, ncA, ncB, ncresult;
        int cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pB = B->elements;
        pR = result->elements;

        if ((nrA != nrB) | (ncA != ncB) | (nrresult != nrA) | (ncresult != ncA)) {
                printf("Wrong dimension");
                return 0;
        }
        #pragma simd vectorlength(8)
        for(cont = 0; cont<(nrA*ncA); cont++) {
                pR[cont] = pA[cont] + pB[cont];
        }
        return 1;


}

// Less of matrix
short int less_m(struct matrix *A, struct matrix *B, struct matrix *result) {

        __data_type *pA, *pB, *pR;
        int nrA, nrB, nrresult, ncA, ncB, ncresult;
        int cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pB = B->elements;
        pR = result->elements;

        if ((nrA != nrB) | (ncA != ncB) | (nrresult != nrA) | (ncresult != ncA)) {
                printf("Wrong dimension");
                return 0;
        }
        #pragma simd vectorlength(8)
        for(cont = 0; cont<(nrA*ncA); cont++) {
                pR[cont] = pA[cont] - pB[cont];
        }
        
        return 1;
}

// Times between two matrices
short int times_m(struct matrix *A, struct matrix *B, struct matrix *result) {

        __data_type *pA, *pB, *pR, aux;
        int nrA, nrB, nrresult, ncA, ncB, ncresult;
        int i, j, cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrB = B->n_row;
        ncB = B->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pB = B->elements;
        pR = result->elements;
        __data_type mem_aux[nrA*ncB];
        if ((ncA != nrB) ||  (nrresult != nrA) || (ncresult != ncB)) {
                printf("Wrong dimension");
                return 0;
        }

        for (i = 0; i<nrA; i++) {
                for (j = 0; j<ncB; j++) 
                {
                        aux = 0;
                        #pragma simd vectorlength(AVX_LEN)
                        for (cont=0; cont<nrB; cont++ ) 
                                aux=aux + pA[ncA*i + cont]*pB[ncB*cont + j];
                        
                        pR[ncB*i + j]= aux;
                }
        }
        return 1;
        
}

// Times between a constant and a matrix
short int timesc_m(__data_type value, struct matrix *A, struct matrix *result) {

        __data_type *pA, *pB, *pR, aux;
        int nrA, nrresult, ncA, ncresult;
        int cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pR = result->elements;

        if ((nrresult != nrA) || (ncresult != ncA)) {
                printf("Wrong dimension");
                return 0;
        }

        #pragma simd vectorlength(AVX_LEN)
        for (cont=0; cont<(nrA*ncA); cont++ ) {
                pR[cont] = value*pA[cont];
        }
        return 1;
}
// Transpose of a matrix
short int transp_m(struct matrix *A, struct matrix *result) {

        __data_type *pA, *pB, *pR;
        int nrA, nrresult, ncA, ncresult;
        int i, j, cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pR = result->elements;

        if ((nrresult != ncA) ||( ncresult != nrA)) {
                printf("Wrong dimension");
                return 0;
        }
        
        for (i=0; i<nrA; i++) {
                for (j=0; j<ncA; j++) {
                        pR[nrA*j+i] = pA[ncA*i+j];
                }
        }
        return 1;
}


// Inverse of a matrix
short int inv_m(struct matrix *A, struct matrix *result) {

        __data_type *pA, *pB, *pR, *pa;
        __data_type sum, aux, *b, *x;
        struct matrix *a;
        int nrA, nrresult, ncA, ncresult;
        int idx2, *idx, mem, flag;
        int i,j,k, cont;
        nrA = A->n_row;
        ncA = A->n_col;
        nrresult = result->n_row;
        ncresult = result->n_col;
        pA = A->elements;
        pR = result->elements;

        a = matrix_m(nrA,ncA);
        pa = a->elements;

       if ((nrresult != nrA) | (ncresult != ncA)) {
                printf("Wrong dimension");
                return 0;
        }
        else {

        for (i = 0; i<nrA; i++) {
                #pragma simd vectorlength(AVX_LEN)
                for (j = 0; j<nrA; j++) {
                        pa[ncA*i+j] = pA[ncA*i+j];
                }
        }
//---------------------------- Partial pivoting --------------------------------
        b = (__data_type*) malloc(nrA*sizeof(__data_type));
        x = (__data_type*) malloc(nrA*sizeof(__data_type));
        idx = (int*) malloc(nrA*sizeof(int));
        for (k = 0; k<nrA; k++)
                idx[k] = k;

        for (i = 0; i<nrA; i++) {
                j = i;
                idx2 = i;
                if (pa[ncA*i+j] == 0) {
                        flag = 1;
                        for (k = i+1; k<nrA; k++ ) {
                                if (fabs(pa[ncA*k+j]) >= TINY && flag == 1) {
                                        mem  = idx[i];
                                        idx[i] = idx[k];
                                        idx[k] = mem;
                                        idx2 = k;
                                        flag = 0;
                                }
                        }
                        if (flag == 1) {
                                for (k = 0; k<nrA; k++) {
                                        if (fabs(pa[ncA*k+j]) > TINY && fabs(pa[ncA*i+k]) > TINY) {
                                                mem = idx[i];
                                                idx[i] = idx[k];
                                                idx[k] = mem;
                                                idx2 = k;
                                                flag = 0;
                                        }
                                }
                        }
                        if (idx2 == i){
                                printf("\n Singular matrix \n \n");
                                pa[ncA*i+j] = TINY;
                        }
                        #pragma simd vectorlength(AVX_LEN)
                        for (k = 0; k<nrA; k++){
                                mem = pa[ncA*i+k];
                                pa[ncA*i+k] = pa[ncA*idx2+k];
                                pa[ncA*idx2+k] = mem;
                        }
                }

        }

//------------------- Crout's algorithm for LU Decomposition -------------------
        for (j = 0; j<nrA; j++) {
                for (i = 0; i<nrA; i++) {
                        if (i<j | i ==j) {
                                sum = pa[ncA*i+j];
                                for (k = 0; k<i; k++) {
                                        sum = sum - pa[ncA*i+k]*pa[ncA*k+j];
                                }
                                pa[ncA*i+j] = sum;
                        }
                        if (i > j) {
                                sum = pa[ncA*i+j];
                                for (k = 0; k<j; k++) {
                                        sum = sum - pa[ncA*i+k]*pa[ncA*k+j];
                                }
                                pa[ncA*i+j] = sum/pa[ncA*j+j];
                        }
                }
        }
//---------------------------- Forward substituion -----------------------------
        for (k = 0; k<nrA; k++) {
                #pragma simd vectorlength(AVX_LEN)
                for (cont = 0; cont<nrA; cont ++ ) {
                        b[cont] = 0;
                }
                b[k] = 1;
                for (i = 0; i<nrA; i++) {
                        sum = b[i];
                        #pragma simd vectorlength(AVX_LEN)
                        for (j = 0; j<i; j++) {
                                sum = sum - pa[ncA*i+j]*x[j];
                        }
                        x[i] = sum;
                }
//---------------------------- Backward substituion ----------------------------
                for (i=(nrA-1); i>=0; i--) {
                        sum = x[i];
                        for (j = i+1; j<nrA; j++) {
                                sum = sum - pa[ncA*i+j]*x[j];
                        }
                        x[i] = sum/pa[ncA*i+i];
                }
                for (cont = 0;  cont<nrA; cont++){
                        pR[ncA*cont+idx[k]] = x[cont];
                }
        }
        delete_m(a);
        free(b);
        free(x);
        free(idx);
        return 1;
        }
}

// Givens of a matrix
void givens_m(struct matrix *A, struct matrix *result) {

        struct matrix *Theta;
        int nrA, ncA, i, j, i2, j2, cont;
        __data_type a, b, c, rho;
        __data_type *pA, *pTheta;
        short int flag;

        nrA = A->n_row;
        ncA = A->n_col;
        pA = A->elements;

        Theta = matrix_m(ncA,ncA);
        pTheta = Theta->elements;

        for (i =0; i<nrA; i++){
                for (j = ncA-1; j>=i+1; j--) {
                        b = pA[ncA*i+j];
                        flag = 0;
                        for (cont = i; cont < j; cont ++) {
                                a = pA[ncA*i+cont];
                                if (fabs(a) >= TINY) {
                                        flag = 1;
                                        break ;
                                }
                        }
                        if (flag ==0) {
                                a = TINY;
                                printf("\n a = 0 \n");
                        }
                        rho = b/a;
                        for (i2 =0; i2<ncA; i2++){
                                for (j2 =0; j2<ncA; j2++){
                                        if (i2 == j2) {
                                                pTheta[ncA*i2+j2] = 1;
                                        }
                                        else {
                                                pTheta[ncA*i2+j2] = 0;
                                        }
                                }
                        }
                        c = 1/sqrt(1 +rho*rho);
                        pTheta[ncA*cont+cont]  = c;
                        pTheta[ncA*cont+j]  = -rho*c;
                        pTheta[ncA*j+cont]  = rho*c;
                        pTheta[ncA*j+j]  = c;

                        times_m(A, Theta, result);
                }
        }
}

void euler2quat_m(struct euler angles, struct matrix *quat) {

        __data_type *q, psi, theta, phi;
        int nrq, ncq;

        psi   =  angles.yaw;
        theta =  angles.pitch;
        phi   =  angles.roll;

        q = quat->elements;

        q[0] = cos(phi/2)*cos(theta/2)*cos(psi/2) + sin(phi/2)*sin(theta/2)*sin(psi/2);
        q[1] = sin(phi/2)*cos(theta/2)*cos(psi/2) - cos(phi/2)*sin(theta/2)*sin(psi/2);
        q[2] = cos(phi/2)*sin(theta/2)*cos(psi/2) + sin(phi/2)*cos(theta/2)*sin(psi/2);
        q[3] = cos(phi/2)*cos(theta/2)*sin(psi/2) - sin(phi/2)*sin(theta/2)*cos(psi/2);
}

short int crossM_f(struct matrix *A, struct matrix *M){
/*
% Descrição: Retorna a matriz antissimétrica M de um vetor da velocidade angular
% Entrada: vetor em R3
% Saída: Matriz antissimétrica de produto
*/
        if(A == NULL || M == NULL) alloc_error();

        if(A->n_col * A->n_row != 3 || M->n_col * M->n_row != 9){
                printf("Error, dimensions differ!\n");
                return 0;
        }
        __data_type v1,v2,v3;

        v1 = A->elements[0];
        v2 = A->elements[1];
        v3 = A->elements[2];

        M->elements[0] =  0; M->elements[1] =-v3; M->elements[2] = v2;
        M->elements[3] = v3; M->elements[4] =  0; M->elements[5] =-v1;
        M->elements[6] =-v2; M->elements[7] = v1; M->elements[8] =  0;

        return 1;
}

//copy a matrix
void mat_cpy(struct matrix *src, struct matrix* dst){

        if(!src) alloc_error();
        
        if(!dst) dst = matrix_m(src->n_row, src->n_col);
        
        if(src->n_col * src->n_row != dst->n_col * dst->n_row)
                length_error();

        __data_type *psrc = src->elements;
        __data_type *pdst = dst->elements; 
        

        int nr = src->n_row;
        int nc = src->n_col;
        int size =  (nr*nc);

        for(int i=0; i < size; i++)
                pdst[i] = psrc[i];
}

void mat_cpyRow(struct matrix *A_src, int r_src ,struct matrix* B_dst, int r_dst){
        
        int ncA = A_src->n_col;
        int nrA = A_src->n_row;
        int ncB = B_dst->n_col;
        int nrB = B_dst->n_row;
        
        if(!B_dst || !A_src) alloc_error(); 

        if(nrA-1 < r_src || nrB-1 < r_dst) bound_error();

        __data_type *p_dst = B_dst->elements;
        __data_type *p_src = A_src->elements;
        
        //Copia pelo total de colunas
        for(int i=0; i< ncA ; i++){
                p_dst[r_dst*ncB + i] = p_src[r_src*ncA + i];
        }
}

//returns the norm of a vector/matrix
__data_type norm_m(struct matrix *v){

        if(!v) alloc_error();
        
        __data_type acc = 0;
        int t_elem = v->n_col * v->n_row;
        for(int i = 0; i < t_elem; i++)
                acc += v->elements[i] * v->elements[i];
    return sqrt(acc);
}

int maxIdx_v(__data_type *A, int size_max){
        
        if(A == NULL) alloc_error();

        //int size = (A->n_col > A->n_row) ? A->n_col : A->n_row; 
        int max_idx = 0;
        __data_type max = A[0];
        
        for (int i = 0; i < size_max; i++) {
                if (max < A[i]) {
                        max = A[i];
                        max_idx = i;
                }
        }
        return max_idx;
}

//C = A + B 
void unsafe_sum(struct matrix *A, struct matrix *B, struct matrix*C){

        if(!A || !B || !C) alloc_error();
        int size = A->n_col*A->n_row;
        for(int i = 0; i < size; i++){
                C->elements[i] = B->elements[i] + A->elements[i];
        }

}

//C = A - B 
void unsafe_less(struct matrix *A, struct matrix *B, struct matrix*C)
{

        if(!A || !B || !C) alloc_error();
        int Asz = A->n_col*A->n_row;
        int Bsz = B->n_col*B->n_row;
        int Csz = C->n_col*C->n_row;

        if(Asz != Bsz || Bsz != Csz) bound_error();
        
        for(int i = 0; i < Asz; i++){
                C->elements[i] = A->elements[i] - B->elements[i];
        }
}


void unsafe_timesc(__data_type v, struct matrix *A, struct matrix *B){
        
        if(!A || !B) alloc_error();
        
        int size = A->n_col * A->n_row;

        for(int i = 0; i < size; i++){
                B->elements[i] = v * A->elements[i];
        }

}


void getpart_mat(struct matrix *src, int i_0, int j_0, struct matrix *dst) {

        int nrA, ncA, nrB, ncB;
        int nr_f, nc_f;
        int i_A,j_A, i_B, j_B;
        __data_type *pA, *pB;

        //matriz que recebe copia
        pA = dst->elements;
        nrA = dst->n_row;
        ncA = dst->n_col;

        //matrix que será copiada
        nrB = src->n_row;
        ncB = src->n_col;
        pB = src->elements;

        //limites de acesso da matriz
        int bound_r = i_0 + nrB;
        int bound_c = j_0 + ncB;

        if(bound_c > ncA || bound_r > nrA) bound_error();

        int g_i, g_j;//local global de acesso a matriz A
        int id_A, id_B;
        for(int i = 0; i < nrB; i++){
                for(int j = 0; j < ncB; j++){
                        g_i = i_0 + i;  g_j = j_0 + j;
                        id_A = g_i*ncA + g_j; id_B = i*ncB + j;
                        pA[id_A] = pB[id_B];



                }

        }

}

void mat_size(struct matrix* A){
    if(!A) printf("Null\n");
    else printf("%d x %d\n", A->n_row, A->n_col);
}

void getpart_vec(struct matrix *src, int init, int end, struct matrix *dst){
        
        
        if(src ==NULL || dst== NULL)
                alloc_error();

        int nrA = src->n_row;
        int ncA = src->n_col;
        int nrB = src->n_row;
        int ncB = src->n_col;


        int max = fmax(nrB,ncB);
        int n_elem =  fmax(nrA,ncA);
        int offset = end- init;
        if (init < 0 || end >= n_elem || offset > max) bound_error();

    for (int i = init, j = 0; i <= end; i++, j++)
        dst->elements[j] = src->elements[i];

}


void geodetic2ecef(struct matrix* llh, struct matrix* xyz){
        __data_type lat = angle2r(llh->elements[0]);
        __data_type lon = angle2r(llh->elements[1]);
        __data_type h = llh->elements[2];

        // World Geodetic System 1984
        // b = 6356752.31424518; SemiminorAxis [m]
        __data_type a = 6378137; // SemimajorAxis [m]
        __data_type e2 = wg84_Eccentricity*wg84_Eccentricity;

        __data_type RN = wg84_SemimajorAxis / (sqrt((1-e2*(sin(lat)*sin(lat)))));

        xyz->elements[0] = (RN + h)*cos(lat)*cos(lon);    // x
        xyz->elements[1] = (RN + h)*cos(lat)*sin(lon);    // y
        xyz->elements[2] = (RN*(1-e2)+h)*sin(lat);        // z
}

void ecef2ned(struct matrix* pe, struct matrix *llh0, struct matrix *pt){
        __data_type lat0= angle2r(llh0->elements[0]);
        __data_type lon0= angle2r(llh0->elements[1]);
        __data_type h0  = llh0->elements[2];

        struct matrix *RTE = matrix_m(3,3);
        struct matrix *pe0 = matrix_m(3,1); 
        struct matrix *aux = matrix_m(3,1);
        __data_type *pRTE  = RTE->elements;

        __data_type slo = sin(lon0);
        __data_type sla = sin(lat0);
        __data_type clo = cos(lon0);
        __data_type cla = cos(lat0);

        pRTE[0] = -sla*clo; pRTE[1] = -sla*slo; pRTE[2] = cla;
        pRTE[3] = -slo    ; pRTE[4] =  clo    ; pRTE[5] = 0;
        pRTE[6] = -cla*clo; pRTE[7] = -cla*slo; pRTE[8] = -sla;

        geodetic2ecef(llh0, pe0);
        less_m(pe, pe0, aux);
        times_m(RTE,aux, pt);
               

}

void geodetic2ned(struct matrix* llh, struct matrix* llh0, struct matrix *ned)
{
        struct matrix *xyz = matrix_m(3,1);
        geodetic2ecef(llh, xyz); 
        ecef2ned(xyz, llh0, ned);
        delete_m(xyz);
}

void save_matrix(struct matrix *M, char *name){
        FILE *file;
        __data_type *p;
        int nr, nc;
        
        p = M->elements; 
        nr = M->n_row;
        nc = M->n_col;

        if((file=fopen(name, "w")) == NULL) file_error();

        fprintf(file, "%d,%d\n", nr,nc);
        for (int i=0; i<nr; i++){
        for (int j=0; j<nc; j++){
                if (j == nc-1){
                        fprintf(file, "%5.17lf", p[i*nc + j]);
                        fprintf(file, "\n");
                }else fprintf(file, "%5.17lf,", (p[i*nc + j]));
                
        }
        }
        fclose(file);
}

struct matrix* read_matrix(char* name)
{
        FILE *file;
        int nr, nc;
        __data_type temp, *p;

        if ((file = fopen(name, "r")) == NULL) file_error();
        int a = fscanf(file, "%d,%d\n", &nr,&nc);

        struct matrix*M = matrix_m(nr, nc);
        if(!M) alloc_error();

        p = M->elements;
        
        for(int i =0; i < nr; i++){
                for(int j=0; j < nc; j++){
                        if (j == nc-1) fscanf(file, "%lf", &temp);
                        else fscanf(file, "%lf,", &temp);
                        
                        p[i*nc+j] = temp;
                }
        }


        fclose(file);
        
        return M;
}

 bool matrix_isEqual(struct matrix *A, struct matrix *B){
        
        if(A == NULL || B == NULL)
                alloc_error();
                // return false;
        int nrA = A->n_row,
            ncA = A->n_col,
            nrB = B->n_row,
            ncB = B->n_col,
            size= nrA * ncA;
        __data_type *pA = A->elements;
        __data_type *pB = B->elements;
        if(nrA != nrB || ncA != ncB){
                printf("Size differ: (%i, %i) != (%i, %i)\n", nrA, ncA, nrB, ncB);
        
        
        }

        bool ret = true;
        for(int i = 0; i < nrA; i++)
                for(int j=0; j<ncA; j++)
                        if(!IsEqual(pA[i*ncA+j], pB[i*ncA+j])){
                                printf("A[%d,%d] = %.8f != %.8f\n", i,j, pA[i*ncA+j], pB[i*ncA+j] );
                                //ret = false;
                                return false;
                                
                        }  
        return ret;
 }