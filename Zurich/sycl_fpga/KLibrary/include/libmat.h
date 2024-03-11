// Prototypes of libmat.c
#ifndef LIBMAT_H
#define LIBMAT_H

#include "structs.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

struct matrix *matrix_m(int n_row, int n_col);
void delete_m(struct matrix *M);
void print_m(struct matrix *M);
 struct matrix* rnd_matrix(int row, int col);
void eye_m(double value, struct matrix *M);
short int equal_m(struct matrix *A, struct matrix *B);
void equalpart_m(struct matrix *A, int i_0, int j_0, struct matrix *B);
void getpart_m(struct matrix *A, int i_0, int j_0, struct matrix *B);
short int sum_m(struct matrix *A, struct matrix *B, struct matrix *result);
short int less_m(struct matrix *A, struct matrix *B, struct matrix *result);
short int times_m(struct matrix *A, struct matrix *B, struct matrix *result);
short int timesc_m(double value, struct matrix *A, struct matrix *result);
short int transp_m(struct matrix *A, struct matrix *result);
short int inv_m(struct matrix *A, struct matrix *result);
void givens_m(struct matrix *A, struct matrix *result);
void euler2quat_m(struct euler angles, struct matrix *quat);


short int crossM_f(struct matrix *A, struct matrix *M);
void mat_cpy(struct matrix *src, struct matrix* dst);
void mat_cpyRow(struct matrix *src, int r_src ,struct matrix* dst, int r_dst);
double norm_m(struct matrix *v);
int maxIdx_v(double *A, int size_max);

void getpart_mat(struct matrix *src, int i_0, int j_0, struct matrix *dst);
void getpart_vec(struct matrix *src, int init, int end, struct matrix *dst);

void mat_size(struct matrix* A);

void unsafe_sum(struct matrix *A, struct matrix *B, struct matrix*C);
void unsafe_timesc(double v, struct matrix *A, struct matrix *B);
void unsafe_less(struct matrix *A, struct matrix *B, struct matrix*C);


void geodetic2ecef(struct matrix* llh, struct matrix* xyz);
void ecef2ned(struct matrix* pe, struct matrix *llh0, struct matrix *pt);

void geodetic2ned(struct matrix* llh, struct matrix* llh0, struct matrix *ned);

void save_matrix(struct matrix *M,const char *name);
struct matrix* read_matrix(char* name);
bool matrix_isEqual(struct matrix *A, struct matrix *B);
#endif /*LIBMAT_H*/