#include <stdio.h>
#include <stdlib.h>

#include "Initialization.h"
#include "libmat.h"
#include "file_reader.h"

/*Here, we are going to check each output from each matrices 
  Obtained in C and Matlab*/

#define check_res(MT1, MT2) (matrix_isEqual(MT1,MT2)) ? printf("MATRIX %s OK !\n\n", #MT1): printf("MATRIX %s NOT OK !\n\n", #MT1);

int main()
{
    
    // FROM C
    struct matrix* output_bias_gyr = read_matrix("output/output_bias_gyr.csv");
    struct matrix* output_vel = read_matrix("output/output_vel.csv");
    struct matrix* output_bias_acc = read_matrix("output/output_bias_acc.csv");
    struct matrix* output_pos = read_matrix("output/output_pos.csv");
    struct matrix* output_qtn = read_matrix("output/output_qtn.csv");
    struct matrix* euler = read_matrix("output/euler.csv");
    
    // FROM Matblab (expected)
    struct matrix* mlab_bias_gyr = read_matrix("output/mlab_bias_gyr.csv");
    struct matrix* mlab_vel = read_matrix("output/mlab_vel.csv");
    struct matrix* mlab_bias_acc = read_matrix("output/mlab_bias_acc.csv");
    struct matrix* mlab_pos = read_matrix("output/mlab_pos.csv");
    struct matrix* mlab_qtn = read_matrix("output/mlab_qtn.csv");
    struct matrix* mlab_euler = read_matrix("DadosZurich/euler.csv");


    check_res(euler, mlab_euler);
    check_res(output_bias_gyr, mlab_bias_gyr);
    check_res(output_vel, mlab_vel);
    check_res(output_bias_acc, mlab_bias_acc);
    check_res(output_pos, mlab_pos);
    check_res(output_qtn, mlab_qtn);
    

    return 0;
}