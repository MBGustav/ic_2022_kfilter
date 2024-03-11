#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstdio>
#include <cassert> /*using assert*/
#include <cstdlib> /*using calloc*/
#include <iostream>
#include <iomanip>
#include <algorithm>

#ifndef MAX_DISPLAY_MATRIX
#define MAX_DISPLAY_MATRIX (8)
#endif
using data_type = double;

#define assertm(exp, msg) assert(((void)msg, exp))


constexpr int sizeCol = COLS_B;
constexpr int sizeRow = ROWS_A;
constexpr int sizeMatrix = sizeCol * sizeRow;


// Declaring a matrix for our use 


typedef struct matrix{
    int row, col; data_type *elem;

//  Constructors and Destructors
    matrix(int kSizeRow, int kSizeCol){allocate(kSizeRow, kSizeCol);}
    matrix(): elem{NULL}{}
    ~matrix(){deallocate();}

    //functions pointer return for data elem
    data_type* data(){return elem;}
    data_type* data() const {return elem;}
    
    // position calculator
    data_type& at(int i, int j){return elem[sizeCol * i + j];}
    data_type& at(int i, int j) const{return elem[sizeCol * i + j];}

    //alocate space to the matrix
    void allocate(int kSizeRow, int kSizeCol){
        
        // check size is capable to handle in FPGA
        assert((kSizeCol <= sizeCol) && (kSizeRow <= sizeRow));
        row = kSizeRow;
        col = kSizeCol;

        // Initializes with zeros
        elem = new data_type[sizeMatrix]();
    }

    // dealloc space in memory
    void deallocate(){delete[] elem;}
}matrix;



void transposed(matrix& matrix_a, matrix& matrix_a_transposed)
{
    assert(matrix_a.row == matrix_a_transposed.col &&
           matrix_a.col == matrix_a_transposed.row);

    for(int i = 0; i < matrix_a.row; i++)
        for(int j = 0; j < matrix_a.col; j++)
            matrix_a_transposed.at(i, j) = matrix_a.at(j, i);
}

void RandMatrix(matrix& matrix_a) {
    int row = matrix_a.row;
    int col = matrix_a.col;
    for(int i = 0; i < row; i++)
        for(int j = 0; j < col; j++)
        matrix_a.at(i, j) = std::round(static_cast<data_type>(rand()) / 
                            static_cast<data_type>(RAND_MAX) * 100.0f);
}

void IdentityMatrix(matrix& matrix_a, data_type scal = 1.0f)
{
    for(int i=0; i<matrix_a.row; i++)
        for (int j = 0; j < matrix_a.col; j++)
            matrix_a.at(i,j) = (i==j) ? scal : 0.0;
}

void displayMatrix(matrix& matrix_a)
{
    std::cout << "\nDisplaying matrix\n";
    for(int row = 0; row < std::min(MAX_DISPLAY_MATRIX, matrix_a.row); row++){
        for(int col = 0; col< std::min(MAX_DISPLAY_MATRIX, matrix_a.col); col++)
        {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << std::setw(8) << matrix_a.at(row, col) << " ";
        }
        std::cout << std::endl;
    }
}

void ListedMatmul(matrix matrices_a[],
                  matrix matrices_b[], 
                  matrix matrices_c[], int num_matrices)
{

  for (int matrix_idx = 0; matrix_idx < num_matrices; matrix_idx++) 
  {
    int cols_b = matrices_b[matrix_idx].col;
    int rows_a = matrices_a[matrix_idx].row;
    int common = matrices_a[matrix_idx].col;
    for (int col = 0; col < cols_b; col++){
      for (int row = 0; row < rows_a; row++){
        float sum = 0;
        for (int k = 0; k < common; k++) {
          sum += matrices_a[matrix_idx].at(k, row)*
                matrices_b[matrix_idx].at(col, k); 
        }
        matrices_c[matrix_idx].at(col,row) = sum;
      }
    }
  }
}
#endif /*MATRIX_HPP*/