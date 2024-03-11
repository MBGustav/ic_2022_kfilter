#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "common_matrix.hpp"



typedef std::vector<data_t> data_vec;

class Matrix{
private:
    int _row, _col, _ld;
    bool is_transposed; 

    data_t *_data;
    int allocated_space;
    void allocate(int row, int col);
    void deallocate();

    void swap_RowCol(){std::swap(this->_row, this->_col);};

public:

    void display();

    //change dimension
    void resize(int row, int col);

    // Returns row value
    int getRow() const;

    // Returns Col value 
    int getCol() const;
    
    // Returns Leading Dimension width
    int get_ld() const;
    
    // Total Amount of (real) data
    int getSize() const;
    
    size_t MemSize( bool=false) const;
    
    // device acess position
    data_t& at(int x, int y);
    data_t& idx(int x);
    data_t& idx(int x) const;
    data_t& at(int x, int y) const;
    data_t&  data();

    void Matrix::scal_multply(data_t val);

    void transposeInPlace(); 

    Matrix Matrix::copy_of();
    
    Matrix transpose();

    data_t norm() const;

    void fill_diag(std::vector<data_t> fill_values, bool reset_values = true);
    
    void zeros();

    Matrix get_collumns(int kCol);
    Matrix get_rows(int kRow);
    
    void fill_matrix(Matrix &in, int offset_x =0, int offset_y = 0);

    void mem_copy_to(Matrix &dest, int initial_pos=0, int final_pos=0);

    void randMatrix(data_t lower_bound, data_t upper_bound);

    void Identity(data_t val = 1.0f);

    Matrix(int row, int col);
    Matrix();
    ~Matrix();

};



#endif /*_MATRIX_HPP_*/