#ifndef __MKLMATRIX_HPP__
#define __MKLMATRIX_HPP__

// #include "Matrix.hpp"
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <cstring>
#include <CL/sycl.hpp>
#include "common_matrix.hpp"


typedef std::vector<data_t> data_vec;

// WHY create an another class ?
// We don't want to use virtual method mainly because the 
// cost of run time execution. Then we copy and add the necessary parameters


class oneAPIMatrix{
private:
    int _row, _col, _ld;
    bool is_transposed; 
    sycl::queue &device_queue;

    data_t *_data;

    void allocate(int row, int col);
    void deallocate();
public:
    sycl::queue &dev(){return this->device_queue;}
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
    data_t& idx(int pos);
    data_t& at(int x, int y) const;
    data_t&  data();

    void transposeInPlace(); 

    oneAPIMatrix transpose();

    void fill_diag(data_vec fill_values, bool reset_values = true);
    void zeros();

    void randMatrix(data_t lower_bound, data_t upper_bound);
    void Identity(data_t val = 1.0f);

    // Constructors and Destructors
    oneAPIMatrix(sycl::queue &device_queue, int row, int col);
    oneAPIMatrix(sycl::queue &device_queue);
    ~oneAPIMatrix();

};



#endif /*__MKLMATRIX_HPP__*/