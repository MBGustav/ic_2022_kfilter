
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <cstring>
#include "common_matrix.hpp"
#include "oneAPIMatrix.hpp"


// template<sycl::queue& queue>
oneAPIMatrix::oneAPIMatrix(sycl::queue &device_queue) :device_queue(device_queue), _data{NULL} {}

// template<sycl::queue& queue>
oneAPIMatrix::oneAPIMatrix(sycl::queue &device_queue, int row, int col) : device_queue(device_queue) ,_row(row), _col(col){allocate(row, col);}

// template<sycl::queue& queue>
oneAPIMatrix::~oneAPIMatrix(){
    deallocate();}

// template<sycl::queue& queue>
void oneAPIMatrix::allocate(int row,int col)
{
    this->_row = row;
    this->_col = col;
    // to make padding efficient
    this->_ld  = ld_padding(col);
    this->is_transposed=false;
    //allocating for shared memory -> move value through page fault
    this->_data = sycl::malloc_shared<data_t>(this->MemSize(false), this->dev());
}

// template<sycl::queue& queue>
void oneAPIMatrix::deallocate(){sycl::free(_data, this->dev());}

// template<sycl::queue& queue>
void oneAPIMatrix::resize(int row, int col)
{
    if (_row == row && _col == col) return;

    if(_data){
        this->deallocate(); //free space
    }

    // And if i just reduce the row and col if smaller? 
    this->allocate(row, col); 
}

// template<sycl::queue& queue>
int oneAPIMatrix::getRow() const {return this->_row;}


// Returns Col value 
// template<sycl::queue& queue>
int oneAPIMatrix::getCol() const { return this->_col;}


// Returns Leading Dimension width
// template<sycl::queue& queue>
int oneAPIMatrix::get_ld() const { return this->_ld;}

// template<sycl::queue& queue>
int oneAPIMatrix::getSize() const { return this->getCol() * this->getRow();}

// template<sycl::queue& queue>
size_t oneAPIMatrix::MemSize(bool bin) const{return get_ld() * getRow() * (bin ? sizeof(data_t) : 1 );}

// template<sycl::queue& queue>
void oneAPIMatrix::display()
{
    std::cout << " showing top left from oneAPIMatrix:\n";
    for (int row = 0; row < std::min(MAX_DISPLAY_MATRIX, this->getRow()); row++) {
    for (int col = 0; col < std::min(MAX_DISPLAY_MATRIX, this->getCol()); col++) {
      // Copy old state of cout
      std::ios oldState(nullptr);
      oldState.copyfmt(std::cout);

      // Edit the output format of cout
      std::cout << std::fixed << std::setprecision(2);

      // Print the results
      std::cout << std::setw(8) << this->at(col, row) << " ";

      // Restore the output format of cout
      std::cout.copyfmt(oldState);
    }
    std::cout << std::endl;
  }
}

// template<sycl::queue& queue>
void oneAPIMatrix::randMatrix(data_t lower_bound, data_t upper_bound)
{
    for(int ii=0; ii< getRow(); ii++)
    for(int jj=0; jj< getCol(); jj++)
        this->at(ii,jj) = 
            static_cast<data_t>(rand()) / 
            static_cast<data_t>(RAND_MAX / (upper_bound - lower_bound)) + lower_bound;
}

// template<sycl::queue& queue>
void oneAPIMatrix::Identity(data_t val)
{
    for(int ii=0; ii< getRow(); ii++)
    for(int jj=0; jj< getCol(); jj++)
        this->at(ii,jj) = ii == jj ? val : 0.0f;
}


// template<sycl::queue& queue>
data_t& oneAPIMatrix::data(){return *_data;}

data_t& oneAPIMatrix::idx(int x){
   #ifndef JUMP_ASSERTION_COMANDS
    assert(x >= 0);
    assert(x <  this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[x];
}
// template<sycl::queue& queue>
data_t& oneAPIMatrix::at(int x, int y){
    #ifndef JUMP_ASSERTION_COMANDS
    assert(idx_matrix(_ld, x, y)>=0);
    assert(idx_matrix(_ld, x, y) < this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[idx_matrix(_ld, x, y)];
}

// template<sycl::queue& queue>
data_t& oneAPIMatrix::at(int x, int y)const{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(idx_matrix(_ld, x, y)>=0);
    assert(idx_matrix(_ld, x, y) < this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[idx_matrix(_ld, x, y)];
}

// template<sycl::queue& queue>
void oneAPIMatrix::transposeInPlace()
{
}

// template<sycl::queue& queue>
oneAPIMatrix oneAPIMatrix::transpose()
{
    oneAPIMatrix result(this->dev(), this->getCol(), this->getRow());

    #pragma omp parallel for
    for(int i=0; i< this->getCol();i++)
        for(int j=0; j< this->getCol();j++)
            result.at(i,j) = this->at(j,i);
    return result;
}

// template<sycl::queue& queue>
void oneAPIMatrix::fill_diag(data_vec fill_values, bool reset_values)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(fill_values.size() == this->getCol());
    assert(fill_values.size() == this->getRow());
    #endif /*JUMP_ASSERTION_COMANDS*/

    if(reset_values) 
        this->zeros();

    for(int i=0; i < this->getCol(); i++)
        this->at(i, i) = fill_values[i];
}

// template<sycl::queue& queue>
void oneAPIMatrix::zeros()
{
    if(!_data)
        alloc_error();

    std::memset(this->_data, 0, this->MemSize(true));

}

