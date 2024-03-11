#include "Matrix.hpp"
#include "common_matrix.hpp"




Matrix::Matrix() : _data{NULL}{}
Matrix::Matrix(int row, int col) : _row(row), _col(col){allocate(row, col);}

Matrix::~Matrix(){deallocate();}

void Matrix::allocate(int row,int col)
{
    this->_row = row;
    this->_col = col;
    // to make padding efficient
    this->_ld  = ld_padding(col);
    
    MALLOC(this->_data, this->MemSize(true));

    this->allocated_space = this->MemSize(false);

    this->is_transposed=false;
}

void Matrix::deallocate()
{
    FREE(_data);
}

void Matrix::resize(int row, int col)
{
    //if same size, nothing to change
    if (_row == row && _col == col) return;

    //if smaller, just change the col and row values#TODO - check ld val

    // if greater, free space then set the new parameters

    if(_data)
        this->deallocate(); //free space

    this->allocate(row, col); 
}

void Matrix::mem_copy_to(Matrix &dest, int initial_pos, int final_pos)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(final_pos >= 0 && initial_pos >= 0);
    assert(final_pos - initial_pos <= dest.MemSize());

    #endif /*JUMP_ASSERTION_COMANDS*/

    memcpy(dest._data, this->_data , this->MemSize(true));
}

int Matrix::getRow() const { return this->_row;}

// Returns Col value 
int Matrix::getCol() const { return this->_col;}

// Returns Leading Dimension width
int Matrix::get_ld() const { return this->_ld;}

int Matrix::getSize() const { return this->getCol() * this->getRow();}

size_t Matrix::MemSize(bool bin) const{return get_ld() * getRow() * (bin ? sizeof(data_t) : 1 );}

void Matrix::display()
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(this->_data != NULL);
    #endif /*JUMP_ASSERTION_COMANDS*/
    
    std::cout << "\n showing top left from Matrix:\n";
    for (int row = 0; row < std::min(MAX_DISPLAY_MATRIX, this->getRow()); row++) {
    for (int col = 0; col < std::min(MAX_DISPLAY_MATRIX, this->getCol()); col++) {
      // Copy old state of cout
      std::ios oldState(nullptr);
      oldState.copyfmt(std::cout);

      // Edit the output format of cout
      std::cout << std::fixed << std::setprecision(2);

      // Print the results
      std::cout << std::setw(8) << this->at(row, col) << " ";

      // Restore the output format of cout
      std::cout.copyfmt(oldState);
    }
    std::cout << std::endl;
  }
}
void Matrix::fill_matrix(Matrix &in, int offset_x, int offset_y)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(in.getRow() + offset_x <= this->getRow());
    assert(in.getCol() + offset_y <= this->getCol());
    assert(offset_x >= 0);
    assert(offset_y >= 0);
    #endif /*JUMP_ASSERTION_COMANDS*/
    
    for(int x=0; x < in.getRow(); x++ )
        for(int y=0; y < in.getCol(); y++)
            this->at(offset_x + x,offset_y + y) = in.at(x,y);
}

void Matrix::randMatrix(data_t lower_bound, data_t upper_bound)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(this->_data != NULL);
    #endif /*JUMP_ASSERTION_COMANDS*/

    for(int ii=0; ii< getRow(); ii++)
    for(int jj=0; jj< getCol(); jj++)
        this->at(ii, jj) = 
            static_cast<data_t>(rand()) / 
            static_cast<data_t>(RAND_MAX / (upper_bound - lower_bound)) + lower_bound;
}

void Matrix::Identity(data_t val)
{
    for(int ii=0; ii< getRow(); ii++)
    for(int jj=0; jj< getCol(); jj++)
        this->at(ii,jj) = ii == jj ? val : 0.0f;
}



data_t& Matrix::data(){return *_data;}
data_t& Matrix::at(int x, int y){
    #ifndef JUMP_ASSERTION_COMANDS
    assert(idx_matrix(_ld, x, y)>=0);
    assert(idx_matrix(_ld, x, y) < this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/

    return _data[idx_matrix(_ld, x, y)];
}

Matrix Matrix::copy_of()
{
    Matrix copy_A(this->getRow(), this->getCol());
    for(int i = 0; i < this->getRow(); i++)
    for(int j = 0; j < this->getCol(); j++)
        copy_A.at(i,j) = this->at(i,j);

    return copy_A;
}

data_t& Matrix::at(int x, int y)const{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(idx_matrix(_ld, x, y)>=0);
    assert(idx_matrix(_ld, x, y) < this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[idx_matrix(_ld, x, y)];
}


data_t& Matrix::idx(int x){
   #ifndef JUMP_ASSERTION_COMANDS
    assert(x >= 0);
    assert(x <  this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[x];
}

data_t& Matrix::idx(int x) const{
   #ifndef JUMP_ASSERTION_COMANDS
    assert(x >= 0);
    assert(x <  this->MemSize(false));
    #endif /*JUMP_ASSERTION_COMANDS*/
    return _data[x];
}

void Matrix::transposeInPlace()
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(_data != NULL);
    #endif /*JUMP_ASSERTION_COMANDS*/

    this->swap_RowCol();

    // If a matrix is squared, we do less operations:
    if(this->getRow() == this->getCol()){
        for(int kRow = 0; kRow < this->getRow(); kRow++)
            for(int kCol = kRow+1; kCol < this->getCol(); kCol++) //we dont need to transp the diagonal (+1)
                std::swap(this->at(kRow, kCol), this->at(kCol, kRow));
        return;
    }

    for(int kRow = 0; kRow < this->getRow(); kRow++)
        for(int kCol = 0; kCol < this->getCol(); kCol++)
            std::swap(this->at(kRow, kCol), this->at(kCol, kRow));


}

void Matrix::scal_multply(data_t val)
{
    for(int idx = 0; idx <this->MemSize(false); idx++)
        this->idx(idx) *= val;
}

Matrix Matrix::transpose()
{
    Matrix result(this->getCol(), this->getRow());

    #pragma omp parallel for
    for(int i=0; i< this->getCol();i++)
        for(int j=0; j< this->getCol();j++)
            result.at(i,j) = this->at(j,i);
    return result;
}

void Matrix::fill_diag(std::vector<data_t> fill_values, bool reset_values)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(fill_values.size() == this->getCol());
    assert(fill_values.size() == this->getRow());
    #endif /*JUMP_ASSERTION_COMANDS*/

    if(reset_values) 
        this->zeros();

    for(int i=0; i < this->getCol(); i++)
        this->at(i,i) = fill_values[i];
}

Matrix Matrix::get_collumns(int kCol)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(kCol >= 0);
    assert(kCol < this->getRow());
    #endif /*JUMP_ASSERTION_COMANDS*/

    Matrix A(1, this->getRow());
    for(int idx=0; idx < this->getRow(); idx++){
        A.idx(idx) = this->at(idx, kCol);
    }
    return A;
}

Matrix Matrix::get_rows(int kRow)
{
    #ifndef JUMP_ASSERTION_COMANDS
    assert(kRow >= 0);
    assert(kRow < this->getCol());
    #endif /*JUMP_ASSERTION_COMANDS*/
    Matrix A(1, this->getCol());
    for(int idx=0; idx < this->getCol(); idx++){
        A.idx(idx) = this->at(kRow, idx);
    }
    return A;

}

data_t Matrix::norm() const
{

    data_t acc  = 0.0f; 
    #pragma omp parallel for reduction(+:acc)
    for(int i = 0; i< this->getSize();i++)
    {
        acc += this->idx(i)*this->idx(i);
    }
    return sqrt(acc);
}

void Matrix::zeros()
{
    if(!_data)
        alloc_error();

    std::memset(this->_data, 0, this->MemSize(true));

}

