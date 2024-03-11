#include <iomanip>
#include <iostream>

#include <sycl/ext/intel/fpga_extensions.hpp>
#include <sycl/sycl.hpp>

#include "matrix.hpp"
#include "exception_handler.hpp"
#include "matmul.hpp"



// Compares num_matrices pairs of matrices; returns true iff they are equal
// given a tolerated error bound.
bool EqualMat(std::vector<double> &c_matrix, std::vector<double> &c_reference,
              int rows, int cols, int num_matrices) {
  int matsize = rows * cols;
  bool passed = true;

  // doubleing-point error threshold value
  constexpr double kEpsilon = 0.01f;

  for (int matrix_idx = 0; matrix_idx < num_matrices; matrix_idx++) {
    for (int col = 0; col < cols; col++) {
      for (int row = 0; row < rows; row++) {
        int idx = matrix_idx * matsize + col * rows + row;
        if (abs(c_matrix[idx] - c_reference[idx]) > kEpsilon) {
          passed = false;
#if DEBUG
          std::cout << "Error: C[" << col << "][" << row << "] = "
                    << c_matrix[idx]
                    << " but REF[" << col << "][" << row << "] = "
                    << c_reference[idx] << std::endl;
#endif
        }
        if (!std::isfinite(c_matrix[idx])) {
          passed = false;
#if DEBUG
          std::cout << "C[" << col << "][" << row << "] = " << c_matrix[idx]
                    << " is not finite" << std::endl;
#endif
        }
      }
    }
  }
  return passed;
}


// Transpose num_matrices matrices in m_matrix and store the results in
// m_transposed.


int main(int argc, char *argv[]) {
  // Matrix paramters specified by build system
  constexpr int kRows = ROWS_A;
  constexpr int kCols = COLS_B;
  
  /*TODO:will be taken off*/
  constexpr int kCommon = COMMON;
  
  constexpr int kTileA = TILE_A;
  constexpr int kTileB = TILE_B;

  // Matrix size - standard for FPGA
  constexpr int kMatsize = kRows * kCommon;
  constexpr bool TransB = true;

  // Repetitions and number of matrices to measure performance
#if FPGA_SIMULATOR
  int repetitions = argc > 1 ? atoi(argv[1]) : 1;
  constexpr int kNumMatrices = 4;
#elif FPGA_HARDWARE
  int repetitions = argc > 1 ? atoi(argv[1]) : 819200;
  constexpr int kNumMatrices = 4;
#else // #if FPGA_EMULATOR
  int repetitions = argc > 1 ? atoi(argv[1]) : 16;
  constexpr int kNumMatrices = 4;
#endif

  try {

#if FPGA_SIMULATOR
  auto selector = sycl::ext::intel::fpga_simulator_selector_v;
#elif FPGA_HARDWARE
  auto selector = sycl::ext::intel::fpga_selector_v;
#else // #if FPGA_EMULATOR
  auto selector = sycl::ext::intel::fpga_emulator_selector_v;
#endif

  // Enable the queue profiling to time the execution
  sycl::property_list queue_properties{
      sycl::property::queue::enable_profiling()};
  sycl::queue q =
      sycl::queue(selector, fpga_tools::exception_handler, queue_properties);

  std::cout << "Running on device: "
            << q.get_device().get_info<sycl::info::device::name>().c_str()
            << std::endl;

  // Create arrays to hold the input and output matrices
  std::vector<double> a_matrix(kMatsize * kNumMatrices);
  std::vector<double> b_matrix(kMatsize * kNumMatrices);
  std::vector<double> c_matrix(kMatsize * kNumMatrices);
  std::vector<double> d_matrix(kMatsize * kNumMatrices);

  // Generate random A and B matrices
  srand(1138);
  // FillRand(a_matrix, kRandMin, kRandMax, kMatsize * kNumMatrices);
  // FillRand(b_matrix, kRandMin, kRandMax, kMatsizeB * kNumMatrices);

  matrix MatrixArray_A[kNumMatrices];
  matrix MatrixArray_B[kNumMatrices];
  matrix MatrixArray_C[kNumMatrices];

  
  matrix C1_ref(4,4), C2_ref(6,6);
  matrix MatrixB1t(4, 4);
  matrix MatrixB2t(6, 6);

  for(int idx_matrix =0;idx_matrix < kNumMatrices; idx_matrix++ ){
  
    // set pointer to the vector offload in FPGA
    MatrixArray_A[idx_matrix].elem = a_matrix.data() + idx_matrix * kMatsize;
    MatrixArray_B[idx_matrix].elem = b_matrix.data() + idx_matrix * kMatsize;
    MatrixArray_C[idx_matrix].elem = c_matrix.data() + idx_matrix * kMatsize;

    //Generate random squared matrices sizes
    int width =  1 + rand()%(6 - 1);
    MatrixArray_A[idx_matrix].row = MatrixArray_A[idx_matrix].col = width;
    MatrixArray_B[idx_matrix].row = MatrixArray_B[idx_matrix].col = width;
    MatrixArray_C[idx_matrix].row = MatrixArray_C[idx_matrix].col = width;

    RandMatrix(MatrixArray_A[idx_matrix]);
    RandMatrix(MatrixArray_B[idx_matrix]);
    RandMatrix(MatrixArray_C[idx_matrix]);
  }



  // Calculate a reference to compare our answer to and store it in c_reference
  // NOTE: since the systolic matrix multiply interprets B as transposed, we
  // need to first transpose b_matrix to b_transposed to use it in the standard
  // MM algorithm
  std::vector<double> c_reference(kMatsize * kNumMatrices);
    std::vector<double> b_transposed;

  if constexpr (TransB){
    b_transposed.resize(kMatsize * kNumMatrices);
    TransposeMat(b_matrix, b_transposed, kCols, kCommon, kNumMatrices);
  }else{
    b_transposed = b_matrix;
  }


  MatmulRef(a_matrix, b_matrix, c_reference, kRows, kCommon, kCols,
            kNumMatrices);

  std::cout << " Matrix A size: " << kRows << " x " << kCommon
            << " (tile: " << kTileA << " x " << kCommon << ")" << std::endl
            << " Matrix B size: " << kCommon << " x " << kCols
            << " (tile: " << kCommon << " x " << kTileB << ")" << std::endl
            << " Systolic array size: " << kTileA << " x " << kTileB << " PEs"
            << std::endl;
  std::cout << "Running matrix multiplication of " << kNumMatrices
            << ((kNumMatrices > 1) ? " matrices " : " matrix ") << repetitions
            << " times" << std::endl;

    
  for(int i =0; i < 2; i++ )

  // Run the matrix multiplication
  MatmulImpl<double, kRows, kCommon, kCols, kTileA, kTileB, TransB, kNumMatrices>
  (q, a_matrix, b_matrix, c_matrix, repetitions, false);

    // MatmulImpl<double, kRows, kCommon, kCols, kTileA, kTileB, TransB, kNumMatrices>(
    //   q, a_matrix, c_matrix, d_matrix, repetitions, true);

#define DEBUG 1
#ifdef DEBUG
  // Print A, B, C and reference matrices
  for (int matrix_idx = 0; matrix_idx < kNumMatrices; matrix_idx++) {
    std::cout << std::endl << matrix_idx << std::endl;

    std::cout << std::endl << "Matrix A("
              << MatrixArray_A[matrix_idx].row <<", "
              << MatrixArray_A[matrix_idx].col <<") " << std::endl;

    std::vector<double> a_vector = {
        a_matrix.begin() + matrix_idx * kMatsize,
        a_matrix.begin() + (matrix_idx + 1) * kMatsize};
    PrintMat(a_vector, kRows, kCommon);

    std::cout << std::endl << "Matrix B("
              << MatrixArray_B[matrix_idx].row <<", "
              << MatrixArray_B[matrix_idx].col <<") " << std::endl;
    std::vector<double> b_vector = {
        b_transposed.begin() + matrix_idx * kMatsize,
        b_transposed.begin() + (matrix_idx + 1) * kMatsize};
    PrintMat(b_vector, kCommon, kCols);

    std::cout << std::endl << "Matrix C reference" << std::endl;
    std::vector<double> c_ref_vector = {
        c_reference.begin() + matrix_idx * kMatsize,
        c_reference.begin() + (matrix_idx + 1) * kMatsize};
    PrintMat(c_ref_vector, kRows, kCols);

    std::cout << std::endl << "Matrix C calculated" << std::endl;
    std::vector<double> c_vector = {
        c_matrix.begin() + matrix_idx * kMatsize,
        c_matrix.begin() + (matrix_idx + 1) * kMatsize};
    PrintMat(c_vector, kRows, kCols);
  }
#endif

  // Verify results
  bool passed = EqualMat(c_matrix, c_reference, kRows, kCols, kNumMatrices);
  std::cout << std::endl << (passed ? "PASSED" : "FAILED") << std::endl;

  return !passed;

  } catch (sycl::exception const &e) {
    std::cerr << "Caught a synchronous SYCL exception: " << e.what()
              << std::endl;
    std::cerr << "   If you are targeting an FPGA hardware, "
                 "ensure that your system is plugged to an FPGA board that is "
                 "set up correctly"
              << std::endl;
    std::cerr << "   If you are targeting the FPGA emulator, compile with "
                 "-DFPGA_EMULATOR"
              << std::endl;

    std::terminate();
  }
} // end of main