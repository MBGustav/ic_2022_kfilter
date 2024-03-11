#ifndef __MATMUL_HPP__
#define __MATMUL_HPP__

#include <iostream>
#include <iomanip>
#include <sycl/ext/intel/ac_types/ac_int.hpp>
#include <sycl/ext/intel/fpga_extensions.hpp>
#include <sycl/sycl.hpp>

#include "memory_transfers.hpp"
// Included from DirectProgramming/C++SYCL_FPGA/include/
#include "streaming_matmul.hpp"


#define MAX_DISPLAY_MATRIX (8)

#if not defined(IS_BSP)
using sycl::ext::intel::experimental::property::usm::buffer_location;
#endif

// Forward declare the kernel and pipe names
// (This prevents unwanted name mangling in the optimization report.)
class FeederA;
class FeederB;
class Matmul;
class Drain;
class APipe;
class BPipe;
class CPipe;
class DonePipe;

// Forward Declaration of functions
using PipeDataA = fpga_tools::NTuple<double, TILE_A>;
using PipeDataB = fpga_tools::NTuple<double, TILE_B>;
using PipeDataC = fpga_tools::NTuple<double, TILE_A>;

// Pipes to communicate the matrices between kernels
using PipeA = sycl::ext::intel::pipe<APipe, PipeDataA, 64>;
using PipeB = sycl::ext::intel::pipe<BPipe, PipeDataB, 64>;
using PipeC = sycl::ext::intel::pipe<CPipe, PipeDataC, 64>;
using PipeDone = sycl::ext::intel::pipe<DonePipe, bool, 64>;


void print_status(sycl::event feeder_a_event){
  sycl::info::event_command_status status = feeder_a_event.get_info<sycl::info::event::exe::command_execution_status>();
  std::cout << "Status de execução do comando: " << static_cast<int>(status) << std::endl;
}


void TransposeMat(std::vector<double> &m_matrix,
                  std::vector<double> &m_transposed, int rows, int cols,
                  int num_matrices);

void FillRand(std::vector<double> &m_matrix, int l_bound, int u_bound,
              int elements);
void ListedMatmul(matrix matrices_a[],
                  matrix matrices_b[], 
                  matrix matrices_c[], int num_matrices);
void PrintMat(std::vector<double> &m_matrix, int rows, int cols);
void MatmulRef(std::vector<double> &a_matrix, std::vector<double> &b_matrix,
               std::vector<double> &c_matrix, int rows_a, int common, int cols_b,
               int num_matrices);

/**
 * Implementation of the matrix multiplication using multiple streaming kernels.
 * Parameterized by datatype, matrix size, and tile size. Exercises the kernels
 * by running multiple repetitions for a set of matrices.
 *
 * Function arguments:
 *  q: device queue
 *  a_matrix: input matrix pointer (given in column-major)
 *  b_matrix: input matrix pointer (given in row-major, i.e., transposed)
 *  c_matrix: output matrix pointer (will be stored in column-major)
 *  repetitions: number of repetitions of the computation to execute
 *
 */
template <typename TT,          // Datatype of the elements of the matrix
          int rows_a,           // Rows of matrix A
          int common,           // Columns of matrix A / rows of matrix B
          int cols_b,           // Columns of matrix B
          int tile_a,           // Tile size for matrix A
          int tile_b,           // Tile size for matrix B
          bool transpose_b,     // Comp. Time - set to transpose matrix B
          int num_matrices>     // Number of pairs of matrices to multiply
void MatmulImpl(sycl::queue &q,            // Device queue
                std::vector<TT> &a_matrix, // Input matrix A
                std::vector<TT> &b_matrix, // Input matrix B
                std::vector<TT> &c_matrix, // Output matrix C = A * B
                int repetitions,           // Number of repetitions
                bool last_operation
) {
  // Number of elements per DDR access
  // NOTE: optimized for single-precision doubleing-point matrices
  constexpr int kElemsPerDDRAccess = 8;

  // Matrix sizes
  constexpr int kMatsizeA = rows_a * common;
  constexpr int kMatsizeB = cols_b * common;
  constexpr int kMatsizeC = rows_a * cols_b;

  // Buffer locations for mmhost interfaces
  constexpr int kBL1 = 1;
  constexpr int kBL2 = 2;
  constexpr int kBL3 = 3;

  // Allocate FPGA DDR memory
#if defined(IS_BSP)
  TT *a = sycl::malloc_device<TT>(kMatsizeA * num_matrices, q);
  TT *b = sycl::malloc_device<TT>(kMatsizeB * num_matrices, q);
  TT *c = sycl::malloc_device<TT>(kMatsizeC * num_matrices, q);
#else
  // malloc_device are not supported when targetting an FPGA part/family
  TT *a = sycl::malloc_shared<TT>(kMatsizeA * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL1)});
  TT *b = sycl::malloc_shared<TT>(kMatsizeB * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL2)});
  TT *c = sycl::malloc_shared<TT>(kMatsizeC * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL3)});
#endif
  std::cout <<"result in\n";
  // TT b_transposed[kMatsizeB * num_matrices];
  if constexpr (transpose_b){
    std::vector<TT> b_transposed(kMatsizeA * num_matrices);
    TransposeMat(b_matrix, b_transposed, cols_b, common, num_matrices);
    q.memcpy(b, b_transposed.data(), kMatsizeB * num_matrices * sizeof(TT)).wait();

  }else{
    q.memcpy(b, b_matrix.data(), kMatsizeB * num_matrices * sizeof(TT)).wait();
  }

  // Copy matrices over
    q.memcpy(a, a_matrix.data(), kMatsizeA * num_matrices * sizeof(TT)).wait();
  
    sycl::event feeder_a_event, feeder_b_event, drain_event;

  
  std::cout << "feedA\n";
  // Producer kernel for matrix A
  feeder_a_event = q.single_task<FeederA>(
      MatrixReadFromDDRToPipeA<TT, kBL1, rows_a, common, cols_b, tile_a, tile_b,
                               kElemsPerDDRAccess, num_matrices, PipeA,
                               PipeDone>{a, repetitions, last_operation});
  
  print_status(feeder_a_event);
      sycl::event::             
  std::cout << "feed\n";
  // Producer kernel for matrix B
  feeder_b_event = q.single_task<FeederB>(
      MatrixReadFromDDRToPipeB<TT, kBL2, rows_a, common, cols_b, tile_a, tile_b,
                               kElemsPerDDRAccess, num_matrices, PipeB>{
          b, repetitions});
  print_status(feeder_b_event);

  std::cout << "kernel\n";
  // Matrix multiply kernel
    q.single_task<Matmul>(
        fpga_linalg::StreamingMatmul<TT, common, tile_a, tile_b, PipeA, PipeB,
                                    PipeC, PipeDone>{});
  
  std::cout << "drain\n";
  // Consumer kernel for matrix C
  drain_event = q.single_task<Drain>(
      MatrixReadPipeToDDR<TT, kBL3, rows_a, cols_b, tile_a, tile_b,
                          kElemsPerDDRAccess, num_matrices, PipeC>{
          c, repetitions});
  
  print_status(drain_event);

    feeder_a_event.wait();
    print_status(feeder_a_event);
    std::cout << "1\n";
    feeder_b_event.wait();
    std::cout << "2\n";
    drain_event.wait();
    std::cout << "3\n";
    q.memcpy(c_matrix.data(), c, kMatsizeC * num_matrices * sizeof(TT)).wait();
  

  Compute the total time the execution lasted
  auto start_time = feeder_a_event.template get_profiling_info<
      sycl::info::event_profiling::command_start>();
  auto end_time = drain_event.template get_profiling_info<
      sycl::info::event_profiling::command_end>();
  double diff = (end_time - start_time) / 1.0e9;
  std::cout << "   Total duration:   " << diff << " s" << std::endl;
  std::cout << "Throughput: " << repetitions * num_matrices / diff * 1e-3
            << "k matrices/s" << std::endl;
  // Copy result matrix back

  // Free USM
  sycl::free(a, q);
  sycl::free(b, q);
  sycl::free(c, q);
}

template <typename TT,          // Datatype of the elements of the matrix
          int rows_a,           // Rows of matrix A
          int common,           // Columns of matrix A / rows of matrix B
          int cols_b,           // Columns of matrix B
          int tile_a,           // Tile size for matrix A
          int tile_b,           // Tile size for matrix B
          bool transpose_b,     // Comp. Time - set to transpose matrix B
          int num_matrices>     // Number of pairs of matrices to multiply
void Multiple(sycl::queue &q,            // Device queue
                std::vector<TT> &a_matrix, // Input matrix A
                std::vector<TT> &b_matrix, // Input matrix B
                std::vector<TT> &c_matrix, // Output matrix C = A * B
                int repetitions            // Number of repetitions
) {
  // Number of elements per DDR access
  // NOTE: optimized for single-precision doubleing-point matrices
  constexpr int kElemsPerDDRAccess = 8;

  // Matrix sizes
  constexpr int kMatsizeA = rows_a * common;
  constexpr int kMatsizeB = cols_b * common;
  constexpr int kMatsizeC = rows_a * cols_b;

  // Buffer locations for mmhost interfaces
  constexpr int kBL1 = 1;
  constexpr int kBL2 = 2;
  constexpr int kBL3 = 3;

  // Allocate FPGA DDR memory
#if defined(IS_BSP)
  TT *a = sycl::malloc_device<TT>(kMatsizeA * num_matrices, q);
  TT *b = sycl::malloc_device<TT>(kMatsizeB * num_matrices, q);
  TT *c = sycl::malloc_device<TT>(kMatsizeC * num_matrices, q);
#else
  // malloc_device are not supported when targetting an FPGA part/family
  TT *a = sycl::malloc_shared<TT>(kMatsizeA * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL1)});
  TT *b = sycl::malloc_shared<TT>(kMatsizeB * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL2)});
  TT *c = sycl::malloc_shared<TT>(kMatsizeC * num_matrices, q,
                                  sycl::property_list{buffer_location(kBL3)});
#endif
  std::cout <<"result in\n";
  // TT b_transposed[kMatsizeB * num_matrices];
  if constexpr (transpose_b){
    std::vector<TT> b_transposed(kMatsizeA * num_matrices);
    TransposeMat(b_matrix, b_transposed, cols_b, common, num_matrices);
    q.memcpy(b, b_transposed.data(), kMatsizeB * num_matrices * sizeof(TT)).wait();

  }else{
    q.memcpy(b, b_matrix.data(), kMatsizeB * num_matrices * sizeof(TT)).wait();
  }

  // Copy matrices over
  q.memcpy(a, a_matrix.data(), kMatsizeA * num_matrices * sizeof(TT)).wait();
  

  using PipeDataA = fpga_tools::NTuple<TT, tile_a>;
  using PipeDataB = fpga_tools::NTuple<TT, tile_b>;
  using PipeDataC = fpga_tools::NTuple<TT, tile_a>;

  // Pipes to communicate the matrices between kernels
  using PipeA = sycl::ext::intel::pipe<APipe, PipeDataA, 64>;
  using PipeB = sycl::ext::intel::pipe<BPipe, PipeDataB, 64>;
  using PipeC = sycl::ext::intel::pipe<CPipe, PipeDataC, 64>;
  using PipeDone = sycl::ext::intel::pipe<DonePipe, bool, 64>;

  // Producer kernel for matrix A
  auto feeder_a_event = q.single_task<FeederA>(
      MatrixReadFromDDRToPipeA<TT, kBL1, rows_a, common, cols_b, tile_a, tile_b,
                               kElemsPerDDRAccess, num_matrices, PipeA,
                               PipeDone>{a, repetitions});
  std::cout << "feed\n";
  // Producer kernel for matrix B
  auto feeder_b_event = q.single_task<FeederB>(
      MatrixReadFromDDRToPipeB<TT, kBL2, rows_a, common, cols_b, tile_a, tile_b,
                               kElemsPerDDRAccess, num_matrices, PipeB>{
          b, repetitions});

  std::cout << "kernel\n";
  // Matrix multiply kernel
  q.single_task<Matmul>(
      fpga_linalg::StreamingMatmul<TT, common, tile_a, tile_b, PipeA, PipeB,
                                   PipeC, PipeDone>{});
std::cout << "drain\n";
  // Consumer kernel for matrix C
  auto drain_event = q.single_task<Drain>(
      MatrixReadPipeToDDR<TT, kBL3, rows_a, cols_b, tile_a, tile_b,
                          kElemsPerDDRAccess, num_matrices, PipeC>{
          c, repetitions});
  // std::cout << "1\n";
  // feeder_a_event.wait();
  // std::cout << "1\n";
  // feeder_b_event.wait();
  // std::cout << "2\n";
  drain_event.wait();
  std::cout << "3\n";

  // Compute the total time the execution lasted
  // auto start_time = feeder_a_event.template get_profiling_info<
  //     sycl::info::event_profiling::command_start>();
  // auto end_time = drain_event.template get_profiling_info<
  //     sycl::info::event_profiling::command_end>();
  // double diff = (end_time - start_time) / 1.0e9;
  // std::cout << "   Total duration:   " << diff << " s" << std::endl;
  // std::cout << "Throughput: " << repetitions * num_matrices / diff * 1e-3
  //           << "k matrices/s" << std::endl;
  std::cout <<"result out\n";
  // Copy result matrix back
  q.memcpy(c_matrix.data(), c, kMatsizeC * num_matrices * sizeof(TT)).wait();

  // Free USM
  sycl::free(a, q);
  sycl::free(b, q);
  sycl::free(c, q);
}



void TransposeMat(std::vector<double> &m_matrix,
                  std::vector<double> &m_transposed, int rows, int cols,
                  int num_matrices) {
  int matsize = rows * cols;

  for (int matrix_idx = 0; matrix_idx < num_matrices; matrix_idx++) {
    for (int row = 0; row < rows; row++) {
      for (int col = 0; col < cols; col++) {
        m_transposed[matrix_idx * matsize + row * cols + col] =
            m_matrix[matrix_idx * matsize + col * rows + row];
      }
    }
  }
}

void FillRand(std::vector<double> &m_matrix, int l_bound, int u_bound,
              int elements) {
  for (int element = 0; element < elements; element++) {
    m_matrix[element] =
        static_cast<double>(rand()) /(static_cast<double>((RAND_MAX) / (u_bound - l_bound))) +
        l_bound;
  }
}



// Output a matrix to the screen (assumes column-major format).
void PrintMat(std::vector<double> &m_matrix, int rows, int cols) {
  for (int row = 0; row < std::min(MAX_DISPLAY_MATRIX, rows); row++) {
    for (int col = 0; col < std::min(MAX_DISPLAY_MATRIX, cols); col++) {
      // Copy old state of cout
      std::ios oldState(nullptr);
      oldState.copyfmt(std::cout);

      // Edit the output format of cout
      std::cout << std::fixed << std::setprecision(2);

      // Print the results
      std::cout << std::setw(8) << m_matrix[col * rows + row] << " ";

      // Restore the output format of cout
      std::cout.copyfmt(oldState);
    }
    std::cout << std::endl;
  }
}

// Multiply num_matrices pairs of matrices from a_matrix and b_matrix and store
// all the results in c_matrix.
void MatmulRef(std::vector<double> &a_matrix, std::vector<double> &b_matrix,
               std::vector<double> &c_matrix, int rows_a, int common, int cols_b,
               int num_matrices) {
  int matsize_a = rows_a * common;
  int matsize_b = cols_b * common;
  int matsize_c = rows_a * cols_b;

  for (int matrix_idx = 0; matrix_idx < num_matrices; matrix_idx++) {
    for (int col = 0; col < cols_b; col++) {
      for (int row = 0; row < rows_a; row++) {
        double sum = 0;
        for (int k = 0; k < common; k++) {
          sum += a_matrix[matrix_idx * matsize_a + k * rows_a + row] *
                 b_matrix[matrix_idx * matsize_b + col * common + k];
        }
        c_matrix[matrix_idx * matsize_c + col * rows_a + row] = sum;
      }
    }
  }
}


#endif /* __MATMUL_HPP__ */