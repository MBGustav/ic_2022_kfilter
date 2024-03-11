#ifndef __COMMON_MATRIX_H__
#define __COMMON_MATRIX_H__

#include <cmath>

#define  DATA_PRECISION == 1

#if defined DATA_PRECISION == 1
#define LD_ALIGN 256
#define LD_BIAS 8
#define data_t float
#else
#define LD_ALIGN 512
#define LD_BIAS 16
#define data_t double
#endif


// MAX WIDTH TO PRINT MATRIX
#define MAX_DISPLAY_MATRIX (10)


// #define JUMP_ASSERTION_COMANDS


// MEMORY ACCESS PATTERN - ROW MAJOR
#define idx_matrix(ld, i, j) (((i) * (ld)) + (j))
#define idx_matrix_t(ld, i, j) (((j) * (ld)) + (i))

#define alloc_error() { \
    fprintf(stderr, "Memory Error:Invalid Access.\n"); \
    exit(EXIT_FAILURE); \
}

#define LessThanOrEq(x, y) ({ \
    constexpr double epsilon = 0.000000000000001; \
    ((x) < (y) + epsilon) || (std::fabs((x) - (y)) < epsilon); \
})



// #TODO: set mkl_malloc or common malloc from libc
// #define MALLOC(x) (data_t *)mkl_malloc((x), 64)
#define MALLOC(x, size) x = (data_t*)malloc((size) * sizeof(size))
#define FREE(ptr) free(ptr)

// PADDING CONFIGURATION
#define HPL_PTR(ptr_, al_) ((((size_t)(ptr_) + (al_)-1) / (al_)) * (al_))
#if defined(PAD_LD)
static inline int ld_padding(int x) {
  int ld;
  ld = HPL_PTR(x, LD_ALIGN); // Rule 1
  if (ld - LD_BIAS >= x)
    ld -= LD_BIAS;
  else
    ld += LD_BIAS; // Rule 2
  return ld;
}
#else
static inline int ld_padding(int x) { return x; }
#endif



// MATH OPERATIONS
#define PI 3.14159265358979
#define TINY (1.0e-200)
#define GRAVITY -9.81
#define Time  1/freq   // s
#define Tau  100  // s
#define T2G (1) /*conversao de unidade: fluxo magnético de Tesla para Gauss */
#define constant_beta (-4.0 * (angle2rad))

#define sin_d(x) (sin((x)*angle2rad))
#define cos_d(x) (cos((x)*angle2rad))
#define atan_2d(x,y) (atan2(angle2r(x),angle2r(y)))

// ############################### Ellipsoid Definition ##########################

// % Informações do elipsóide de referência para o Sistema Geodésico Mundial de 
// % 1984 com unidade em metros
#define wg84_SemimajorAxis (6378137) //(meters)
#define wg84_Eccentricity  (0.0818191908426215) //(meters)

#define angle2rad (PI/180)
#define rad2angle(x) (180/PI)
#define angle2r(x)   ((x) * (PI/180))
#define r2angle(x)   ((x) * (180/PI))
#define EPSILON (0.000000005f)
#define IsEqual(x,y) (( ((x)-(y))*((x)-(y))  <= EPSILON * 10000)? 1 : 0)


#define LEFT_INTERVAL  0.00000000001
#define RIGHT_INTERVAL 0.000000001

#define FLOAT_COMPARE_WITH_INTERVAL(left, right) \
    ((fabs((right) - (left)) <= LEFT_INTERVAL) || (fabs((left) - (right)) <= RIGHT_INTERVAL))
    
#define FLT_EQ(x,y) (( ((x)-(y))*((x)-(y))  <= EPSILON)? 1 : 0)




#endif /*__COMMON_MATRIX_H__*/
