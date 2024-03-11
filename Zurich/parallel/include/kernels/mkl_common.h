#ifndef _MKL_COMMON_HPP_
#define _MKL_COMMON_HPP_

#include <stdlib.h>
#include <iostream>
#include <algorithm> 
#include <list>
#include <map>
#include <type_traits>



#ifdef MKL_ONEAPI
    #include <sycl/sycl.hpp>
    #include "oneapi/mkl.hpp"
#else
    #include "mkl.h"
#endif

//usr libs
#include "structs.h"


#define DATA_PRECISION (2)

//data define
#if DATA_PRECISION == 1
#define __data_type double
#else
#define __data_type float
#endif


//performance padding
#if DATA_PRECISION == 1
#define LD_ALIGN 256
#define LD_BIAS 8
#else
#define LD_ALIGN 512
#define LD_BIAS 16
#endif


//Macros for matrices
#define HPL_PTR(ptr_, al_) ((((size_t)(ptr_) + (al_)-1) / (al_)) * (al_))
#define index(i, j, ld) (((j) * (ld)) + (i))
#define MALLOC(x) (__data_type *)mkl_malloc((x), 64);
#define MALLOC_CHECK(p)                                                        \
  if (p == NULL) {                                                             \
    fprintf(stderr, "%s:%d: memory allocation error\n", __FILE__, __LINE__);   \
    exit(EXIT_FAILURE);                                                       \
  }


// padding for optimal performance ??
inline int nice_ld(int x)
{
    x = std::max(x, 1);
    x *= sizeof(__data_type);
    x = (x + 511) & ~511;
    x += 256;
    x /= sizeof(__data_type);
    return x;
}


#if defined(PAD_LD)
static inline int getld(int x) {
  int ld;
  ld = HPL_PTR(x, LD_ALIGN); // Rule 1
  if (ld - LD_BIAS >= x)
    ld -= LD_BIAS;
  else
    ld += LD_BIAS; // Rule 2
  return ld;
}
#else
static inline int getld(int x) { return x; }
#endif

// ###############################################
// #Matrix Declarations for MKL
// (We consider the archs as collumn -major)
// ###############################################

struct mkl_matrix
{
    #ifdef MKL_ONEAPI
    sycl::queue q;

    mkl_matrix(sycl::queue &main_queue) : q{main_queue}, data{NULL} {};
    void deallocate(){if (data) free(data, q.get_context());}
    #else
    void deallocate(){if (data) free(data);}
    #endif

    __data_type *data;
    int row, col, ldw;
    int size;
    mkl_matrix() : data{NULL} {};
    
    void allocate(int _row, int _col)
    {
        if(_row <1 || _col<1 ) 
        {
            fprintf(stderr, "%s:%d: memory allocation error\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }

        row = _row;
        col = _col;
        ldw = getld(col); // leading dimension
        size = ldw*row;
        // deallocate();
        #ifdef MKL_ONEAPI
        data = (__data_type *)sycl::malloc_shared(sizeof(__data_type) * size, q.get_device(), q.get_context());
        #else
        data = MALLOC(sizeof(__data_type) * size);
        #endif
        MALLOC_CHECK(data);
    }

    void rnd_matrix(bool rand = true)
    {
        for(int i=0; i < row; i++)
            #pragma simd vectorlength(8)
            for(int j=0; j < col; j++)
                data[index(i, j, ldw)] = (rand ? (__data_type) (std::rand() %100) : i*col + j);
    }

    void printMat(){
        int i, j;
        printf("Printing top left:\n");
        for (i = 0; i < std::min(6, row); i++) {
            printf("\n");
            for (j = 0; j < std::min(6, col); j++)
                printf("  %010.4lf", data[index(i, j, ldw)]);
                
        }
        printf("\n\n");
        
    }

    __data_type* data_NoPadding()
    {
        __data_type *data_nopad = (__data_type*) malloc(row * col * sizeof(__data_type));
        //copy data changing leading dimension
        for(int i = 0; i < row; i++)
        for(int j = 0; j < col; j++)
            data_nopad[index(i,j,col)] = data[index(i,j,ldw)];
            
        return data_nopad;
    }
    
    __data_type& begin(){return *data;}
    
    // ~mkl_matrix() { deallocate(); }
};

// ###############################################


// ###############################################
// Data Management to gemm, potri and other 
// operations
// ###############################################

typedef struct gemm_conf{
    int m, k, n; //matrix parameters
    int lda, ldb,ldc; 
    // int sizeA, sizeB, sizeC; // matriz size
#ifdef MKL_ONEAPI
    oneapi::mkl::transpose trA, trB;
#else
    CBLAS_TRANSPOSE trA, trB;
#endif
}gemm_conf;

static inline gemm_conf set_gemm_conf(struct mkl_matrix &A, struct mkl_matrix &B, bool transA, bool transB)
{   

    int m  = (transA) ? A.col : A.row;
    int n  = (transB) ? B.row : B.col;
    int k  = (transA) ? A.row : A.col;
    
    int lda = A.ldw;
    int ldb = B.ldw;
    int ldc = getld(n); //if PAD, switch ldC output

    // std::cout << "ld a,b,c: " << lda <<", " << ldb << ", " << ldc <<std::endl;
    // std::cout << "m,n,k   : " << m <<", " << n << ", " << k <<std::endl;
    // std::cout << "sizes   : " << sizeA << ", " << sizeB << ", " << sizeC << std::endl;
#ifdef MKL_ONEAPI
    oneapi::mkl::transpose trA = (transA) ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::nontrans;
    oneapi::mkl::transpose trB = (transB) ? oneapi::mkl::transpose::trans : oneapi::mkl::transpose::nontrans;
#else
    CBLAS_TRANSPOSE trA = (transA) ? CblasTrans : CblasNoTrans ;
    CBLAS_TRANSPOSE trB = (transB) ? CblasTrans : CblasNoTrans ;
#endif
    gemm_conf conf = {.m=m , .k=k, .n=n, 
                    .lda=lda, .ldb=ldb, .ldc=ldc,
                    .trA=trA, .trB=trB};
    return conf;
}


#ifdef MKL_ONEAPI
enum my_sycl_device_types {cpu_device, gpu_device};

std::map<my_sycl_device_types, std::string>
sycl_device_names = { {cpu_device, "CPU"},
                      {gpu_device, "GPU"}  };

//
// Users can set environment flags like SYCL_DEVICES_{all,cpu,gpu} to specify
// which devices to run on
//
void set_list_of_devices(std::list<my_sycl_device_types> & list_of_devices)
{
#if defined(SYCL_DEVICES_all)
    list_of_devices.push_back(cpu_device);
    list_of_devices.push_back(gpu_device);
#else
#if defined(SYCL_DEVICES_cpu)
    list_of_devices.push_back(cpu_device);
#endif
#if defined(SYCL_DEVICES_gpu)
    list_of_devices.push_back(gpu_device);
#endif
#endif
}

void get_sycl_device(sycl::device &my_dev, bool & my_dev_is_found, my_sycl_device_types & desired_sycl_device)
{

    my_dev_is_found = true;

    try {

        switch (desired_sycl_device) {
            case cpu_device:
                my_dev = sycl::device(sycl::cpu_selector_v);
                break;
            case gpu_device:
                my_dev = sycl::device(sycl::gpu_selector_v);
                break;
        }

    } catch (...) {
        my_dev_is_found = false;
    }

}

std::list<my_sycl_device_types> get_list_of_found_devices(bool failDueToMissingDevice)
{
    std::list<my_sycl_device_types> listOfDevices;
    std::list<my_sycl_device_types> listOfFoundDevices;
    set_list_of_devices(listOfDevices);

    failDueToMissingDevice = false;
    for (auto &deviceType : listOfDevices) {
        sycl::device myDev;
        bool myDevIsFound = false;
        get_sycl_device(myDev, myDevIsFound, deviceType);
        if (myDevIsFound) {
            listOfFoundDevices.push_back(deviceType);
        } else {
#ifdef FAIL_ON_MISSING_DEVICES
        std::cout << "No " << sycl_device_names[deviceType] << " devices found; Fail on missing devices is enabled.\n";
        failDueToMissingDevice = true;
#else
        std::cout << "No " << sycl_device_names[deviceType] << " devices found; skipping " << sycl_device_names[deviceType] << " tests.\n";
#endif
        }
    }
    return listOfFoundDevices;
}

#endif




#endif //_MKL_COMMON_HPP_