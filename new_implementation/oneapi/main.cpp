#include "Matrix.hpp"
// #include "oneAPIMatrix.hpp"
#include <CL/sycl.hpp>
#include "KalmanFilter.hpp"
// #include "Quaternion.hpp"
#include "Linalg.hpp"
int main()
{
    
    try{
        sycl::queue Q(sycl::gpu_selector_v);
    } catch(...){
        std::cout << "GPU Not Found, trying to use default device...\n";
        sycl::queue Q(sycl::default_selector_v);    
    }
    sycl::queue Q(sycl::default_selector_v);
    std::cout << "Running on: "
              << Q.get_device().get_info<sycl::info::device::name>()
              << std::endl;
    

    std::cout << "oneMKL DPC++ GEMM benchmark\n"
              << "---------------------------\n"
              << "Device:                  " << Q.get_device().get_info<sycl::info::device::name>()                          << std::endl
              << "Core/EU count:           " << Q.get_device().get_info<sycl::info::device::max_compute_units>()             << std::endl
              << "Maximum clock frequency: " << Q.get_device().get_info<sycl::info::device::max_clock_frequency>() << " MHz" << std::endl;
    
    Matrix A(2,6), B(2,6), C(6,6), D(6,2); 
    
    // oneAPIMatrix B(Q, 4, 4);
    
    A.randMatrix(1,3);
    B.randMatrix(1,3);
    // C.resize(6,6);
    A.display();
    B.display();

    // naive::axpy(A, B, 1.0, 1.0);
    naive::gemm<true, false>(A,B, C, 1.0, 0.0);
    naive::gemm<false, true>(C, B, D, 1.0, 0.0);
    C.display();

    C.transposeInPlace();
    // C.display();
    naive::inverse(C, D);

    
    // Matrix F; F.fill_diag({1,4,3,35});

    D.mem_copy_to(A);
    A.display();
     
    Matrix in(1,3);
    in.randMatrix(1,3);
    std::cout <<"in\n";
    in.display();
    A.zeros();
    A.fill_matrix(in, 2, 2);
    A.display();

       





    

    return 0; 
}