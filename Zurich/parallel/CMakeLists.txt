cmake_minimum_required(VERSION 3.10)
project(SeuProjeto)

set(REPORT_DIRECTORY report)
set(SRCDIR include)
set(INCDIR include)

# set(CMAKE_CXX_COMPILER gcc) 
set(DATA_PRECISION 1 CACHE STRING "Data precision: 1 for double, 0 for float")
set(CMAKE_C_COMPILER icpx)
set(CMAKE_CXX_COMPILER icpx)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -D_PROFILING")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -xHost -qopenmp -Wno-writable-strings")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl -xHost -qopenmp -Wno-writable-strings")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_PROFILING -DDATA_PRECISION=${DATA_PRECISION}")

set(EXEC main)
set(BENCHMARK bench)

set(ONEAPI_ROOT opt/intel/oneapi)
set(APM ${ONEAPI_ROOT}/advisor/latest/pythonapi/advisor/)
set(MKL_COPTS "-DMKL_ILP64 -I${MKLROOT}/include")
set(MKL_LIBS "-L${MKLROOT}/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lsycl -lOpenCL -lpthread -lm -ldl")
set(DPCPP_OPTS "${MKL_COPTS} -fsycl-device-code-split=per_kernel ${MKL_LIBS} -fopenmp-targets=spir64")
set(MKL_PAR "${MKL_COPTS} -fsycl-device-code-split=per_kernel ${MKL_LIBS}")
set(MKL_SEQ "${MKL_COPTS} -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")

file(GLOB C_FILES "${SRCDIR}/*.c")
file(GLOB H_FILES "${INCDIR}/*.h")



add_executable(${EXEC} ${C_FILES} ${H_FILES} ${EXEC}.cpp)
target_compile_options(${EXEC} PRIVATE ${DPCPP_OPTS})

add_executable(UnitTest UnitTest.cpp ${C_FILES})
target_include_directories(UnitTest PRIVATE ${INCDIR})

add_executable(CheckResult CheckResult.c ${C_FILES})
target_compile_options(CheckResult PRIVATE ${CMAKE_C_FLAGS} ${CFLAGS})
target_include_directories(CheckResult PRIVATE ${INCDIR})

add_custom_target(OutputCompare COMMAND python3 plot_output.py)

add_executable(testbench testbench-omp.cpp ${C_FILES})
target_compile_options(testbench PRIVATE ${DPCPP_OPTS} -DDATA_PRECISION=${DATA_PRECISION})
target_include_directories(testbench PRIVATE 
                            ${INCDIR}
                            ${INCDIR}/kernels
                            ${BLAS_LIBRARIES}
                             )


# Adicione outros alvos ou ajustes conforme necessário

# Adicione a seção de geração de relatórios do gnu profiler conforme necessário
# add_custom_target(TimeExecution_gnu
#     COMMAND mkdir -p ${REPORT_DIRECTORY}
#     COMMAND gcc ${CMAKE_C_FLAGS} -D_PROFILING -pg -I ${INCDIR} ${C_FILES} -lm ${EXEC}.c -o ${EXEC}
#     COMMAND ./${EXEC}
#     COMMAND gprof ./${EXEC} gmon.out -b > ${REPORT_DIRECTORY}/report.txt
#     COMMAND gprof2dot ${REPORT_DIRECTORY}/report.txt > ${REPORT_DIRECTORY}/report.dot
#     COMMAND dot -Tpng -o ${REPORT_DIRECTORY}/report.png ${REPORT_DIRECTORY}/report.dot
#     COMMAND python3 graphic.py fun_time.out
#     COMMAND rm gmon.out
# )

# Alvo para limpar
add_custom_target(clean_all
    COMMAND rm -rf ${EXEC} test *.out UnitTest ${REPORT_DIRECTORY} ${BENCHMARK} ${BENCHMARK}_par ${BENCHMARK}_seq
)