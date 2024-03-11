cmake_minimum_required (VERSION 3.7)

set(CMAKE_CXX_COMPILER icpx)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-g -O3 -fsycl -Wno-write-strings -w ")

project(test CXX)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include_directories(common)
include_directories(Matrix/include)
# include_directories(oneAPIMatrix/include)
include_directories(KalmanFilter/include)

add_subdirectory(Matrix)
# add_subdirectory(oneAPIMatrix)
add_subdirectory(KalmanFilter)

add_executable(main main.cpp)
# target_sources(main PRIVATE oneAPIMatrix/oneAPIMatrix.cpp)

target_link_libraries(main PRIVATE 
                    Matrix 
                    # oneAPIMatrix
                    KalmanFilter
)