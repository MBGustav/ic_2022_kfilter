
set(CMAKE_CXX_COMPILER icpx)

cmake_minimum_required (VERSION 3.7)

project(Matmul CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


file(GLOB MLAB_FILES "${CMAKE_SOURCE_DIR}/../DadosMatLab/*.csv")
file(GLOB CSV_FILES  "${CMAKE_SOURCE_DIR}/../DadosZurich/*.csv")


file(COPY ${MLAB_FILES} 
DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/output"
)

file(COPY ${CSV_FILES} 
DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/DadosZurich"
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

add_subdirectory(KLibrary)
# intel sample
add_subdirectory (matmul_test)
# kalman testing
# add_subdirectory(kalman_filter)
