set(MTX_LIB Matrix)

set(MTX_SOURCE
    "Matrix.cpp"
    "include/Linalg.hpp"
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(${MTX_LIB} STATIC ${MTX_SOURCE})
target_include_directories(${MTX_LIB} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
