set(LIB_NAME oneAPIMatrix)

find_package(MKL CONFIG REQUIRED)
message(STATUS "${MKL_IMPORTED_TARGETS}")


set(SOURCE_FILE
    "oneAPIMatrix.cpp"
)
include_directories(include)

add_library(${LIB_NAME} STATIC ${SOURCE_FILE})

target_sources(${LIB_NAME} PRIVATE oneAPIMatrix.cpp)
target_compile_options(${LIB_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(${LIB_NAME} PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(${LIB_NAME} PUBLIC $<LINK_ONLY:MKL::MKL>)

