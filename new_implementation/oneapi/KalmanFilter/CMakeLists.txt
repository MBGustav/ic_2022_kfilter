
set(KF_LIB KalmanFilter)

set(KF_SOURCE
    "KalmanFilter.cpp"
    "include/Quaternion.hpp"
    "include/KalmanFilter.hpp"
)
include_directories(include)

add_library(${KF_LIB} STATIC ${KF_SOURCE})
