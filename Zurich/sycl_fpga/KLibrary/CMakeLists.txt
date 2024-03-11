
set(KLIBRARY_NAME KLibrary)

# Configuração para encontrar os arquivos de cabeçalho
set(LIBRARY_SOURCE
    "${CMAKE_CURRENT_LIST_DIR}/src/err_handler.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/Euler.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/file_reader.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/Initialization.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/kalman_filter.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/libmat.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/src/Quaternions.cpp"
)

include_directories(include)

add_library(${KLIBRARY_NAME} ${LIBRARY_SOURCE})
