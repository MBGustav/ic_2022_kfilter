



set(SOURCES main.cpp)

# target_link_libraries(${SOURCE_MAIN} KLibrary)
# Create a executable

# Include folder with .hpp or .h (prototypes)
include_directories(${CMAKE_SOURCE_DIR}/KLibrary/include)

set(WIDTH 16)
set(COMMON 16)
set(TILE_A 8)
set(TILE_B 8)
set(COMPILE_FLAGS "-Wall -DWIDTH=${WIDTH} -DCOMMON=${COMMON} -DTILE_B=${TILE_B} -DTILE_A=${TILE_A}")

add_executable(main ${SOURCES})
target_link_libraries(main  KLibrary)
set_target_properties(main PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")

