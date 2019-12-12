cmake_minimum_required(VERSION 3.1)
project(curcuma)


include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/external/eigen)

add_executable(curcuma 
        src/tools/geometry.h
        src/capabilities/rmsd.cpp
        src/core/molecule.cpp
        src/main.cpp)
set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)
install(TARGETS curcuma RUNTIME DESTINATION bin)
