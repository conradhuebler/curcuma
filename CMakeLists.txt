cmake_minimum_required(VERSION 3.1)
project(curcuma)


include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/external/eigen)

add_executable(curcuma
        src/capabilities/confscan.cpp
        src/capabilities/docking.cpp
        src/capabilities/rmsd.cpp
        src/capabilities/nebdocking.cpp
        src/core/molecule.cpp
        src/core/pseudoff.cpp
        src/tools/geometry.h
        src/tools/general.h
        src/main.cpp)
set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)
install(TARGETS curcuma RUNTIME DESTINATION bin)
