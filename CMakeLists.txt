cmake_minimum_required(VERSION 3.1)
project(curcuma)


add_subdirectory(external/xtb)

include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_SOURCE_DIR}/external/eigen)

set_target_properties(xtb_shared PROPERTIES
  IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libxtb_shared.so"
)


add_executable(curcuma
        src/capabilities/optimiser/LBFGSInterface.h
        src/capabilities/confscan.cpp
        src/capabilities/docking.cpp
        src/capabilities/rmsd.cpp
        src/capabilities/nebdocking.cpp
        src/capabilities/simplemd.cpp
        src/core/molecule.cpp
        src/core/xtbinterface.cpp
        src/core/pseudoff.cpp
        src/tools/geometry.h
        src/tools/general.h
        src/main.cpp)
set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)


target_link_libraries(curcuma xtb_shared)

install(TARGETS curcuma RUNTIME DESTINATION bin)
