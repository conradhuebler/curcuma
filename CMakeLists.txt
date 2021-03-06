cmake_minimum_required(VERSION 3.15)
project(curcuma)



option (COMPILE_XTB
    "Compile XTB and use those libs" OFF)

option (LINK_XTB
    "Link to prebuilt XTB library" OFF)

option (C17
    "Use C++ 17 Standard" ON)

add_subdirectory(${PROJECT_SOURCE_DIR}/external/fmt EXCLUDE_FROM_ALL)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/fmt/include ${PROJECT_SOURCE_DIR}/external/eigen ${PROJECT_SOURCE_DIR}/external/LBFGSpp/include ${PROJECT_SOURCE_DIR}/external/CppNumericalSolvers/include  ${CMAKE_CURRENT_BINARY_DIR}/src/)
include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest tag within the working branch
execute_process(
  COMMAND git describe --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_TAG
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(UNIX) # dont work on windows systems, compilation fails unfortranatly ...
    execute_process(
        COMMAND git log -1 --format=%aD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif(UNIX)

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")
add_definitions("-DGIT_COMMIT_DATE=${GIT_COMMIT_DATE}")
add_definitions("-DGIT_COMMIT_TAG=${GIT_COMMIT_TAG}")

if(CMAKE_COMPILER_IS_GNUCXX)
    set (_UNIX  ON)


endif()
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    set (openMP ON)
endif()

IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -pedantic \
                                           -Wextra -Wcast-align -Wcast-qual  -Wchar-subscripts  \
                                           -Wcomment -Wdisabled-optimization \
                                           -Wformat  -Wformat=2 -Wformat-nonliteral -Wformat-security\
                                           -Wformat-y2k  -Wimport  -Winit-self  -Winline -Winvalid-pch\
                                           -Wunsafe-loop-optimizations  -Wmissing-braces\
                                           -Wmissing-field-initializers -Wmissing-format-attribute  \
                                           -Wmissing-include-dirs -Wmissing-noreturn -Wpacked  -Wparentheses\
                                           -Wpointer-arith -Wredundant-decls -Wsequence-point\
                                           -Wsign-compare  -Wstack-protector -Wstrict-aliasing\
                                           -Wstrict-aliasing=2 -Wswitch -Wsuggest-override\
                                           -Wtrigraphs  -Wuninitialized -Wunknown-pragmas  -Wunreachable-code\
                                           -Wunused -Wunused-function  -Wunused-label  -Wunused-parameter\
                                           -Wunused-value  -Wunused-variable  -Wvariadic-macros\
                                           -Wvolatile-register-var  -Wwrite-strings -Wdeprecated-declarations\
                                           -Wno-error=unused-local-typedefs -Wno-error=enum-compare -Wno-narrowing -Werror=return-type -g -rdynamic")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


if(COMPILE_XTB)
    add_subdirectory(external/xtb)
    set_target_properties(lib-xtb-shared PROPERTIES
      IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/libxtb.so"
    )
    set(USE_XTB on)
endif()

if(LINK_XTB)
    if(COMPILE_XTB)
        message(WARNING "Actually, I will link curcuma to the currently building xtb not a pre-built one")
    else()
        if(DEFINED XTB_DIR)
            # Check does not work, strange
            #if(EXISTS '${XTB_DIR}/libxtb.so.6')
            #    message("XTB found :-)")
            #else()
            #    message("Not XTB at given location")
            #endif()
        else()
            message(FATAL_ERROR "Please set XTB_DIR as variable")
        endif()

    endif()
    set(USE_XTB on)
endif()




if(NOT EXISTS src/json.hpp)
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.9.1/json.hpp src/json.hpp
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=7804b38146921d03374549c9e2a5e3acda097814c43caf2b96a0278e58df26e0
       TLS_VERIFY ON)
endif()


configure_file (
  "${PROJECT_SOURCE_DIR}/src/global_config.h.in"
  "${PROJECT_BINARY_DIR}/src/global_config.h"
)

configure_file(
  "${CMAKE_SOURCE_DIR}/src/version.h.in"
  "${CMAKE_BINARY_DIR}/src/version.h"
)

set(curcuma_core_SRC
        src/capabilities/analysenciplot.cpp
        src/capabilities/curcumamethod.cpp
        src/capabilities/curcumaopt.cpp
        src/capabilities/confscan.cpp
        src/capabilities/confstat.cpp
        src/capabilities/docking.cpp
        src/capabilities/analysenciplot.cpp
        src/capabilities/nebdocking.cpp
        src/capabilities/pairmapper.cpp
        src/capabilities/rmsd.cpp
        src/capabilities/rmsdtraj.cpp
        src/capabilities/simplemd.cpp
        src/core/molecule.cpp
        src/core/xtbinterface.cpp
        src/core/pseudoff.cpp
        src/tools/geometry.h
        src/tools/general.h
        )
    add_library(curcuma_core  ${curcuma_core_SRC})

    add_executable(curcuma
            src/main.cpp)

if(C17)
    set_property(TARGET curcuma_core PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma_core PROPERTY CXX_STANDARD 14)
endif(C17)

if(C17)
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)
endif(C17)

set(target_dir "bin")
set_property(TARGET curcuma PROPERTY LINK_FLAGS_RELEASE -s)
set_property(TARGET curcuma PROPERTY LINK_FLAGS_RELWITHDEBINFO -s)

if(USE_XTB)
    if(COMPILE_XTB)
        target_link_libraries(curcuma_core lib-xtb-shared pthread dl fmt::fmt-header-only)
    else()
        target_link_libraries(curcuma_core ${XTB_DIR}/libxtb.so pthread dl fmt::fmt-header-only)
    endif()
else()
    target_link_libraries(curcuma_core pthread dl fmt::fmt-header-only)
endif()

    target_link_libraries(curcuma fmt::fmt-header-only)
    target_link_libraries(curcuma curcuma_core)

add_executable(curcuma_helper
        src/helpers/main.cpp)

    target_link_libraries(curcuma_helper fmt::fmt-header-only)
    target_link_libraries(curcuma_helper curcuma_core)
    
set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELEASE -s)
set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELWITHDEBINFO -s)

if(C17)
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 14)
endif(C17)

if(USE_XTB)
    add_executable(xtb_helper
        src/helpers/xtb_helper.cpp
        )
    if(C17)
        set_property(TARGET xtb_helper PROPERTY CXX_STANDARD 17)
    else()
        set_property(TARGET xtb_helper PROPERTY CXX_STANDARD 14)
    endif(C17)

    if(COMPILE_XTB)
        target_link_libraries(xtb_helper lib-xtb-shared pthread dl)
    else()
        target_link_libraries(xtb_helper ${XTB_DIR}/libxtb.so pthread dl)
    endif()
endif()

add_subdirectory(test_cases)

install(TARGETS curcuma RUNTIME DESTINATION bin)
install(TARGETS curcuma_helper RUNTIME DESTINATION bin)
