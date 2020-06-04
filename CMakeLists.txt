cmake_minimum_required(VERSION 3.15)
project(curcuma)



option (COMPILE_XTB
    "Compile XTB and use those libs" OFF)

option (LINK_XTB
    "Link to prebuilt XTB library" OFF)

option (C17
    "Use C++ 17 Standard" ON)

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/eigen ${CMAKE_CURRENT_BINARY_DIR}/src/)
include_directories(${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
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
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.7.3/json.hpp src/json.hpp
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA512=4ecbbdd2c5e88c897096670cfdaa7ec00483ac9ed6e8ac33be23b05f4da70f213e10c4b381f5e8799619a24a58f417f04dd442d1d59a2e0bfca3385007e620d5
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

add_executable(curcuma
        src/capabilities/optimiser/LBFGSInterface.h
        src/capabilities/curcumamethod.cpp
        src/capabilities/confscan.cpp
        src/capabilities/docking.cpp
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
        src/tools/parser.h
        src/main.cpp)
if(C17)
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)
endif(C17)

if(USE_XTB)
    if(COMPILE_XTB)
        target_link_libraries(curcuma lib-xtb-shared pthread dl)
    else()
        target_link_libraries(curcuma ${XTB_DIR}/libxtb.so pthread dl)
    endif()
else()
    target_link_libraries(curcuma  pthread dl)
endif()


add_executable(curcuma_helper
        src/helpers/main.cpp)

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

install(TARGETS curcuma RUNTIME DESTINATION bin)
install(TARGETS curcuma_helper RUNTIME DESTINATION bin)
