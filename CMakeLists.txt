cmake_minimum_required(VERSION 3.18)
project(curcuma)

option (USE_XTB
    "Compile XTB and use as library" OFF)

option (USE_TBLITE
    "Compile TBLITE and use as library" OFF)

option (USE_D3
    "Compile cpp-D4 and use as library" OFF)

option (USE_D4
    "Compile cpp-D4 and use as library" OFF)

option (C17
    "Use C++ 17 Standard" ON)

option (Helper
    "Helper apps for interfaces" OFF)

option (WriteMoreInfo
    "Write statistic files with more info" OFF)

add_subdirectory(${PROJECT_SOURCE_DIR}/external/fmt EXCLUDE_FROM_ALL)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#if(NOT DEFINED ${CMAKE_BUILD_TYPE})
#    set(CMAKE_BUILD_TYPE RelWithDebInfo)
#endif()



include_directories(SYSTEM  ${PROJECT_SOURCE_DIR}/external/fmt/include
                            ${PROJECT_SOURCE_DIR}/external/LBFGSpp/include
                            ${PROJECT_SOURCE_DIR}/external/ripser
                            ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0
                            ${PROJECT_SOURCE_DIR}/external/cpp-d4/include
                            ${PROJECT_SOURCE_DIR}/external/)
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

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
                                           -Wunsafe-loop-optimizations  -Wno-missing-braces\
                                           -Wmissing-field-initializers -Wmissing-format-attribute  \
                                           -Wmissing-include-dirs -Wmissing-noreturn -Wpacked  -Wparentheses\
                                           -Wpointer-arith -Wno-redundant-decls -Wsequence-point\
                                           -Wsign-compare  -Wstack-protector -Wstrict-aliasing\
                                           -Wstrict-aliasing=2 -Wswitch -Wsuggest-override\
                                           -Wtrigraphs  -Wuninitialized -Wunknown-pragmas  -Wunreachable-code\
                                           -Wunused -Wunused-function  -Wunused-label  -Wunused-parameter\
                                           -Wunused-value  -Wunused-variable  -Wvariadic-macros\
                                           -Wvolatile-register-var  -Wwrite-strings -Wdeprecated-declarations\
                                           -Wno-error=unused-local-typedefs -Wno-error=enum-compare -Wno-narrowing -Werror=return-type -Wno-sign-compare -g")
set (GCC ON)
if(WIN32) # Check if we are on Windows
else()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
endif(WIN32)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/json.hpp)
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.9.1/json.hpp ${PROJECT_SOURCE_DIR}/external/json.hpp
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=7804b38146921d03374549c9e2a5e3acda097814c43caf2b96a0278e58df26e0
       TLS_VERIFY ON)
endif()


if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip)
  file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=1ccaabbfe870f60af3d6a519c53e09f3dcf630207321dffa553564a8e75c4fc8
       TLS_VERIFY ON)
endif()

file(ARCHIVE_EXTRACT INPUT ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip DESTINATION ${PROJECT_SOURCE_DIR}/external)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/global_config.h.in"
  "${PROJECT_BINARY_DIR}/src/global_config.h"
)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.h.in"
  "${PROJECT_BINARY_DIR}/src/version.h"
)

set(curcuma_core_SRC
        src/capabilities/persistentdiagram.cpp
        src/capabilities/analysenciplot.cpp
        src/capabilities/curcumamethod.cpp
        src/capabilities/curcumaopt.cpp
        src/capabilities/confscan.cpp
        src/capabilities/confsearch.cpp
        src/capabilities/confstat.cpp
        src/capabilities/docking.cpp
        src/capabilities/analysenciplot.cpp
        src/capabilities/nebdocking.cpp
        src/capabilities/pairmapper.cpp
        src/capabilities/rmsd.cpp
        src/capabilities/rmsdtraj.cpp
        src/capabilities/simplemd.cpp
        src/core/hessian.cpp
        src/core/energycalculator.cpp
        src/core/molecule.cpp
        #src/core/pseudoff.cpp
        src/core/eigen_uff.cpp
        src/tools/formats.h
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
    add_subdirectory(external/xtb)

    if(NOT TARGET mctc-lib)
      set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
      message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
      include(FetchContent)
      FetchContent_Declare(
        "mctc-lib"
        GIT_REPOSITORY "${mctc-lib-url}"
        GIT_TAG "HEAD"
      )
      FetchContent_MakeAvailable("mctc-lib")
    endif()
    include_directories(SYSTEM ${PROJECT_BINARY_DIR}/_deps/tblite-src/include/)

    if(HELPERS)
	add_executable(xtb_helper
        	src/helpers/xtb_helper.cpp
            )

    if(C17)
        set_property(TARGET xtb_helper PROPERTY CXX_STANDARD 17)
    else()
        set_property(TARGET xtb_helper PROPERTY CXX_STANDARD 14)
    endif(C17)

    target_link_libraries(xtb_helper lib-xtb-shared pthread dl)
    endif()
    set(curcuma_xtb_SRC
        src/core/xtbinterface.cpp
    )
    add_library(curcuma_xtb  ${curcuma_xtb_SRC})

    target_link_libraries(curcuma_xtb PUBLIC "mctc-lib" lib-xtb-shared gfortran)
    target_link_libraries(curcuma_core curcuma_xtb)

endif()

if(USE_TBLITE)

    if(NOT TARGET mctc-lib)
      set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
      message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
      include(FetchContent)
      FetchContent_Declare(
        "mctc-lib"
        GIT_REPOSITORY "${mctc-lib-url}"
        GIT_TAG "HEAD"
      )
      FetchContent_MakeAvailable("mctc-lib")
    endif()
    if(NOT TARGET tblite-lib)
    add_subdirectory(external/tblite)
    include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/tblite/include)
    elseif()
 #       include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/tblite/include)
    endif()
    if(HELPERS)
    add_executable(tblite_helper
           src/helpers/tblite_helper.cpp
    )
    if(C17)
        set_property(TARGET tblite_helper PROPERTY CXX_STANDARD 17)
    else()
        set_property(TARGET tblite_helper PROPERTY CXX_STANDARD 14)
    endif(C17)
    target_link_libraries(tblite_helper PUBLIC tblite-lib gfortran pthread)
    endif()
    set(curcuma_tblite_SRC
        src/core/tbliteinterface.cpp
    )
    add_library(curcuma_tblite  ${curcuma_tblite_SRC})

    target_link_libraries(curcuma_tblite PUBLIC "mctc-lib" tblite-lib gfortran)
    target_link_libraries(curcuma_core curcuma_tblite)
endif()

if(USE_D3)
    if(NOT TARGET s-dftd3)
        add_subdirectory(external/simple-dftd3)
    endif()
    if(HELPERS)
        add_executable(dftd3_helper
           src/helpers/dftd3_helper.cpp
    )
    target_link_libraries(dftd3_helper PUBLIC s-dftd3 gfortran pthread) 
    endif()
    set(curcuma_D3_SRC
           src/core/dftd3interface.cpp
    )
    add_library(curcuma_D3  ${curcuma_D3_SRC})
    target_link_libraries(curcuma_D3 s-dftd3)
    target_link_libraries(curcuma_core curcuma_D3 gfortran pthread)
endif()

if(USE_D4)
    if(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)
  	find_package(LAPACKE REQUIRED)
	find_package(CBLAS CONFIG REQUIRED)
    else()
	include_directories(${INCLUDE_DIR})
    endif(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)

    add_subdirectory(external/cpp-d4)
    if(HELPERS)
    add_executable(dftd4_helper
        src/helpers/dftd4_helper.cpp
        )
    if(C17)
        set_property(TARGET dftd4_helper PROPERTY CXX_STANDARD 17)
    else()
        set_property(TARGET dftd4_helper PROPERTY CXX_STANDARD 14)
    endif(C17)
    if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    	link_directories(${LIBS_DIR})
        target_link_libraries(dftd4_helper PRIVATE libcpp_d4 cblas lapacke pthread)
    else()
    target_link_libraries(dftd4_helper PRIVATE libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} pthread)
    endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    endif()

    set(curcuma_d4_SRC
        src/core/dftd4interface.cpp
    )
    add_library(curcuma_d4  ${curcuma_d4_SRC})
    if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 cblas lapacke )
    else()
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} )
    endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    target_link_libraries(curcuma_core curcuma_d4)
endif()

#target_link_libraries(curcuma_core pthread fmt::fmt-header-only finitediff::finitediff)
target_link_libraries(curcuma_core pthread fmt::fmt-header-only )

if(WIN32) # Check if we are on Windows
else()
     target_link_libraries(curcuma_core dl )
endif(WIN32)

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

add_subdirectory(test_cases)
include(CTest)
add_test(NAME AAAbGal_mtemplate COMMAND AAAbGal mtemplate WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_mhybrid COMMAND AAAbGal mhybrid WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_template COMMAND AAAbGal template WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_hybrid COMMAND AAAbGal hybrid WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_incremental COMMAND AAAbGal incr WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

set_tests_properties(AAAbGal_incremental PROPERTIES TIMEOUT 300)


install(TARGETS curcuma RUNTIME DESTINATION bin)
install(TARGETS curcuma_helper RUNTIME DESTINATION bin)
