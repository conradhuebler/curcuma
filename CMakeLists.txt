cmake_minimum_required(VERSION 3.18)
project(curcuma)

option (USE_XTB
    "Compile XTB and use as library" OFF)

option (USE_TBLITE
    "Compile TBLITE and use as library" ON)
 
option(USE_ULYSSES
    "Compile Ulysses and use as library" ON)

option (USE_D3
    "Compile cpp-D4 and use as library" ON)

option (USE_D4
    "Compile cpp-D4 and use as library" OFF)

option (C17
    "Use C++ 17 Standard" ON)

option (Helper
    "Helper apps for interfaces" OFF)

option (WriteMoreInfo
    "Write statistic files with more info" OFF)

option (USE_Plumed
    "Enable interface to plumed for metadynamics")

option(USE_AVX2
    "Enable AVX2 instructions" OFF)

option(USE_BLAS
    "Enable BLAS" OFF)

option(USE_MKL
    "Enable MKL" OFF)

option(CURCUMA_DEBUG
    "Enable debug output and timing" OFF)

option(CURCUMA_COLOR_DEFAULT
    "Enable colors by default" ON)

set(CURCUMA_DEBUG_LEVEL 2 CACHE STRING "Debug verbosity level 0-3")


option(USE_PCH
    "Enable Precompiled Headers for faster compilation" ON)

if(USE_MKL)
    if(NOT DEFINED MKL_ROOT)
        set(MKL_ROOT $ENV{MKL_ROOT})
    endif()
    message(STATUS "MKL_ROOT set to ${MKL_ROOT}")
    if( MKL_ROOT STREQUAL "")
        message(FATAL_ERROR "MKL_ROOT not set")
        message(FATAL_ERROR "Please set MKL_ROOT to the root directory of your MKL installation")
        message(FATAL_ERROR "For example: export MKL_ROOT=/opt/intel/oneapi/mkl/latest")
        message(FATAL_ERROR "Or use -DMKL_ROOT=/opt/intel/oneapi/mkl/latest")
    endif()
    set(MKLROOT ${MKL_ROOT})

endif()



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
                            ${PROJECT_SOURCE_DIR}/external/
                            ${PROJECT_SOURCE_DIR}/external/ulysses-main/src/
                            )


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

if(USE_AVX2)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

endif()

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/json.hpp)
  file(DOWNLOAD https://github.com/nlohmann/json/releases/download/v3.9.1/json.hpp ${PROJECT_SOURCE_DIR}/external/json.hpp
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=7804b38146921d03374549c9e2a5e3acda097814c43caf2b96a0278e58df26e0
       TLS_VERIFY ON)
endif()

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/stb_image_write.h)
  file(DOWNLOAD https://raw.githubusercontent.com/syoyo/tinygltf/refs/heads/release/stb_image_write.h ${PROJECT_SOURCE_DIR}/external/stb_image_write.h
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=bdbaa2c548ac9bde2913b6b46c4d315a04d2624c81ee1553390eb8fb2de2835d
       TLS_VERIFY ON)
endif()

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip)
  file(DOWNLOAD https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip
       TIMEOUT 60  # seconds
       EXPECTED_HASH SHA256=eba3f3d414d2f8cba2919c78ec6daab08fc71ba2ba4ae502b7e5d4d99fc02cda
       TLS_VERIFY ON)
endif()
file(ARCHIVE_EXTRACT INPUT ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0.zip DESTINATION ${PROJECT_SOURCE_DIR}/external)

if(NOT EXISTS ${PROJECT_SOURCE_DIR}/external/ulysses.zip)
  file(DOWNLOAD https://gitlab.com/siriius/ulysses/-/archive/main/ulysses-main.zip  ${PROJECT_SOURCE_DIR}/external/ulysses.zip
         TIMEOUT 60  # seconds
         EXPECTED_HASH SHA256=5f6f39e6fb12538ab30f9aef94a76d8e283151b30f1b6eef80f80b5bed46477a
         )
endif()

file(ARCHIVE_EXTRACT INPUT ${PROJECT_SOURCE_DIR}/external/ulysses.zip DESTINATION ${PROJECT_SOURCE_DIR}/external)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/global_config.h.in"
  "${PROJECT_BINARY_DIR}/src/global_config.h"
)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.h.in"
  "${PROJECT_BINARY_DIR}/src/version.h"
)

set(curcuma_core_SRC
        src/capabilities/optimiser/lbfgs.cpp
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
        src/capabilities/rmsd/rmsd_strategies.cpp
        src/capabilities/rmsd/rmsd_costmatrix.cpp
        src/capabilities/rmsd/rmsd_assignment.cpp
        src/capabilities/rmsdtraj.cpp
        src/capabilities/simplemd.cpp
        src/capabilities/hessian.cpp
        src/capabilities/qmdfffit.cpp
        src/core/qm_methods/eht.cpp
        src/core/qm_methods/eht_parameters.cpp
        src/core/qm_methods/gfnff.cpp
        src/core/qm_methods/gfnff_advanced.cpp
        src/core/qm_methods/qm_driver.cpp
        src/core/qm_methods/ulyssesinterface.cpp
        src/core/energy_calculators/method_factory.cpp
        src/core/energy_calculators/qm_methods/eht_method.cpp
        src/core/energy_calculators/qm_methods/gfnff_method.cpp
        src/core/energy_calculators/qm_methods/xtb_method.cpp
        src/core/energy_calculators/qm_methods/tblite_method.cpp
        src/core/energy_calculators/qm_methods/ulysses_method.cpp
        src/core/energy_calculators/qm_methods/dispersion_method.cpp
        src/core/energy_calculators/ff_methods/forcefield_method.cpp
        src/core/energycalculator.cpp
        src/core/molecule.cpp
        src/core/xyz_comment_parser.cpp
        src/core/fileiterator.cpp
        src/core/forcefieldthread.cpp
        src/core/forcefield.cpp
        src/core/forcefieldfunctions.h
        src/core/forcefieldgenerator.cpp
        src/core/curcuma_logger.cpp
        src/tools/formats.h
        src/tools/geometry.h
        src/tools/general.h
        src/core/orcainterface.h
        src/core/orcainterface.cpp
)

add_library(curcuma_core  ${curcuma_core_SRC})

# Enable precompiled headers for faster compilation
if(USE_PCH)
    target_precompile_headers(curcuma_core PRIVATE src/pch.h)
    message(STATUS "Precompiled headers enabled for curcuma_core")
    
    # Additional PCH configuration for Ulysses-heavy builds
    if(USE_ULYSSES)
        message(STATUS "Ulysses PCH module enabled - significant compilation speedup expected")
    endif()
    
    # Additional PCH configuration for external libraries
    if(USE_D3 OR USE_D4)
        message(STATUS "DFT dispersion correction PCH modules enabled")
    endif()
endif()

    set(km_SRC
    src/capabilities/c_code/hungarian.c
    src/capabilities/c_code/interface.c
)
add_library(km ${km_SRC})
if(USE_ULYSSES)

set(ulysess_SRC
    src/core/qm_methods/interface/ulysses.cpp
)
add_executable(ulysses_helper
src/helpers/ulysses_helper.cpp
)
else()
set(ulysess_SRC
    src/core/qm_methods/eht.cpp
)
endif(USE_ULYSSES)
add_library(ulysses_lib ${ulysess_SRC})


if(USE_BLAS)
    target_link_libraries(ulysses_lib blas)
endif()
if(USE_MKL)
    target_link_libraries(ulysses_lib -m64 -I${MKLROOT}/include  -Wl,--no-as-needed -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core)

endif()
if(USE_ULYSSES)

target_link_libraries(ulysses_helper PUBLIC ulysses_lib)
endif()



target_link_libraries(curcuma_core km ulysses_lib)

add_executable(curcuma
        src/main.cpp
)

# Enable precompiled headers for main executable
if(USE_PCH)
    target_precompile_headers(curcuma REUSE_FROM curcuma_core)
endif()

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
        src/core/qm_methods/xtbinterface.cpp
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
        src/core/qm_methods/tbliteinterface.cpp
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
           src/core/qm_methods/dftd3interface.cpp
    )
    add_library(curcuma_D3  ${curcuma_D3_SRC})
    target_link_libraries(curcuma_D3 s-dftd3)
    target_link_libraries(curcuma_core curcuma_D3 gfortran pthread)
endif()

if(USE_D4)
    if(NOT DEFINED D4LIBS_DIR AND NOT DEFINED D4INCLUDE_DIR)
        #find_package(LAPACKE REQUIRED)
        #find_package(CBLAS CONFIG REQUIRED)
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
        src/core/qm_methods/dftd4interface.cpp
    )
    add_library(curcuma_d4  ${curcuma_d4_SRC})
    if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 cblas lapacke )
    else()
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} )
    endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    target_link_libraries(curcuma_core curcuma_d4)
endif()

if(USE_Plumed)
   target_link_libraries(curcuma_core ${PROJECT_SOURCE_DIR}/external/plumed2/src/lib/libplumed.so)
endif()

target_link_libraries(curcuma_core pthread fmt::fmt-header-only )
if(USE_BLAS)
    target_link_libraries(curcuma_core blas)
endif()
if(USE_MKL)
    target_link_libraries(curcuma_core -m64 -I${MKLROOT}/include  -Wl,--no-as-needed -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core)

endif()
if(WIN32) # Check if we are on Windows
else()
     target_link_libraries(curcuma_core dl )
endif(WIN32)

target_link_libraries(curcuma fmt::fmt-header-only)

target_link_libraries(curcuma curcuma_core)
add_executable(curcuma_helper
        src/helpers/main.cpp)
        
add_executable(cli_helper
        src/helpers/cli_test.cpp)

add_executable(parallel_scf
        src/helpers/parallel_scf.cpp)



target_link_libraries(curcuma_helper fmt::fmt-header-only)
target_link_libraries(curcuma_helper curcuma_core)

# Enable precompiled headers for helper executables
if(USE_PCH)
    target_precompile_headers(curcuma_helper REUSE_FROM curcuma_core)
    target_precompile_headers(cli_helper REUSE_FROM curcuma_core)
    target_precompile_headers(parallel_scf REUSE_FROM curcuma_core)
endif()

set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELEASE -s)
set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELWITHDEBINFO -s)

if(C17)
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 14)
endif(C17)

add_subdirectory(test_cases)
include(CTest)

target_sources(Continuous
  PRIVATE
    src/core/qmdff_par.h
)
add_test(NAME molecule_comprehensive COMMAND test_molecule)

add_test(NAME AAAbGal_dtemplate COMMAND AAAbGal dtemplate WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_free COMMAND AAAbGal free WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_template COMMAND AAAbGal template WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME AAAbGal_subspace COMMAND AAAbGal subspace WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
#add_test(NAME AAAbGal_incremental COMMAND AAAbGal incr WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

#set_tests_properties(AAAbGal_incremental PROPERTIES TIMEOUT 300)

add_test(NAME confscan_subspace COMMAND confscan_test subspace WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_free COMMAND confscan_test free WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_dtemplate COMMAND confscan_test dtemplate WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_template COMMAND confscan_test template WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_molalign COMMAND confscan_test molalign WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX1 COMMAND confscan_test sLX1 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX2 COMMAND confscan_test sLX2 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX2Reset COMMAND confscan_test sLX2Reset WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX20 COMMAND confscan_test sLX20 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX20Reset COMMAND confscan_test sLX20Reset WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

add_test(NAME confstat COMMAND confstat_tests WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

# Claude Generated: Comprehensive test for Molecule class before refactoring

install(TARGETS curcuma RUNTIME DESTINATION bin)
install(TARGETS curcuma_helper RUNTIME DESTINATION bin)
