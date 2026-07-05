cmake_minimum_required(VERSION 3.18)
project(curcuma)

# Claude Generated 2026 - MSVC doesn't define M_PI (and friends) from <cmath>
# unless _USE_MATH_DEFINES is set before the first include; 20+ files across
# the codebase use M_PI directly (gcc/clang define it unconditionally), so
# set this project-wide instead of patching every translation unit.
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
    # Claude Generated 2026 - the codebase uses non-ASCII UTF-8 characters in a few
    # identifiers/comments (e.g. "double Δqraw[6];" in xtb_gradient.cpp). Without
    # /utf-8, MSVC falls back to the system codepage to interpret the source file,
    # misreading the multi-byte UTF-8 sequence as unrelated characters (observed:
    # "error C3872: 'U+201d' is not allowed in an identifier" for what is actually Δ).
    add_compile_options(/utf-8)
endif()

option (USE_XTB
    "Compile XTB and use as library" OFF)

option (USE_TBLITE
    "Compile TBLITE and use as library" OFF)

option (USE_GFNFF
    "Compile external GFN-FF and use as library" OFF)
 
option(USE_ULYSSES
    "Compile Ulysses and use as library" OFF)

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

option (USE_Plumed
    "Enable interface to plumed for metadynamics")

option(USE_MARCH_NATIVE
    "Enable -march=native (best performance, NOT portable)" OFF)

option(USE_AVX2
    "Enable AVX2 + FMA instructions (Haswell+, not portable)" OFF)

option(USE_AVX512
    "Enable AVX-512 instructions (Skylake-X+, not portable)" OFF)

option(USE_BLAS
    "Enable BLAS" OFF)

option(USE_MKL
    "Enable MKL" OFF)

# BLAS/LAPACK discovery for Eigen optimization
# Claude Generated (March 2026): Enable BLAS/LAPACK for native QM methods
if(USE_BLAS AND NOT USE_MKL)
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if(BLAS_FOUND)
        message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
    else()
        message(WARNING "BLAS requested via USE_BLAS but not found - using Eigen internal implementations")
    endif()
    if(LAPACK_FOUND)
        message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
    endif()
endif()

option(CURCUMA_DEBUG
    "Enable debug output and timing" OFF)

option(CURCUMA_COLOR_DEFAULT
    "Enable colors by default" ON)

set(CURCUMA_DEBUG_LEVEL 2 CACHE STRING "Debug verbosity level 0-3")


option(USE_PCH
    "Enable Precompiled Headers for faster compilation" ON)

option(USE_CUDA
    "Enable CUDA-accelerated GPU (gfnff + gfn1/gfn2 -gpu cuda, requires CUDA toolkit)" OFF)

# Claude Generated (2026-06): CUDA-accelerated native xTB (gfn1/gfn2 -gpu cuda).
# Defaults to USE_CUDA. Controls only which xTB GPU sources are added; the runtime
# code gates on USE_CUDA, the factory offers the GPU path when USE_CUDA_XTB is set.
option(USE_CUDA_XTB
    "Enable CUDA-accelerated native xTB SCF/integrals/gradient (gfn1/gfn2 -gpu cuda)" ${USE_CUDA})

# Claude Generated (June 2026): ROCm/HIP-accelerated GFN-FF + native xTB (-gpu rocm).
# Requires the ROCm toolkit (hipcc, hipBLAS, rocSOLVER). Separate HIP source tree
# under .../rocm/ mirrors the CUDA .../cuda/ kernels; runtime code gates on USE_ROCM
# (gfnff) / USE_ROCM_XTB (native xTB), exactly like the CUDA flags.
option(USE_ROCM
    "Enable ROCm/HIP-accelerated GPU (gfnff/gfn1/gfn2 -gpu rocm, requires ROCm toolkit)" OFF)
option(USE_ROCM_XTB
    "Enable ROCm/HIP-accelerated native xTB SCF/integrals/gradient (gfn1/gfn2 -gpu rocm)" ${USE_ROCM})
option(USE_ROCM_GFNFF
    "Enable ROCm/HIP-accelerated GFN-FF energy/gradient (gfnff -gpu rocm)" ${USE_ROCM})

# Claude Generated (June 2026): Vulkan-compute-accelerated GFN-FF + native xTB
# (-gpu vulkan). Hand-written SPIR-V compute shaders (no BLAS/LAPACK on Vulkan).
# Requires the Vulkan SDK (loader + headers + glslc/glslangValidator). Runtime code
# gates on USE_VULKAN (gfnff) / USE_VULKAN_XTB (native xTB).
option(USE_VULKAN
    "Enable Vulkan-compute-accelerated GPU (gfnff/gfn1/gfn2 -gpu vulkan, requires Vulkan SDK)" OFF)
option(USE_VULKAN_XTB
    "Enable Vulkan-compute-accelerated native xTB SCF/integrals/gradient (gfn1/gfn2 -gpu vulkan)" ${USE_VULKAN})

# ---------------------------------------------------------------------------
# Unified GPU backend switches (Jun 2026)
# ---------------------------------------------------------------------------
# Dissolve the per-method sub-flags: the umbrella USE_CUDA / USE_ROCM / USE_VULKAN is the
# single switch — turning one on enables BOTH the native-xTB (gfn1/gfn2) AND the GFN-FF GPU
# path for that backend. The per-method sub-flags are kept only for the source-gating logic
# below (and explicit/legacy `-DUSE_*_XTB=...` invocation) and are forced to follow their
# umbrella here so the two can never diverge — previously a sub-flag only DEFAULTED to the
# umbrella on the FIRST configure and silently stuck afterwards. (Plain set() shadows the
# cache for this configure and propagates to the GPU blocks + test subdirectory below.)
if(USE_CUDA)
    set(USE_CUDA_XTB ON)    # CUDA GFN-FF is gated directly by USE_CUDA (no separate flag)
endif()
if(USE_ROCM)
    set(USE_ROCM_XTB   ON)
    set(USE_ROCM_GFNFF ON)
endif()
if(USE_VULKAN)
    set(USE_VULKAN_XTB ON)  # + USE_VULKAN_GFNFF once the Vulkan GFN-FF stack (V-AP6) lands
endif()

# WP-D Stage B (May 2026): hand-rolled AVX2 exp(-x^2) for GFN-FF dcn step 3 hot loop.
# Opt-in. Default OFF keeps the build bit-identical to pre-WP. Requires AVX2 (set
# USE_AVX2 or USE_MARCH_NATIVE alongside); CMake errors out if AVX2 cannot be enabled.
option(USE_GFNFF_FAST_EXP
    "Enable SIMD-vectorized exp() in GFN-FF dcn hot loop (WP-D Stage B, requires AVX2)" OFF)
option(GFNFF_CN_DCN_FUSION
    "Enable fused CN+DCN single-pass (WP-D Stage D, no SIMD dependency)" OFF)
if(GFNFF_CN_DCN_FUSION)
    add_compile_definitions(GFNFF_CN_DCN_FUSION=1)
endif()

if(USE_MKL)
    # Support both MKL_ROOT and MKLROOT (environment and CMake variables)
    if(NOT MKL_ROOT)
        if(DEFINED ENV{MKL_ROOT})
            set(MKL_ROOT $ENV{MKL_ROOT})
        elseif(DEFINED ENV{MKLROOT})
            set(MKL_ROOT $ENV{MKLROOT})
        endif()
    endif()
    if(NOT MKL_ROOT AND DEFINED MKLROOT)
        set(MKL_ROOT ${MKLROOT})
    endif()

    message(STATUS "MKL_ROOT set to ${MKL_ROOT}")
    if(NOT MKL_ROOT OR MKL_ROOT STREQUAL "")
        message(FATAL_ERROR "MKL_ROOT not set. Please set MKL_ROOT or MKLROOT to the root directory of your MKL installation.\n"
                            "Examples:\n"
                            "  export MKL_ROOT=/opt/intel/oneapi/mkl/latest\n"
                            "  cmake -DMKL_ROOT=/opt/intel/oneapi/mkl/latest ..")
    endif()
    set(MKLROOT ${MKL_ROOT})
endif()



# External dependencies via FetchContent (replaces git submodules)
# SOURCE_DIR set to external/<name> to preserve existing include paths without any source changes
include(FetchContent)

FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG        676c2a107e9575663f79095d9b2c9d15168e1285
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/fmt
)
FetchContent_Declare(LBFGSpp
    GIT_REPOSITORY https://github.com/conradhuebler/LBFGSpp
    GIT_TAG        d16c89da2a490307791764aeb2778fa0c1458960
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/LBFGSpp
)
FetchContent_Declare(CxxThreadPool
    GIT_REPOSITORY https://github.com/conradhuebler/CxxThreadPool.git
    GIT_TAG        140625c509f245185ec486b5f5424b318c8ccef6
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/CxxThreadPool
)
# ripser: header-only, no CMakeLists.txt at root — use Populate to avoid add_subdirectory
FetchContent_Declare(ripser
    GIT_REPOSITORY https://github.com/conradhuebler/ripser
    GIT_TAG        dcb13c22db2017ef80e68bdf9bed6ec3b77ad70d
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/ripser
)
# cpp-d4: headers always needed, library only when USE_D4=ON
FetchContent_Declare(cpp-d4
    GIT_REPOSITORY https://github.com/conradhuebler/cpp-d4
    GIT_TAG        95882c7203e3b847528d0746a92b184f4cef4dc3
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/cpp-d4
)
# Conditional deps — declared here, FetchContent_MakeAvailable called in if(USE_*) blocks only
FetchContent_Declare(xtb
    GIT_REPOSITORY https://github.com/grimme-lab/xtb
    GIT_TAG        26b28010e805f7d1aeeef39813feb473e69cc4be
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/xtb
)
FetchContent_Declare(tblite
    GIT_REPOSITORY https://github.com/tblite/tblite
    GIT_TAG        9ca469a77cb9433dc6b0469ba5f6effc48a1bda8
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/tblite
)
FetchContent_Declare(gfnff
    GIT_REPOSITORY https://github.com/pprcht/gfnff.git
    GIT_TAG        0491df2ff1fe3cd80fc6a8630fea16f9b6840996
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/gfnff
)
FetchContent_Declare(simple-dftd3
    GIT_REPOSITORY https://github.com/dftd3/simple-dftd3
    GIT_TAG        3a50cbdf8566b3caa8760499ae21746ec9a5e535
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/simple-dftd3
)
FetchContent_Declare(plumed2
    GIT_REPOSITORY https://github.com/plumed/plumed2
    GIT_TAG        d4f7b550416d5aba2ca44ce6c2bca981282acea7
    SOURCE_DIR     ${CMAKE_SOURCE_DIR}/external/plumed2
)

# fmt: full MakeAvailable to get the fmt::fmt-header-only CMake target
FetchContent_MakeAvailable(fmt)

# Header-only libs: Populate only — no add_subdirectory needed (avoids old cmake_minimum_required issues)
# CMP0169 was introduced in CMake 3.30 — only set if available
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.30")
    cmake_policy(SET CMP0169 OLD)
endif()
foreach(_fc_dep LBFGSpp CxxThreadPool ripser cpp-d4)
    string(TOLOWER "${_fc_dep}" _fc_dep_lower)
    FetchContent_GetProperties(${_fc_dep} POPULATED ${_fc_dep_lower}_populated)
    if(NOT ${_fc_dep_lower}_populated)
        FetchContent_Populate(${_fc_dep})
    endif()
endforeach()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#if(NOT DEFINED ${CMAKE_BUILD_TYPE})
#    set(CMAKE_BUILD_TYPE RelWithDebInfo)
#endif()



include_directories(SYSTEM  ${CMAKE_SOURCE_DIR}/external/fmt/include
                            ${CMAKE_SOURCE_DIR}/external/LBFGSpp/include
                            ${CMAKE_SOURCE_DIR}/external/ripser
                            ${PROJECT_SOURCE_DIR}/external/eigen-3.4.0
                            ${CMAKE_SOURCE_DIR}/external/cpp-d4/include
                            ${PROJECT_SOURCE_DIR}/external/
                            ${CMAKE_SOURCE_DIR}/external/ulysses-main/src/
                            )


include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

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
# Claude Generated (March 2026): Remove comma from RFC 2822 date — nvcc splits -D values at commas
string(REPLACE "," " " GIT_COMMIT_DATE "${GIT_COMMIT_DATE}")
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
    # Claude Generated 2026 - OpenMP_CXX_FLAGS is only the compiler flag (e.g. Apple
    # Clang's "-Xclang -fopenmp"); it does NOT carry omp.h's include path. On Linux/
    # GCC omp.h ships with the compiler itself so this was never missed, but Apple
    # Clang has no OpenMP runtime at all -- Homebrew's libomp is keg-only, and without
    # this, any TU that includes omp.h (directly, or via Eigen's OpenMP integration)
    # fails with "fatal error: 'omp.h' file not found".
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    # Claude Generated 2026 - similarly, OpenMP_EXE_LINKER_FLAGS alone is not enough
    # to actually link the OpenMP runtime on Apple Clang: on GCC, -fopenmp at link
    # time makes the driver pull in libgomp itself, but Apple Clang's "-Xclang
    # -fopenmp" carries no such magic -- the actual "-lomp" comes from
    # OpenMP_CXX_LIBRARIES, which the old flags-only approach here never added.
    # Without it, every executable failed at link with "symbol(s) not found for
    # architecture arm64" (_omp_get_thread_num, _omp_set_num_threads, ...).
    link_libraries(${OpenMP_CXX_LIBRARIES})
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

# CPU feature flags — safe default: none (portable). Enable for target machine only.
include(CheckCXXCompilerFlag)
if(USE_MARCH_NATIVE)
    check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
        message(STATUS "CPU: -march=native enabled")
    else()
        message(WARNING "-march=native not supported by compiler")
    endif()
elseif(USE_AVX512)
    check_cxx_compiler_flag("-mavx512f" COMPILER_SUPPORTS_AVX512)
    if(COMPILER_SUPPORTS_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl -mavx512dq -mfma")
        message(STATUS "CPU: AVX-512 enabled")
    else()
        message(WARNING "AVX-512 not supported by compiler")
    endif()
elseif(USE_AVX2)
    check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
        message(STATUS "CPU: AVX2 + FMA enabled")
    else()
        message(WARNING "AVX2 not supported by compiler")
    endif()
endif()

# Eigen ABI alignment (Jun 2026): pin EIGEN_MAX_ALIGN_BYTES to ONE value derived from the
# active CPU-vectorization switch, applied IDENTICALLY to the g++ build (here) and the
# hipcc/nvcc-compiled GPU TUs below. A mismatch (e.g. g++ AVX-512 → 64 vs hipcc clang
# default → 16) gives Eigen objects a different allocator alignment AND a different sizeof
# for over-aligned fixed-size types across the TU boundary → aligned-load SIGSEGV / silent
# corruption (this was the GFN-FF gradient-norm crash on multi-fragment systems). Both sides
# reference ${CURCUMA_EIGEN_ALIGN}, so they can never diverge. 0 = leave Eigen's default (16)
# on both sides (no override) — only set it when a vectorization switch is actually on.
set(CURCUMA_EIGEN_ALIGN 0)
if(USE_MARCH_NATIVE OR USE_AVX512)
    set(CURCUMA_EIGEN_ALIGN 64)
elseif(USE_AVX2)
    set(CURCUMA_EIGEN_ALIGN 32)
endif()
if(NOT CURCUMA_EIGEN_ALIGN EQUAL 0)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_MAX_ALIGN_BYTES=${CURCUMA_EIGEN_ALIGN}")
    message(STATUS "Eigen: EIGEN_MAX_ALIGN_BYTES=${CURCUMA_EIGEN_ALIGN} (pinned for g++/GPU-TU ABI parity)")
endif()

# WP-D Stage B (May 2026): USE_GFNFF_FAST_EXP requires AVX2. Either USE_AVX2 or
# USE_MARCH_NATIVE (or USE_AVX512) must be set so __AVX2__ is defined.
if(USE_GFNFF_FAST_EXP)
    if(NOT (USE_AVX2 OR USE_AVX512 OR USE_MARCH_NATIVE))
        message(FATAL_ERROR "USE_GFNFF_FAST_EXP=ON requires AVX2: enable USE_AVX2, USE_AVX512, or USE_MARCH_NATIVE.")
    endif()
    add_compile_definitions(GFNFF_FAST_EXP=1)
    message(STATUS "GFN-FF: fast SIMD exp() enabled in dcn hot loop (WP-D Stage B)")
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

# Ulysses: fetched via FetchContent.
# UPDATE_DISCONNECTED prevents the "not a git repository" CI failure when the
# source tree was already populated as a plain directory (no .git present).
# GIT_SHALLOW reduces network traffic on first download.
FetchContent_Declare(ulysses
    GIT_REPOSITORY      https://gitlab.com/siriius/ulysses.git
    GIT_TAG             main
    GIT_SHALLOW         TRUE
    UPDATE_DISCONNECTED TRUE
    SOURCE_DIR          ${CMAKE_SOURCE_DIR}/external/ulysses-main
)
if(USE_ULYSSES)
    FetchContent_MakeAvailable(ulysses)
else()
    # Headers needed for the include path; only fetch if not already present.
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/ulysses-main/src")
        FetchContent_GetProperties(ulysses POPULATED ulysses_populated)
        if(NOT ulysses_populated)
            FetchContent_Populate(ulysses)
        endif()
    endif()
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/src/global_config.h.in"
  "${PROJECT_BINARY_DIR}/src/global_config.h"
)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/version.h.in"
  "${PROJECT_BINARY_DIR}/src/version.h"
)

set(curcuma_core_SRC
        # Optimization algorithms - Claude Generated (Educational-First Design)
        src/capabilities/optimisation/rf_solver.cpp
        src/capabilities/optimisation/lbfgs.cpp
        src/capabilities/optimizer_interface.cpp
        src/capabilities/optimizer_driver.cpp
        src/capabilities/optimizer_factory.cpp
        src/capabilities/ancopt_optimizer.cpp
        src/capabilities/lbfgspp_optimizer.cpp
        src/capabilities/native_optimizer_adapters.cpp
        src/capabilities/persistentdiagram.cpp
        src/capabilities/analysenciplot.cpp
        src/capabilities/curcumamethod.cpp
        # Native LBFGS optimizer - Claude Generated (disabled, superseded by native_optimizer_adapters)
        # src/capabilities/native_lbfgs_optimizer.cpp
        src/capabilities/confscan.cpp
        src/capabilities/confsearch.cpp
        src/capabilities/shared_bias_pool.cpp
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
        src/capabilities/md_diagnostics.cpp
        src/capabilities/analysis.cpp
        src/capabilities/analysis_output.cpp
        src/capabilities/handlers/general_properties_handler.cpp
        src/capabilities/handlers/scattering_handler.cpp
        src/capabilities/handlers/rdf_handler.cpp
        src/capabilities/handlers/topology_handler.cpp
        src/capabilities/trajectory_statistics.cpp
        src/capabilities/casino.cpp
        src/capabilities/scnp_parser.cpp
        src/capabilities/trajectoryanalysis.cpp
        src/capabilities/polymerbuild.cpp
        src/capabilities/tda_engine.cpp
        src/capabilities/hessian.cpp
        src/capabilities/qmdfffit.cpp
        src/core/energy_calculators/qm_methods/eht.cpp
        src/core/energy_calculators/qm_methods/eht_parameters.cpp
        src/core/energy_calculators/qm_methods/xtb_native.cpp
        src/core/energy_calculators/qm_methods/xtb_h0.cpp
        src/core/energy_calculators/qm_methods/xtb_coulomb.cpp
        src/core/energy_calculators/qm_methods/xtb_thirdorder.cpp
        src/core/energy_calculators/qm_methods/xtb_scf.cpp
        src/core/energy_calculators/qm_methods/xtb_fragment_scf.cpp
        src/core/energy_calculators/qm_methods/xtb_dc.cpp
        src/core/energy_calculators/qm_methods/xtb_sparse.cpp
        src/core/energy_calculators/qm_methods/native_eigensolver.cpp
        src/core/energy_calculators/qm_methods/xtb_multipole.cpp
        src/core/energy_calculators/qm_methods/xtb_response.cpp
        src/core/energy_calculators/qm_methods/xtb_gradient.cpp
        src/core/energy_calculators/qm_methods/nddo_params.cpp
        src/core/energy_calculators/ff_methods/gfnff_method.cpp
        src/core/energy_calculators/ff_methods/gfnff_advanced.cpp
        src/core/energy_calculators/ff_methods/gfnff_torsions.cpp
        src/core/energy_calculators/ff_methods/gfnff_inversions.cpp
        src/core/energy_calculators/ff_methods/ff_workspace.cpp
        src/core/energy_calculators/ff_methods/ff_workspace_gfnff.cpp
        src/core/energy_calculators/ff_methods/ff_workspace_uff.cpp
        src/core/energy_calculators/ff_methods/alpb_solvation.cpp
        src/core/energy_calculators/ff_methods/eeq_solver.cpp
        src/core/energy_calculators/ff_methods/huckel_solver.cpp
        src/core/energy_calculators/qm_methods/qm_driver.cpp
        src/core/energy_calculators/qm_methods/ulyssesinterface.cpp
        src/core/energy_calculators/method_factory.cpp
        src/core/energy_calculators/qm_methods/eht_method.cpp
        src/core/energy_calculators/qm_methods/native_xtb_method.cpp
        src/core/energy_calculators/qm_methods/nddo.cpp
        src/core/energy_calculators/qm_methods/nddo_method.cpp
        src/core/energy_calculators/qm_methods/gfnff_method.cpp
        src/core/energy_calculators/qm_methods/xtb_method.cpp
        src/core/energy_calculators/qm_methods/tblite_method.cpp
        src/core/energy_calculators/qm_methods/ulysses_method.cpp
        src/core/energy_calculators/qm_methods/dispersion_method.cpp
        src/core/energy_calculators/ff_methods/forcefield_method.cpp
        src/core/energycalculator.cpp
        src/core/solvation/gbsa.cpp
        src/core/molecule.cpp
        src/core/xyz_comment_parser.cpp
        src/core/fileiterator.cpp
        src/core/energy_calculators/ff_methods/forcefieldthread.cpp
        src/core/energy_calculators/ff_methods/forcefield.cpp
        src/core/energy_calculators/ff_methods/forcefieldfunctions.h
        src/core/energy_calculators/ff_methods/forcefieldgenerator.cpp
        src/core/energy_calculators/ff_methods/d3param_generator.cpp
        src/core/energy_calculators/ff_methods/d3_reference_cn.cpp
        src/core/energy_calculators/ff_methods/d3_reference_c6.cpp
        src/core/energy_calculators/ff_methods/cn_calculator.cpp
        src/core/energy_calculators/dispersion/d4param_generator.cpp
        src/core/energy_calculators/dispersion/d4_evaluator.cpp
        src/core/energy_calculators/dispersion/d4_charge_model.cpp
        src/core/energy_calculators/dispersion/d4_ncoord.cpp
        src/core/energy_calculators/ff_methods/cg_potentials.cpp
        src/core/energy_calculators/ff_methods/fast_exp.cpp
        src/core/curcuma_logger.cpp
        src/core/citation_database.cpp
        src/core/citation_registry.cpp
        src/tools/formats.h
        src/tools/geometry.h
        src/tools/general.h
        src/core/energy_calculators/qm_methods/orcainterface.h
        src/core/energy_calculators/qm_methods/orcainterface.cpp
        src/core/energy_calculators/qm_methods/orca_method.h
        src/core/energy_calculators/qm_methods/orca_method.cpp
        src/tools/trajectory_writer.cpp
        src/tools/bmt_utils.cpp
        src/core/imagewriter.cpp
        src/core/parameter_registry.cpp
        src/core/config_manager.cpp
        src/core/periodic_table.cpp
        src/core/functional_groups.cpp
)

set(curcuma_cap_SRC
        src/capabilities/optimisation/rf_solver.cpp
        src/capabilities/optimisation/lbfgs.cpp
        src/capabilities/persistentdiagram.cpp
        src/capabilities/analysenciplot.cpp
        src/capabilities/curcumamethod.cpp
        src/capabilities/confscan.cpp
        src/capabilities/confsearch.cpp
        src/capabilities/shared_bias_pool.cpp
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
        src/capabilities/md_diagnostics.cpp
        src/capabilities/analysis.cpp
        src/capabilities/analysis_output.cpp
        src/capabilities/handlers/general_properties_handler.cpp
        src/capabilities/handlers/scattering_handler.cpp
        src/capabilities/handlers/rdf_handler.cpp
        src/capabilities/handlers/topology_handler.cpp
        src/capabilities/trajectory_statistics.cpp
        src/capabilities/casino.cpp
        src/capabilities/scnp_parser.cpp
        src/capabilities/trajectoryanalysis.cpp
        src/capabilities/polymerbuild.cpp
        src/capabilities/tda_engine.cpp
        src/capabilities/hessian.cpp
        src/capabilities/qmdfffit.cpp
        src/core/imagewriter.cpp
)

add_library(curcuma_core  ${curcuma_core_SRC})

# Claude Generated (October 2025): Add feature flag compile definitions to curcuma_core target
# These are required for method_factory.cpp #ifdef checks to work properly
if(USE_GFNFF)
    target_compile_definitions(curcuma_core PRIVATE USE_GFNFF)
endif()

if(USE_XTB)
    target_compile_definitions(curcuma_core PRIVATE USE_XTB)
endif()

if(USE_TBLITE)
    target_compile_definitions(curcuma_core PRIVATE USE_TBLITE)
endif()

if(USE_ULYSSES)
    target_compile_definitions(curcuma_core PRIVATE USE_ULYSSES)
endif()

if(USE_D3)
    target_compile_definitions(curcuma_core PRIVATE USE_D3)
endif()

if(USE_D4)
    target_compile_definitions(curcuma_core PRIVATE USE_D4)
endif()

# Claude Generated (March 2026): BLAS/MKL compile definitions for native QM methods
# Enables EIGEN_USE_BLAS/EIGEN_USE_LAPACKE optimizations in curcuma_eigen_config.h
if(USE_BLAS)
    target_compile_definitions(curcuma_core PRIVATE USE_BLAS)
endif()

if(USE_MKL)
    target_compile_definitions(curcuma_core PRIVATE USE_MKL)
endif()

add_library(curcuma_cap  ${curcuma_cap_SRC})

# Enable precompiled headers for faster compilation
if(USE_PCH)
    target_precompile_headers(curcuma_core PRIVATE src/pch.h)
    target_precompile_headers(curcuma_cap PRIVATE src/pch.h)

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

# Claude Generated 2026 - km (Kuhn-Munkres/Hungarian assignment, c_code/hungarian.c)
# is dead code: nothing in the codebase includes its public interface.c::assign()
# outside c_code/ itself (RMSD reordering uses the header-only Eigen Munkres in
# rmsd/munkres.h instead). It also relies on C99 VLAs throughout (Boolean col_covered[m],
# a flexible array member, etc.), which MSVC's cl.exe has never supported. Skip it
# entirely on MSVC rather than porting a vendored library nothing calls.
if(NOT MSVC)
    set(km_SRC
        src/capabilities/c_code/hungarian.c
        src/capabilities/c_code/interface.c
    )
    add_library(km ${km_SRC})
    target_link_libraries(curcuma_core km)
endif()

if(USE_ULYSSES)
    add_library(ulysses_lib
        src/core/energy_calculators/qm_methods/interface/ulysses.cpp
    )
    target_compile_definitions(ulysses_lib PRIVATE USE_ULYSSES)
    if(USE_BLAS)
        target_link_libraries(ulysses_lib ${BLAS_LIBRARIES})
    endif()
    if(USE_MKL)
        target_include_directories(ulysses_lib PRIVATE ${MKLROOT}/include)
        target_link_libraries(ulysses_lib -m64 -Wl,--no-as-needed -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core)
    endif()

    add_executable(ulysses_helper
        src/helpers/ulysses_helper.cpp
    )
    target_link_libraries(ulysses_helper PUBLIC ulysses_lib curcuma_core)

    target_link_libraries(curcuma_core ulysses_lib)
endif()

add_executable(curcuma
        src/main.cpp
)

# Enable precompiled headers for main executable
if(USE_PCH)
    target_precompile_headers(curcuma REUSE_FROM curcuma_core)
    target_precompile_headers(curcuma REUSE_FROM curcuma_cap)

endif()

if(C17)
    set_property(TARGET curcuma_core PROPERTY CXX_STANDARD 17)
    set_property(TARGET curcuma_cap PROPERTY CXX_STANDARD 17)
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma_core PROPERTY CXX_STANDARD 14)
    set_property(TARGET curcuma_cap PROPERTY CXX_STANDARD 14)
    set_property(TARGET curcuma PROPERTY CXX_STANDARD 14)
endif(C17)

set(target_dir "bin")
set_property(TARGET curcuma PROPERTY LINK_FLAGS_RELEASE -s)
set_property(TARGET curcuma PROPERTY LINK_FLAGS_RELWITHDEBINFO -s)

if(USE_XTB)
    FetchContent_MakeAvailable(xtb)

    if(NOT TARGET mctc-lib)
      set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
      message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
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
        src/core/energy_calculators/qm_methods/xtbinterface.cpp
    )
    add_library(curcuma_xtb  ${curcuma_xtb_SRC})

    # Claude Generated (October 2025): Add USE_XTB compile definition
    target_compile_definitions(curcuma_xtb PRIVATE USE_XTB)

    target_link_libraries(curcuma_xtb PUBLIC "mctc-lib" lib-xtb-shared gfortran)
    target_link_libraries(curcuma_core curcuma_xtb)

endif()

if(USE_TBLITE)

    if(NOT TARGET mctc-lib)
      set("mctc-lib-url" "https://github.com/grimme-lab/mctc-lib")
      message(STATUS "Retrieving mctc-lib from ${mctc-lib-url}")
      FetchContent_Declare(
        "mctc-lib"
        GIT_REPOSITORY "${mctc-lib-url}"
        GIT_TAG "HEAD"
      )
      FetchContent_MakeAvailable("mctc-lib")
    endif()
    if(NOT TARGET tblite-lib)
    FetchContent_MakeAvailable(tblite)
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
        src/core/energy_calculators/qm_methods/tbliteinterface.cpp
    )
    add_library(curcuma_tblite  ${curcuma_tblite_SRC})

    # Claude Generated (October 2025): Add USE_TBLITE compile definition
    target_compile_definitions(curcuma_tblite PRIVATE USE_TBLITE)

    target_link_libraries(curcuma_tblite PUBLIC "mctc-lib" tblite-lib gfortran)
    target_link_libraries(curcuma_core curcuma_tblite)
endif()

if(USE_GFNFF)
    # External GFN-FF library integration
    if(NOT TARGET gfnff-lib)
        FetchContent_MakeAvailable(gfnff)
        include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/external/gfnff/include)
    endif()
    
    if(HELPERS)
        add_executable(gfnff_helper
               src/helpers/gfnff_helper.cpp
        )

        target_compile_definitions(gfnff_helper PRIVATE USE_GFNFF)
        if(C17)
            set_property(TARGET gfnff_helper PROPERTY CXX_STANDARD 17)
        else()
            set_property(TARGET gfnff_helper PROPERTY CXX_STANDARD 14)
        endif(C17)
        target_link_libraries(gfnff_helper PUBLIC gfnff-lib gfortran pthread)
    endif()
    
    set(curcuma_gfnff_SRC
        src/core/energy_calculators/qm_methods/gfnffinterface.cpp
        src/core/energy_calculators/qm_methods/external_gfnff_method.cpp
    )
    add_library(curcuma_gfnff ${curcuma_gfnff_SRC})

    # Claude Generated (October 2025): Add USE_GFNFF compile definition (CRITICAL FIX)
    target_compile_definitions(curcuma_gfnff PRIVATE USE_GFNFF)

    target_link_libraries(curcuma_gfnff PUBLIC gfnff-lib gfortran)
    target_link_libraries(curcuma_core curcuma_gfnff)

    # Claude Generated (May 2026, ICX-build OpenMP fix): when the C++ side is built with
    # Intel ICX, libiomp5 is the active OpenMP runtime, but gfnff-lib (built via the
    # FetchContent gfortran path) pulls in libgomp transitively. Two OpenMP runtimes in
    # one process corrupt each other's heap (we hit double-free in __kmp_invoke_microtask
    # on polymer MD). libiomp5 ships GOMP_* compat symbols (263 vs 192 in libgomp), so it
    # is a runtime drop-in replacement. We can't drop libgomp from the link line because
    # gfnff-lib carries it as a direct DT_NEEDED entry. Solution: a small wrapper script
    # that LD_PRELOADs libiomp5 before exec — this forces all GOMP_* calls to resolve
    # against libiomp5 and the libgomp instance becomes inert.
    set(_curcuma_need_iomp5_wrapper FALSE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        get_filename_component(_icx_bin_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
        get_filename_component(_icx_root_dir "${_icx_bin_dir}" DIRECTORY)
        find_library(_curcuma_iomp5
            NAMES iomp5
            HINTS "${_icx_root_dir}/lib" "${_icx_root_dir}/lib/intel64"
            NO_DEFAULT_PATH)
        if(_curcuma_iomp5)
            set(_curcuma_need_iomp5_wrapper TRUE)
        else()
            message(WARNING "ICX detected but libiomp5 not found near ${CMAKE_CXX_COMPILER}. "
                           "If you USE_GFNFF, polymer MD will crash with double-free in "
                           "__kmp_invoke_microtask.")
        endif()
    endif()
endif()

if(USE_D3)
    if(NOT TARGET s-dftd3)
        FetchContent_MakeAvailable(simple-dftd3)
    endif()
    if(HELPERS)
        add_executable(dftd3_helper
           src/helpers/dftd3_helper.cpp
    )
    target_link_libraries(dftd3_helper PUBLIC s-dftd3 gfortran pthread) 
    endif()
    set(curcuma_D3_SRC
           src/core/energy_calculators/qm_methods/dftd3interface.cpp
    )
    add_library(curcuma_D3  ${curcuma_D3_SRC})

    # Claude Generated (October 2025): Add USE_D3 compile definition
    target_compile_definitions(curcuma_D3 PRIVATE USE_D3)

    target_link_libraries(curcuma_D3 s-dftd3)
    target_link_libraries(curcuma_core curcuma_D3 gfortran pthread)
endif()

if(USE_D4)
    # Claude Generated (November 2025): Auto-discover LAPACK/BLAS for DFT-D4
    # First try to find system-installed libraries
    if(NOT DEFINED D4LIBS_DIR OR NOT DEFINED D4INCLUDE_DIR)
        find_package(LAPACK QUIET)
        find_package(BLAS QUIET)

        # Try to find LAPACKE C interface if available
        find_library(LAPACKE_LIBRARY lapacke
                     HINTS /usr/lib /usr/lib64 /usr/local/lib
                           ${LAPACK_ROOT}/lib ${BLAS_ROOT}/lib)

        if(LAPACK_FOUND AND BLAS_FOUND)
            set(LAPACKE_LIBRARIES "${LAPACKE_LIBRARY}" CACHE STRING "LAPACKE library")
            set(CBLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE STRING "CBLAS library")
            message(STATUS "DFT-D4: Found LAPACK/BLAS - ${LAPACK_LIBRARIES}")
            set(D4_AUTODISCOVERED TRUE)
        else()
            message(WARNING "DFT-D4: LAPACK/BLAS not found via find_package(). Will use manual D4LIBS_DIR/D4INCLUDE_DIR if set.")
        endif()
    else()
        message(STATUS "DFT-D4: Using manual library paths - D4LIBS_DIR=${D4LIBS_DIR}")
        include_directories(${D4INCLUDE_DIR})
        set(D4_AUTODISCOVERED FALSE)
    endif()

    FetchContent_MakeAvailable(cpp-d4)
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
        src/core/energy_calculators/qm_methods/dftd4interface.cpp
    )
    add_library(curcuma_d4  ${curcuma_d4_SRC})

    # Claude Generated (October 2025): Add USE_D4 compile definition
    target_compile_definitions(curcuma_d4 PRIVATE USE_D4)

    if(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 cblas lapacke )
    else()
	    target_link_libraries(curcuma_d4 PUBLIC libcpp_d4 ${LAPACKE_LIBRARIES}  ${CBLAS_LIBRARIES} )
    endif(DEFINED D4LIBS_DIR AND DEFINED D4INCLUDE_DIR)
    target_link_libraries(curcuma_core curcuma_d4)
endif()

if(USE_Plumed)
    # Build plumed2 from source via autotools (FetchContent only fetches, does not add_subdirectory)
    FetchContent_GetProperties(plumed2)
    if(NOT plumed2_POPULATED)
        FetchContent_Populate(plumed2)
    endif()

    set(PLUMED2_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/plumed2)
    set(PLUMED2_LIB        ${PLUMED2_SOURCE_DIR}/src/lib/libplumed.so)

    include(ExternalProject)
    ExternalProject_Add(plumed2_build
        SOURCE_DIR        ${PLUMED2_SOURCE_DIR}
        BINARY_DIR        ${PLUMED2_SOURCE_DIR}
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env "CXXFLAGS=-fpermissive -O3" ${PLUMED2_SOURCE_DIR}/configure --disable-python
        BUILD_COMMAND     make -C src/lib -j
        INSTALL_COMMAND   ""
        BUILD_BYPRODUCTS  ${PLUMED2_LIB}
    )

    add_custom_target(plumed2_patch_rpath
        COMMAND patchelf --set-rpath '$$ORIGIN' ${PLUMED2_LIB}
        DEPENDS plumed2_build
        COMMENT "Patching RPATH of libplumed.so"
    )

    add_library(plumed UNKNOWN IMPORTED GLOBAL)
    set_target_properties(plumed PROPERTIES
        IMPORTED_LOCATION             ${PLUMED2_LIB}
        INTERFACE_INCLUDE_DIRECTORIES ${PLUMED2_SOURCE_DIR}/src/wrapper
    )
    add_dependencies(plumed plumed2_build plumed2_patch_rpath)

    target_link_libraries(curcuma_core plumed)
    target_compile_definitions(curcuma_core PRIVATE USE_Plumed)

    message(STATUS "PLUMED2: external project configured, library expected at ${PLUMED2_LIB}")
endif()

# Claude Generated (March 2026): CUDA-accelerated GFN-FF (gfnff -gpu cuda)
# Requires: CUDA toolkit, compute capability >= 6.0 (Pascal / GTX 10xx)
# Build: cmake -DUSE_CUDA=ON && make -j4
if(USE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_STANDARD 17)
    # Target: Pascal(6.0), Volta(7.0), Turing(7.5), Ampere(8.0/8.6), Ada(8.9)
    set(CMAKE_CUDA_ARCHITECTURES "60;70;75;80;86;89;120" CACHE STRING "CUDA architectures")

    # CUDA source files - compiled as part of curcuma_core
    set(CUDA_SOURCES
        src/core/energy_calculators/ff_methods/cuda/gfnff_kernels.cu
        src/core/energy_calculators/ff_methods/cuda/ff_workspace_gpu.cu
        src/core/energy_calculators/ff_methods/cuda/eeq_solver_gpu.cu
        src/core/energy_calculators/ff_methods/cuda/gpu_utils.cpp
        src/core/energy_calculators/qm_methods/gfnff_gpu_method.cpp
    )
    # Claude Generated (2026-06): native xTB GPU kernels (compiled by nvcc).
    # The host-side wrapper (xtb_gpu_method.cpp) is added below as a normal source.
    if(USE_CUDA_XTB)
        list(APPEND CUDA_SOURCES
            src/core/energy_calculators/qm_methods/cuda/xtb_gpu_context.cu)
    endif()
    set_property(SOURCE ${CUDA_SOURCES} PROPERTY LANGUAGE CUDA)
    set_property(SOURCE ${CUDA_SOURCES} PROPERTY CUDA_SEPARABLE_COMPILATION ON)

    # Add CUDA sources to curcuma_core
    target_sources(curcuma_core PRIVATE ${CUDA_SOURCES})
    target_include_directories(curcuma_core PRIVATE ${PROJECT_SOURCE_DIR})
    target_link_libraries(curcuma_core CUDA::cudart CUDA::cusolver CUDA::cublas)
    target_compile_definitions(curcuma_core PUBLIC USE_CUDA)
    # curcuma_core is created before this block, so its CUDA_ARCHITECTURES target
    # property is still empty (CMP0104). Set it explicitly from the global list.
    # CUDA_RESOLVE_DEVICE_SYMBOLS: resolve device code within the static library so
    # the final curcuma binary does not need a device-link step and links as plain PIC.
    set_target_properties(curcuma_core PROPERTIES
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
        CUDA_RESOLVE_DEVICE_SYMBOLS ON)

    # Claude Generated (2026-06): native xTB GPU host wrapper + feature flag.
    if(USE_CUDA_XTB)
        target_sources(curcuma_core PRIVATE
            src/core/energy_calculators/qm_methods/xtb_gpu_method.cpp)
        target_compile_definitions(curcuma_core PUBLIC USE_CUDA_XTB)
        message(STATUS "CUDA native xTB enabled: gfn1/gfn2 -gpu cuda available")
    endif()

    message(STATUS "CUDA support enabled: gfnff -gpu cuda available")
    message(STATUS "CUDA architectures: ${CMAKE_CUDA_ARCHITECTURES}")
endif()

# Claude Generated (June 2026): ROCm/HIP backend (gfnff/gfn1/gfn2 -gpu rocm).
# Build: cmake -DUSE_ROCM=ON -DUSE_ROCM_XTB=ON -DCMAKE_PREFIX_PATH=/opt/rocm && make -j4
# The HIP context (device handshake + rocSOLVER eigensolver + the Stage-2 resident-SCF
# __global__ kernels) is compiled by hipcc into a plain relocatable object via a custom
# command, then linked into curcuma_core as an EXTERNAL_OBJECT with g++. This builds the
# device kernels while keeping the HIP-only --offload-arch / --hip-link flags off the GNU
# link of curcuma_core + the executables (which otherwise flips the link to ld.lld and
# drops the GNU OpenMP runtime libgomp). curcuma_core links libamdhip64 (+ rocSOLVER /
# rocBLAS when present); no enable_language(HIP) / CMake HIP machinery is used.
if(USE_ROCM)
    find_library(ROCM_HIP_RT amdhip64 REQUIRED PATHS /opt/rocm/lib /opt/rocm/lib64)
    find_package(rocblas QUIET)
    find_package(rocsolver QUIET)

    target_compile_definitions(curcuma_core PUBLIC USE_ROCM)
    if(USE_ROCM_XTB)
        find_program(HIPCC_EXECUTABLE hipcc REQUIRED PATHS /opt/rocm/bin)
        if(NOT DEFINED ROCM_GPU_ARCH)
            set(ROCM_GPU_ARCH "gfx1150" CACHE STRING "HIP offload arch (rocminfo | grep gfx)")
        endif()

        # hipcc compile flags + defines. HAVE_ROCSOLVER gates the eigensolver/resident code.
        set(_hipcc_defs -DUSE_ROCM_XTB -D__HIP_PLATFORM_AMD__)
        # Match the g++ build's Eigen alignment (CURCUMA_EIGEN_ALIGN, set with the march/AVX
        # switch above) so Eigen objects crossing this hipcc TU <-> g++ boundary share one ABI.
        if(NOT CURCUMA_EIGEN_ALIGN EQUAL 0)
            list(APPEND _hipcc_defs -DEIGEN_MAX_ALIGN_BYTES=${CURCUMA_EIGEN_ALIGN})
        endif()
        set(_rocm_wrap_defs "")
        if(rocsolver_FOUND AND rocblas_FOUND)
            find_library(ROCM_ROCSOLVER rocsolver REQUIRED PATHS /opt/rocm/lib /opt/rocm/lib64)
            find_library(ROCM_ROCBLAS   rocblas   REQUIRED PATHS /opt/rocm/lib /opt/rocm/lib64)
            list(APPEND _hipcc_defs -DHAVE_ROCSOLVER -DHAVE_ROCBLAS)
            list(APPEND _rocm_wrap_defs "HAVE_ROCSOLVER")
            message(STATUS "ROCm Stage 1/2 enabled (rocSOLVER dsygvd + HIP kernels)")
        else()
            message(STATUS "ROCm Stage 0 only: install rocsolver + rocblas for the GPU eigensolver")
        endif()

        # hipcc -c <ctx>.hip -> .o  (host + device code in one relocatable object).
        set(_hip_ctx_src ${PROJECT_SOURCE_DIR}/src/core/energy_calculators/qm_methods/rocm/xtb_hip_context.hip)
        set(_hip_ctx_hdr ${PROJECT_SOURCE_DIR}/src/core/energy_calculators/qm_methods/rocm/xtb_hip_context.h)
        set(_hip_ctx_obj ${CMAKE_CURRENT_BINARY_DIR}/xtb_hip_context.hip.o)
        add_custom_command(
            OUTPUT  ${_hip_ctx_obj}
            COMMAND ${HIPCC_EXECUTABLE} -c -fPIC -O2 -std=c++17 --offload-arch=${ROCM_GPU_ARCH}
                    ${_hipcc_defs} -I${CMAKE_SOURCE_DIR} -I${PROJECT_SOURCE_DIR} -I/opt/rocm/include
                    ${_hip_ctx_src} -o ${_hip_ctx_obj}
            DEPENDS ${_hip_ctx_src} ${_hip_ctx_hdr}
            COMMENT "Compiling HIP xTB context with hipcc (arch ${ROCM_GPU_ARCH})"
            VERBATIM)
        add_custom_target(curcuma_hip_ctx DEPENDS ${_hip_ctx_obj})
        set_source_files_properties(${_hip_ctx_obj} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)

        # Host wrapper (g++) + the prebuilt HIP object, both into curcuma_core.
        target_sources(curcuma_core PRIVATE
            src/core/energy_calculators/qm_methods/xtb_hip_method.cpp
            ${_hip_ctx_obj})
        add_dependencies(curcuma_core curcuma_hip_ctx)
        target_include_directories(curcuma_core PRIVATE /opt/rocm/include)
        target_link_libraries(curcuma_core ${ROCM_HIP_RT})
        target_compile_definitions(curcuma_core PUBLIC USE_ROCM_XTB)
        if(rocsolver_FOUND AND rocblas_FOUND)
            target_link_libraries(curcuma_core ${ROCM_ROCSOLVER} ${ROCM_ROCBLAS})
        endif()
        if(_rocm_wrap_defs)
            set_source_files_properties(
                src/core/energy_calculators/qm_methods/xtb_hip_method.cpp
                PROPERTIES COMPILE_DEFINITIONS "${_rocm_wrap_defs}")
        endif()

        message(STATUS "ROCm native xTB enabled: gfn1/gfn2 -gpu rocm available")
    endif()

    # ---- ROCm GFN-FF (gfnff -gpu rocm) -------------------------------------
    # Single-TU HIP port of cuda/{gfnff_kernels,ff_workspace_gpu,eeq_solver_gpu,gpu_utils}.
    # Same EXTERNAL_OBJECT pattern as the xTB block: hipcc -c -> one .o, g++ links it.
    # USE_D3/USE_D4 are NOT passed here — the .hip pulls them from this build's generated
    # global_config.h (-I${PROJECT_BINARY_DIR}), so the device TU and curcuma_core never
    # disagree on the dispersion backend. GFN-FF uses its own D4 (dispersion/), so configure
    # with -DUSE_D4=ON -DUSE_D3=OFF.
    if(USE_ROCM_GFNFF)
        if(NOT HIPCC_EXECUTABLE)
            find_program(HIPCC_EXECUTABLE hipcc REQUIRED PATHS /opt/rocm/bin)
        endif()
        if(NOT DEFINED ROCM_GPU_ARCH)
            set(ROCM_GPU_ARCH "gfx1150" CACHE STRING "HIP offload arch (rocminfo | grep gfx)")
        endif()

        # FMT_HEADER_ONLY: curcuma_core links fmt::fmt-header-only (INTERFACE define), but the
        # hipcc custom command does not inherit that target's compile defs, so set it here too —
        # otherwise the device TU references the non-inline fmt::vprint and the g++ link fails.
        set(_hip_ff_defs -DUSE_ROCM_GFNFF -D__HIP_PLATFORM_AMD__ -DFMT_HEADER_ONLY)
        # Match the g++ Eigen alignment (see the xTB block / the march-derived CURCUMA_EIGEN_ALIGN)
        # so workspace matrices read in the g++ wrapper (e.g. the gradient norm) share one ABI.
        if(NOT CURCUMA_EIGEN_ALIGN EQUAL 0)
            list(APPEND _hip_ff_defs -DEIGEN_MAX_ALIGN_BYTES=${CURCUMA_EIGEN_ALIGN})
        endif()
        if(rocsolver_FOUND AND rocblas_FOUND)
            if(NOT ROCM_ROCSOLVER)
                find_library(ROCM_ROCSOLVER rocsolver REQUIRED PATHS /opt/rocm/lib /opt/rocm/lib64)
                find_library(ROCM_ROCBLAS   rocblas   REQUIRED PATHS /opt/rocm/lib /opt/rocm/lib64)
            endif()
            list(APPEND _hip_ff_defs -DHAVE_ROCSOLVER -DHAVE_ROCBLAS)
        endif()

        set(_hip_ff_src ${PROJECT_SOURCE_DIR}/src/core/energy_calculators/ff_methods/rocm/gfnff_rocm.hip)
        set(_hip_ff_obj ${CMAKE_CURRENT_BINARY_DIR}/gfnff_rocm.hip.o)
        add_custom_command(
            OUTPUT  ${_hip_ff_obj}
            COMMAND ${HIPCC_EXECUTABLE} -c -fPIC -O2 -std=c++17 --offload-arch=${ROCM_GPU_ARCH}
                    ${_hip_ff_defs}
                    -I${CMAKE_SOURCE_DIR} -I${PROJECT_SOURCE_DIR} -I${PROJECT_BINARY_DIR}
                    -I${PROJECT_SOURCE_DIR}/external
                    -I${PROJECT_SOURCE_DIR}/external/fmt/include
                    -I${PROJECT_SOURCE_DIR}/external/eigen-3.4.0
                    -I${PROJECT_SOURCE_DIR}/external/cpp-d4/include
                    -I${PROJECT_SOURCE_DIR}/external/LBFGSpp/include
                    -I${PROJECT_SOURCE_DIR}/external/ripser
                    -I/opt/rocm/include
                    ${_hip_ff_src} -o ${_hip_ff_obj}
            DEPENDS ${_hip_ff_src}
            COMMENT "Compiling HIP GFN-FF kernels with hipcc (arch ${ROCM_GPU_ARCH})"
            VERBATIM)
        add_custom_target(curcuma_hip_gfnff DEPENDS ${_hip_ff_obj})
        set_source_files_properties(${_hip_ff_obj} PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)

        target_sources(curcuma_core PRIVATE
            src/core/energy_calculators/qm_methods/gfnff_hip_method.cpp
            ${_hip_ff_obj})
        add_dependencies(curcuma_core curcuma_hip_gfnff)
        target_include_directories(curcuma_core PRIVATE /opt/rocm/include)
        target_link_libraries(curcuma_core ${ROCM_HIP_RT})
        if(rocsolver_FOUND AND rocblas_FOUND)
            target_link_libraries(curcuma_core ${ROCM_ROCSOLVER} ${ROCM_ROCBLAS})
        endif()
        target_compile_definitions(curcuma_core PUBLIC USE_ROCM_GFNFF)
        # The g++-compiled host wrapper needs the HIP platform macro for <hip/hip_runtime.h>
        # (single host-side hipMemcpy for the GPU-Schur charge download).
        set(_ff_wrap_defs "__HIP_PLATFORM_AMD__")
        if(rocsolver_FOUND AND rocblas_FOUND)
            list(APPEND _ff_wrap_defs "HAVE_ROCSOLVER")
        endif()
        set_source_files_properties(
            src/core/energy_calculators/qm_methods/gfnff_hip_method.cpp
            PROPERTIES COMPILE_DEFINITIONS "${_ff_wrap_defs}")
        message(STATUS "ROCm GFN-FF enabled: gfnff -gpu rocm available")
    endif()

    message(STATUS "ROCm: rocBLAS=${rocblas_FOUND} rocSOLVER=${rocsolver_FOUND}")
endif()

# Claude Generated (June 2026): Vulkan compute backend (gfnff/gfn1/gfn2 -gpu vulkan).
# Build: cmake -DUSE_VULKAN=ON -DUSE_VULKAN_XTB=ON && make -j4
# Stage 0 wires the host Vulkan context (instance/device/queue handshake) + native-xTB
# pass-through wrapper. The hand-written SPIR-V compute shaders under vulkan/shaders/
# are compiled+embedded starting at Stage 1 (the eigensolver path); not referenced here.
if(USE_VULKAN)
    find_package(Vulkan REQUIRED)

    set(VULKAN_SOURCES
        src/core/energy_calculators/qm_methods/vulkan/vk_context.cpp)
    if(USE_VULKAN_XTB)
        list(APPEND VULKAN_SOURCES
            src/core/energy_calculators/qm_methods/vulkan/xtb_vulkan_context.cpp
            src/core/energy_calculators/qm_methods/xtb_vulkan_method.cpp)
    endif()

    target_sources(curcuma_core PRIVATE ${VULKAN_SOURCES})
    target_include_directories(curcuma_core PRIVATE ${PROJECT_SOURCE_DIR} ${Vulkan_INCLUDE_DIRS})
    target_link_libraries(curcuma_core Vulkan::Vulkan)
    target_compile_definitions(curcuma_core PUBLIC USE_VULKAN)

    if(USE_VULKAN_XTB)
        target_compile_definitions(curcuma_core PUBLIC USE_VULKAN_XTB)
        message(STATUS "Vulkan native xTB enabled: gfn1/gfn2 -gpu vulkan available")
    endif()
    message(STATUS "Vulkan support enabled: -gpu vulkan available (Stage 0)")
endif()

target_link_libraries(curcuma_core pthread fmt::fmt-header-only )
target_link_libraries(curcuma_cap pthread fmt::fmt-header-only )

# Link BLAS/LAPACK for curcuma_core.
# EIGEN_USE_BLAS: per-file in eeq_solver.cpp (avoids Eigen 3.4 complex-type conflict).
# EIGEN_USE_LAPACKE: per-target on curcuma_core — safe because no TU uses complex Eigen types
#   alongside LAPACKE. Enables LAPACK dsyevd (divide-and-conquer) inside SelfAdjointEigenSolver
#   for ~3-5x speedup on matrices with N > 300 (Tier M/L ANC regeneration fallback).
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
# WP1 (2026-06): when USE_MKL, link ONLY the explicit threaded MKL below
# (mkl_gnu_thread). find_package(BLAS) resolves to the SEQUENTIAL MKL
# (libmkl_sequential), which — linked alongside the threaded layer — loads first and
# wins, forcing every BLAS/LAPACK call (incl. the SCF eigensolve) single-threaded.
if(BLAS_FOUND AND NOT USE_MKL)
    target_link_libraries(curcuma_core ${BLAS_LIBRARIES})
    target_compile_definitions(curcuma_core PUBLIC HAVE_BLAS)
    message(STATUS "BLAS linked for curcuma_core")
endif()
if(LAPACK_FOUND AND NOT USE_MKL)
    target_link_libraries(curcuma_core ${LAPACK_LIBRARIES})
    message(STATUS "LAPACK linked for curcuma_core")
endif()

# LAPACKE C interface — needed for EIGEN_USE_LAPACKE. Re-use the LAPACKE_LIBRARY already
# discovered by the DFT-D4 search above, or do a fresh search if that block was skipped.
if(NOT LAPACKE_LIBRARY)
    find_library(LAPACKE_LIBRARY lapacke
                 HINTS /usr/lib /usr/lib64 /usr/local/lib)
endif()
if(LAPACKE_LIBRARY)
    target_link_libraries(curcuma_core ${LAPACKE_LIBRARY})
    target_compile_definitions(curcuma_core PUBLIC HAVE_LAPACKE)
    # EIGEN_USE_LAPACKE is NOT set globally here — it breaks TUs that use local variables
    # named 'I' (C99 complex macro conflict with LAPACKE headers). Define it per-file instead:
    #   rf_solver.cpp, lbfgs.cpp — both safe (no variable named 'I').
    # See ancopt_optimizer.cpp: renamed Matrix3d I -> inertia to fix the same issue.
    # Source files guard their EIGEN_USE_LAPACKE with #ifdef HAVE_LAPACKE so the build
    # works without LAPACKE installed.
    message(STATUS "LAPACKE found: ${LAPACKE_LIBRARY} — linking, EIGEN_USE_LAPACKE applied per-file")
else()
    message(STATUS "LAPACKE not found — SelfAdjointEigenSolver uses Eigen built-in path")
endif()


# NOT USE_MKL: the threaded MKL block below already provides BLAS/LAPACK symbols, and
# linking the find_package (sequential MKL) result here would override its threading (WP1).
if(USE_BLAS AND NOT USE_MKL)
    if(BLAS_FOUND)
        target_link_libraries(curcuma_core ${BLAS_LIBRARIES})
        if(LAPACK_FOUND)
            target_link_libraries(curcuma_core ${LAPACK_LIBRARIES})
        endif()
    else()
        # Fallback to hardcoded library name if find_package failed but library exists
        target_link_libraries(curcuma_core blas)
    endif()
endif()
if(USE_MKL)
    target_include_directories(curcuma_core PRIVATE ${MKLROOT}/include)
    target_include_directories(curcuma_cap PRIVATE ${MKLROOT}/include)
    target_link_libraries(curcuma_core -m64 -Wl,--no-as-needed -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core)
endif()
if(WIN32) # Check if we are on Windows
else()
     target_link_libraries(curcuma_core dl )
endif(WIN32)

target_link_libraries(curcuma fmt::fmt-header-only)

target_link_libraries(curcuma curcuma_core)

# Claude Generated (May 2026, ICX + USE_GFNFF OpenMP fix): when both Intel ICX (libiomp5)
# and gfortran-built gfnff-lib (libgomp DT_NEEDED) end up in the same binary, the two
# OpenMP runtimes corrupt each other's heap. We hit double-free in __kmp_invoke_microtask
# on polymer MD. libiomp5 ships GOMP_* compat symbols (263 vs 192 in libgomp), so it is a
# runtime drop-in replacement. Use patchelf to rewrite curcuma's DT_NEEDED so it loads
# libiomp5 instead of libgomp.so.1 — single binary, no wrapper, no LD_PRELOAD needed.
if(_curcuma_need_iomp5_wrapper)
    find_program(PATCHELF_EXECUTABLE patchelf)
    if(PATCHELF_EXECUTABLE)
        add_custom_command(TARGET curcuma POST_BUILD
            COMMAND "${PATCHELF_EXECUTABLE}" --replace-needed libgomp.so.1 libiomp5.so
                    "$<TARGET_FILE:curcuma>"
            COMMENT "ICX + GFN-FF: patchelf replace-needed libgomp.so.1 → libiomp5.so "
                    "to avoid double-OMP-runtime heap corruption"
            VERBATIM)
        message(STATUS "ICX + GFN-FF: patchelf will rewrite libgomp→libiomp5 in curcuma "
                       "post-build (avoids double-runtime heap corruption).")
    else()
        message(WARNING "ICX + USE_GFNFF detected but patchelf not found. The curcuma "
                        "binary will crash with double-free in MD parallel regions. "
                        "Install patchelf, OR run with: "
                        "LD_PRELOAD=${_curcuma_iomp5} ./curcuma ...")
    endif()
endif()
add_executable(curcuma_helper
        src/helpers/main.cpp)
        
# cli_helper disabled - requires cli_parser.h (native-gfnff WIP)
# add_executable(cli_helper
#         src/helpers/cli_test.cpp)

add_executable(parallel_scf
        src/helpers/parallel_scf.cpp)
target_link_libraries(parallel_scf curcuma_core)

target_link_libraries(curcuma_helper fmt::fmt-header-only)
target_link_libraries(curcuma_helper curcuma_core)

# Enable precompiled headers for helper executables
if(USE_PCH)
    target_precompile_headers(curcuma_helper REUSE_FROM curcuma_core)
    # target_precompile_headers(cli_helper REUSE_FROM curcuma_core)
    target_precompile_headers(parallel_scf REUSE_FROM curcuma_core)
endif()

set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELEASE -s)
set_property(TARGET curcuma_helper PROPERTY LINK_FLAGS_RELWITHDEBINFO -s)

if(C17)
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 17)
else()
    set_property(TARGET curcuma_helper PROPERTY CXX_STANDARD 14)
endif(C17)

include(CTest)
add_subdirectory(test_cases)

target_sources(Continuous
  PRIVATE
    src/core/energy_calculators/ff_methods/qmdff_par.h
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
# Claude Generated (June 2026): The molalign scenario shells out to the external
# 'molalign' binary. Only register it when that binary is found, otherwise the
# test fails on machines without the optional dependency (false negative).
find_program(MOLALIGN_BIN molalign)
if(MOLALIGN_BIN)
    add_test(NAME confscan_molalign COMMAND confscan_test molalign WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
else()
    message(STATUS "molalign binary not found - skipping confscan_molalign test")
endif()
add_test(NAME confscan_sLX1 COMMAND confscan_test sLX1 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX2 COMMAND confscan_test sLX2 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX2Reset COMMAND confscan_test sLX2Reset WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX20 COMMAND confscan_test sLX20 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)
add_test(NAME confscan_sLX20Reset COMMAND confscan_test sLX20Reset WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

add_test(NAME confstat COMMAND confstat_tests WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

# Claude Generated: Comprehensive test for Molecule class before refactoring
add_test(NAME energy_methods COMMAND test_energy_methods WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test_cases)

# Claude Generated: Parameter I/O system comprehensive tests (October 2025)
add_test(
    NAME parameter_io_tests
    COMMAND bash ${CMAKE_SOURCE_DIR}/test_cases/test_parameter_io.sh
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)

# --- Auto-generated Parameter Registry (Plan 4.1) ---

# 1. Define the parser executable
add_executable(curcuma_param_parser scripts/param_parser/main.cpp)

# 2. Get all header files in the project
file(GLOB_RECURSE ALL_PROJECT_HEADERS "src/*.h" "src/**/*.h")

# Claude Generated (October 2025): Ensure generated directory exists before writing output file
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")

# 3. Define the custom command to run the parser
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h
    COMMAND curcuma_param_parser
            --output ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h
            --inputs ${ALL_PROJECT_HEADERS}
    DEPENDS curcuma_param_parser ${ALL_PROJECT_HEADERS}
    COMMENT "Generating C++ parameter registry from source headers"
)

# 4. Create a target for the generation step
add_custom_target(GenerateParams DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/generated/parameter_registry.h)

# 5. Make the main library depend on the generated file
add_dependencies(curcuma_core GenerateParams)

# 6. Add the generated directory to the include paths
include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)
target_include_directories(curcuma_core SYSTEM BEFORE PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/generated)

# --- End of Auto-generated Parameter Registry ---

# Claude Generated 2026 - only install curcuma's own standalone executables when
# curcuma is the top-level project (its own release/ build), not when it's
# consumed as a FetchContent dependency (e.g. by qurcuma, which excludes both
# targets from the default build via EXCLUDE_FROM_ALL). Otherwise the consuming
# project's own `install`/`cpack` aggregates this rule and fails trying to
# install a curcuma_helper binary that was never built.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    install(TARGETS curcuma RUNTIME DESTINATION bin)
    install(TARGETS curcuma_helper RUNTIME DESTINATION bin)
endif()
