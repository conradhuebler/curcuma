#!/usr/bin/env bash
# Claude Generated (Sep 2026)
#
# Find cuSOLVERMp / cuBLASMp / NCCL on this machine and print the cmake flags for them.
#
# Background: curcuma can split ONE large molecule's eigensolve over several GPUs, which needs
# those three libraries. They are OPTIONAL - without them single-GPU runs, the batch spread over
# GPUs (-gpu_devices) and the distributed density (-gpu_density_devices) all still work - but the
# fallback for that one case is cusolverMg, measured 15x slower than a single GPU on polymer_2x.
#
# What is where, and why this script exists:
#   cusolverMg    IS part of the CUDA toolkit         -> a plain `module load cuda` has it
#   cuSOLVERMp    is NOT part of the CUDA toolkit     -> HPC SDK, standalone tarball, conda, pip
#   cuBLASMp      is NOT part of the CUDA toolkit     -> same
#   NCCL          is NOT part of the CUDA toolkit, but many clusters put it in the CUDA module
#                                                        or offer it as its own module
#
# Usage:  bash scripts/find_mgpu_libs.sh            (searches the usual places)
#         bash scripts/find_mgpu_libs.sh /extra/prefix ...
set -u

echo "=== curcuma: looking for cuSOLVERMp / cuBLASMp / NCCL ==="
echo

# ---------------------------------------------------------------- search roots
roots=()
for v in CUDA_HOME CUDA_ROOT CUDA_PATH CUDATOOLKIT_HOME NVHPC_ROOT NVHPC_HOME HPCSDK_HOME \
         NCCL_HOME NCCL_ROOT CUSOLVERMP_ROOT CUBLASMP_ROOT CONDA_PREFIX VIRTUAL_ENV; do
    [ -n "${!v:-}" ] && roots+=("${!v}")
done
for d in /usr /usr/local/cuda /opt/cuda /opt/nvidia/hpc_sdk "$HOME/.local"; do
    [ -d "$d" ] && roots+=("$d")
done
# python wheels, wherever the interpreter puts them
if command -v python3 >/dev/null 2>&1; then
    sp=$(python3 -c "import site,sys; print(next(iter(site.getsitepackages()),''))" 2>/dev/null)
    [ -n "$sp" ] && [ -d "$sp/nvidia" ] && roots+=("$sp/nvidia")
fi
for extra in "$@"; do roots+=("$extra"); done

if [ ${#roots[@]} -eq 0 ]; then
    echo "No search roots. Load your modules first (module load cuda / nvhpc), or pass a prefix."
    exit 1
fi

echo "Search roots:"
printf '  %s\n' "${roots[@]}" | sort -u
echo

# ------------------------------------------------------------------ the search
# find_one <pretty name> <library glob> <header name> -> sets FOUND_LIB / FOUND_INC
find_one() {
    local name="$1" libglob="$2" header="$3"
    FOUND_LIB=""; FOUND_INC=""
    local r
    for r in $(printf '%s\n' "${roots[@]}" | sort -u); do
        [ -z "$FOUND_LIB" ] && FOUND_LIB=$(find "$r" -maxdepth 6 -name "$libglob" 2>/dev/null | head -1)
        [ -z "$FOUND_INC" ] && FOUND_INC=$(find "$r" -maxdepth 6 -name "$header" 2>/dev/null | head -1)
        [ -n "$FOUND_LIB" ] && [ -n "$FOUND_INC" ] && break
    done
    # the loader's cache catches system packages that live outside the roots
    if [ -z "$FOUND_LIB" ] && command -v ldconfig >/dev/null 2>&1; then
        FOUND_LIB=$(ldconfig -p 2>/dev/null | awk -v g="${libglob%%.*}" '$1 ~ g {print $NF; exit}')
    fi
    if [ -n "$FOUND_LIB" ] || [ -n "$FOUND_INC" ]; then
        echo "$name:"
        echo "    library: ${FOUND_LIB:-NOT FOUND}"
        echo "    header : ${FOUND_INC:-NOT FOUND (a runtime-only package has no header - curcuma needs it)}"
    else
        echo "$name: not found"
    fi
}

find_one "cuSOLVERMp" "libcusolverMp.so*" "cusolverMp.h";  MP_LIB="$FOUND_LIB"; MP_INC="$FOUND_INC"
find_one "cuBLASMp"   "libcublasmp.so*"   "cublasmp.h";    BMP_LIB="$FOUND_LIB"; BMP_INC="$FOUND_INC"
find_one "NCCL"       "libnccl.so*"       "nccl.h";        NCCL_LIB="$FOUND_LIB"; NCCL_INC="$FOUND_INC"
find_one "cusolverMg (CUDA toolkit, fallback)" "libcusolverMg.so*" "cusolverMg.h"
echo

# ------------------------------------------------------------------- the verdict
root_of() {   # <found file> <subdir to strip: lib or include>
    [ -z "$1" ] && return
    local d; d=$(dirname "$1")
    case "$d" in */lib|*/lib64|*/include) dirname "$d" ;; *) echo "$d" ;; esac
}
if [ -n "$MP_LIB" ] && [ -n "$MP_INC" ] && [ -n "$BMP_LIB" ] && [ -n "$BMP_INC" ] \
   && [ -n "$NCCL_LIB" ] && [ -n "$NCCL_INC" ]; then
    echo "All three found. Configure with:"
    echo
    echo "  cmake .. -DUSE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=<your cc, e.g. 90 for H100/H200> \\"
    echo "        -DCUSOLVERMP_ROOT=$(root_of "$MP_LIB") \\"
    echo "        -DCUBLASMP_ROOT=$(root_of "$BMP_LIB") \\"
    echo "        -DNCCL_ROOT=$(root_of "$NCCL_LIB") \\"
    echo "        -DCURCUMA_REQUIRE_MULTI_GPU_EIGENSOLVER=ON"
    echo
    echo "Then check the summary:  cmake .. 2>&1 | grep '=== curcuma multi-GPU eigensolver'"
else
    echo "Not complete. This is not an error - curcuma builds and runs without them, and only"
    echo "the split of ONE molecule's eigensolve is unavailable (use -gpu_eigensolver_devices"
    echo "none there; single GPU, -gpu_devices batches and -gpu_density_devices are unaffected)."
    echo
    echo "To get them, in the order that is least work on a cluster:"
    echo "  1. module avail 2>&1 | grep -iE 'nvhpc|hpc.sdk|nccl|cuda'      # an nvhpc module has ALL THREE"
    echo "     module load nvhpc; bash $0                                  # then run this again"
    echo "  2. conda install -c nvidia libcusolvermp-dev libcublasmp-dev nccl"
    echo "  3. pip install nvidia-cusolvermp-cu13 nvidia-cublasmp-cu13 nvidia-nccl-cu13"
    echo "     (note: 'nvidia-cusolver' is the SINGLE-GPU library and not what is needed here,"
    echo "      and the suffix-less 'nvidia-nccl' on PyPI is a placeholder warning package)"
    echo "  4. tarballs from developer.nvidia.com/cusolvermp-downloads"
fi
