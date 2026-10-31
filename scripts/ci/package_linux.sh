#!/bin/bash
# package_linux.sh <build dir> <out dir> <isa>: pack the curcuma executable and each GPU plugin
# into its own archive, all with the top directory curcuma-<isa>/ so that unpacking the plugin
# archives next to the base archive merges them into one flat directory.
#   curcuma-linux-<isa>.tar.gz          executable and the non-glibc libraries it needs
#   curcuma-plugin-cuda-<isa>.tar.gz    libcurcuma_cuda.so and the CUDA libraries (if built)
#   curcuma-plugin-rocm-<isa>.tar.gz    libcurcuma_rocm.so, ROCm libraries, rocBLAS kernels (if built)
# Flat layout because the core dlopen()s libcurcuma_<backend>.so from the executable's directory
# (gpu_plugin.cpp) and rocBLAS looks for rocblas/library next to its own .so. A plugin only fits
# the executable of the same <isa> (shared Eigen alignment). Vulkan is not packaged.
# Claude Generated.
set -euo pipefail
BUILD="$1"; OUT="$2"; ISA="$3"
STAGE="$OUT/stage"
TOP="curcuma-$ISA"
rm -rf "$STAGE"; mkdir -p "$STAGE/base/$TOP" "$OUT"

# Libraries that must come from the target system (glibc and the loader).
SYSTEM='^(libc|libm|libdl|librt|libpthread|libutil|libresolv|libnsl|libcrypt|ld-linux-x86-64|linux-vdso)\.so'

# bundle <dir>: copy every non-glibc library the files in <dir> need, unless the base
# directory already has it.
bundle() {
    local dir="$1"
    for round in 1 2 3; do
        for f in "$dir"/*.so* "$dir/curcuma"; do
            [ -e "$f" ] || continue
            ldd "$f" | awk '/=> \//{print $3}' | while read -r lib; do
                base=$(basename "$lib")
                [[ "$base" =~ $SYSTEM ]] && continue
                [ -e "$STAGE/base/$TOP/$base" ] || [ -e "$dir/$base" ] || cp -L "$lib" "$dir/$base"
            done
        done
    done
}

cp "$BUILD/curcuma" "$STAGE/base/$TOP/"
bundle "$STAGE/base/$TOP"

for b in cuda rocm; do
    [ -f "$BUILD/libcurcuma_$b.so" ] || continue
    d="$STAGE/$b/$TOP"; mkdir -p "$d"
    cp "$BUILD/libcurcuma_$b.so" "$d/"
    bundle "$d"
done

# rocBLAS kernels (RDNA only; the gfx9 kernels are the bulk of the library and unused).
if [ -d "$STAGE/rocm/$TOP" ] && [ -d /opt/rocm/lib/rocblas/library ]; then
    mkdir -p "$STAGE/rocm/$TOP/rocblas"
    cp -r /opt/rocm/lib/rocblas/library "$STAGE/rocm/$TOP/rocblas/"
    rm -f "$STAGE"/rocm/"$TOP"/rocblas/library/*gfx9* "$STAGE"/rocm/"$TOP"/rocblas/library/*gfx8* "$STAGE"/rocm/"$TOP"/rocblas/library/*gfx7*
fi

for f in "$STAGE"/*/"$TOP"/*.so* "$STAGE/base/$TOP/curcuma"; do
    patchelf --set-rpath '$ORIGIN' "$f" 2>/dev/null || true
done

tar -C "$STAGE/base" -czf "$OUT/curcuma-linux-$ISA.tar.gz" "$TOP"
for b in cuda rocm; do
    [ -d "$STAGE/$b" ] && tar -C "$STAGE/$b" -czf "$OUT/curcuma-plugin-$b-$ISA.tar.gz" "$TOP"
done
for d in "$STAGE"/*; do du -sh "$d"; done
du -ah "$STAGE" | sort -rh | head -n 12 || true
ls -l "$OUT"/*.tar.gz
