#!/bin/bash
# package_linux.sh <build dir> <target dir>: copy curcuma, the GPU plugins and every
# non-glibc library they need into one flat directory and point their RUNPATH at it.
# Flat layout because the core dlopen()s libcurcuma_<backend>.so from the executable's
# directory (gpu_plugin.cpp) and rocBLAS looks for rocblas/library next to its own .so.
# Claude Generated.
set -euo pipefail
BUILD="$1"; DEST="$2"
rm -rf "$DEST"; mkdir -p "$DEST"
cp "$BUILD/curcuma" "$DEST/"
for p in libcurcuma_cuda.so libcurcuma_rocm.so libcurcuma_vulkan.so; do
    [ -f "$BUILD/$p" ] && cp "$BUILD/$p" "$DEST/"
done

# Libraries that must come from the target system (glibc and the loader).
SYSTEM='^(libc|libm|libdl|librt|libpthread|libutil|libresolv|libnsl|libcrypt|ld-linux-x86-64|linux-vdso)\.so'

for round in 1 2 3; do
    for f in "$DEST"/*.so* "$DEST/curcuma"; do
        [ -e "$f" ] || continue
        ldd "$f" | awk '/=> \//{print $3}' | while read -r lib; do
            base=$(basename "$lib")
            [[ "$base" =~ $SYSTEM ]] && continue
            [ -e "$DEST/$base" ] || cp -L "$lib" "$DEST/$base"
        done
    done
done

# rocBLAS kernels for the cards the plugin was built for (RDNA only; the gfx9 kernels
# are the bulk of the library and unused).
if [ -d /opt/rocm/lib/rocblas/library ]; then
    mkdir -p "$DEST/rocblas"
    cp -r /opt/rocm/lib/rocblas/library "$DEST/rocblas/"
    rm -f "$DEST"/rocblas/library/*gfx9* "$DEST"/rocblas/library/*gfx8* "$DEST"/rocblas/library/*gfx7*
fi

for f in "$DEST"/*.so* "$DEST/curcuma"; do
    patchelf --set-rpath '$ORIGIN' "$f" 2>/dev/null || true
done
du -sh "$DEST"
du -ah "$DEST" | sort -rh | head -n 15 || true
