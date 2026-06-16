#!/usr/bin/env bash
# Regenerate the embedded SPIR-V (*.spv.inc) from the GLSL compute shaders.
# Run after editing any *.comp. Requires glslc (Vulkan SDK / shaderc).
# The committed *.spv.inc are #included by spirv_kernels.h, so the curcuma build
# itself does NOT need glslc.
set -euo pipefail
cd "$(dirname "$0")"
for s in angles col row vec gemm scale_cols fock popband \
         cn self_energy overlap_h0 gamma \
         grad_rep grad_coulomb grad_pulay multipole_ints \
         fock_multipole multipole_moments \
         angles_f32 col_f32 row_f32 vec_f32 \
         tri_matvec tri_rank2 tri_applyl tri_house tri_kw; do
    glslc --target-env=vulkan1.1 -mfmt=c -fshader-stage=compute "$s.comp" -o "$s.spv.inc"
    echo "compiled $s.comp -> $s.spv.inc"
done
echo "done. Rebuild curcuma to pick up the new SPIR-V."
