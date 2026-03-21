/*
 * <GFN-FF GPU Struct-of-Arrays (SoA) Data Structures>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): SoA layout for GPU-coalesced memory access.
 * Converts GFN-FF AoS parameter structs to GPU-friendly SoA with RAII cudaMalloc.
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#pragma once

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>
#include <vector>

// Forward-declare GFN-FF parameter structs (defined in forcefieldthread.h / gfnff_parameters.h)
struct Bond;
struct Angle;
struct Dihedral;
struct Inversion;
struct GFNFFDispersion;
struct GFNFFRepulsion;
struct GFNFFCoulomb;

// ============================================================================
// CudaBuffer<T>: RAII wrapper around cudaMalloc / cudaFree
//
// Claude Generated (March 2026): GPU memory owner, non-copyable, moveable.
// ============================================================================
template <typename T>
struct CudaBuffer {
    T*  ptr = nullptr;
    int n   = 0;

    CudaBuffer()                               = default;
    CudaBuffer(const CudaBuffer&)              = delete;
    CudaBuffer& operator=(const CudaBuffer&)   = delete;

    CudaBuffer(CudaBuffer&& o) noexcept : ptr(o.ptr), n(o.n) {
        o.ptr = nullptr; o.n = 0;
    }
    CudaBuffer& operator=(CudaBuffer&& o) noexcept {
        if (this != &o) { free(); ptr = o.ptr; n = o.n; o.ptr = nullptr; o.n = 0; }
        return *this;
    }
    ~CudaBuffer() { free(); }

    /// Allocate count elements on GPU (frees previous allocation)
    void alloc(int count) {
        free();
        if (count <= 0) return;
        cudaError_t err = cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != cudaSuccess)
            throw std::runtime_error(std::string("cudaMalloc failed: ") + cudaGetErrorString(err));
        n = count;
    }

    /// Upload count elements from host to GPU (auto-allocates if needed)
    void upload(const T* host_data, int count, cudaStream_t stream = nullptr) {
        if (count <= 0) return;
        if (n < count) alloc(count);
        cudaError_t err = cudaMemcpyAsync(ptr, host_data, count * sizeof(T),
                                          cudaMemcpyHostToDevice, stream);
        if (err != cudaSuccess)
            throw std::runtime_error(std::string("cudaMemcpyAsync H2D failed: ") + cudaGetErrorString(err));
    }

    /// Upload from std::vector
    void upload(const std::vector<T>& v, cudaStream_t stream = nullptr) {
        upload(v.data(), static_cast<int>(v.size()), stream);
    }

    /// Download count elements from GPU to host
    void download(T* host_data, int count, cudaStream_t stream = nullptr) const {
        if (count <= 0 || !ptr) return;
        cudaError_t err = cudaMemcpyAsync(host_data, ptr, count * sizeof(T),
                                          cudaMemcpyDeviceToHost, stream);
        if (err != cudaSuccess)
            throw std::runtime_error(std::string("cudaMemcpyAsync D2H failed: ") + cudaGetErrorString(err));
    }

    /// Zero-fill count elements on GPU
    void zero(int count, cudaStream_t stream = nullptr) {
        if (count <= 0 || !ptr) return;
        cudaMemsetAsync(ptr, 0, count * sizeof(T), stream);
    }

    void free() {
        if (ptr) { cudaFree(ptr); ptr = nullptr; n = 0; }
    }

    operator T*()             { return ptr; }
    operator const T*() const { return ptr; }
    bool empty() const        { return n == 0 || ptr == nullptr; }
};

// ============================================================================
// DispersionSoA: GFNFFDispersion AoS → SoA for coalesced GPU access
// ============================================================================
struct DispersionSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> C6, r4r2ij, r0_sq, zetac6, r_cut;
    int n = 0;

    void upload(const std::vector<GFNFFDispersion>& v, cudaStream_t stream = nullptr);
};

// ============================================================================
// RepulsionSoA: GFNFFRepulsion AoS → SoA (bonded and nonbonded share struct)
// ============================================================================
struct RepulsionSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> alpha, repab, r_cut;
    int n = 0;

    void upload(const std::vector<GFNFFRepulsion>& v, cudaStream_t stream = nullptr);
};

// ============================================================================
// CoulombSoA: GFNFFCoulomb AoS → SoA (charges come from per-step d_charges)
// ============================================================================
struct CoulombSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> gamma_ij, r_cut;
    int n = 0;

    void upload(const std::vector<GFNFFCoulomb>& v, cudaStream_t stream = nullptr);
};

// ============================================================================
// BondSoA: Bond AoS → SoA (includes CN-dependent r0 parameters for GFN-FF)
// ============================================================================
struct BondSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> r0;         ///< Static r0_ij (used when CN unavailable)
    CudaBuffer<double> fc;         ///< Force constant
    CudaBuffer<double> alpha;      ///< Gaussian exponent
    CudaBuffer<double> rabshift;   ///< GFN-FF topology radius shift
    CudaBuffer<double> ff;         ///< EN-correction factor
    CudaBuffer<double> cnfak_i;    ///< CN sensitivity of atom i's covalent radius
    CudaBuffer<double> cnfak_j;    ///< CN sensitivity of atom j's covalent radius
    CudaBuffer<double> r0_base_i;  ///< Base r0 for atom i (r0_gfnff[Z-1])
    CudaBuffer<double> r0_base_j;  ///< Base r0 for atom j
    int n = 0;

    void upload(const std::vector<Bond>& v, cudaStream_t stream = nullptr);
};

// ============================================================================
// AngleSoA: Angle AoS → SoA (atom types needed for rcov distance damping)
// ============================================================================
struct AngleSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k;
    CudaBuffer<int>    ati, atj, atk;  ///< Atomic numbers for rcov lookup
    CudaBuffer<double> fc, theta0;
    int n = 0;

    void upload(const std::vector<Angle>& v,
                const std::vector<int>& atom_types,
                cudaStream_t stream = nullptr);
};

// ============================================================================
// DihedralSoA: Dihedral AoS → SoA (standard + extra combined, atom types for damping)
// ============================================================================
struct DihedralSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k, idx_l;
    CudaBuffer<int>    ati, atj, atk, atl;  ///< Atomic numbers for rcov lookup
    CudaBuffer<double> V;         ///< Barrier height
    CudaBuffer<double> phi0;      ///< Reference angle
    CudaBuffer<double> n_period;  ///< Periodicity (stored as double, typically 1-3)
    CudaBuffer<int>    is_nci;    ///< NCI torsion flag (uses reduced cutoff)
    int n = 0;

    void upload(const std::vector<Dihedral>& standard,
                const std::vector<Dihedral>& extra,
                const std::vector<int>& atom_types,
                cudaStream_t stream = nullptr);
};

// ============================================================================
// InversionSoA: Inversion AoS → SoA
// ============================================================================
struct InversionSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k, idx_l;
    CudaBuffer<double> fc, omega0, C0, C1, C2;
    CudaBuffer<int>    potential_type;
    int n = 0;

    void upload(const std::vector<Inversion>& v, cudaStream_t stream = nullptr);
};

#endif // USE_CUDA
