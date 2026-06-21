#include "hip/hip_runtime.h"
/*
 * <GFN-FF GPU Struct-of-Arrays (SoA) Data Structures>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (March 2026): SoA layout for GPU-coalesced memory access.
 * Converts GFN-FF AoS parameter structs to GPU-friendly SoA with RAII hipMalloc.
 *
 * Reference: Spicher/Grimme J. Chem. Theory Comput. 2020 (GFN-FF)
 */

#pragma once

#ifdef USE_ROCM_GFNFF

#include <hip/hip_runtime.h>
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
struct GFNFFSTorsion;
struct ATMTriple;
struct GFNFFBatmTriple;
struct GFNFFHalogenBond;
struct GFNFFHydrogenBond;

// ============================================================================
// CudaBuffer<T>: RAII wrapper around hipMalloc / hipFree
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

    /// Allocate count elements only if the current buffer is smaller — reuse the
    /// existing allocation across calls (e.g. MD/opt geometry steps) instead of
    /// free+realloc churn. Claude Generated.
    void ensure(int count) {
        if (count > n) alloc(count);
    }

    /// Allocate count elements on GPU (frees previous allocation)
    void alloc(int count) {
        free();
        if (count <= 0) return;
        hipError_t err = hipMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T));
        if (err != hipSuccess)
            throw std::runtime_error(std::string("hipMalloc failed: ") + hipGetErrorString(err));
        n = count;
    }

    /// Upload count elements from host to GPU (auto-allocates if needed)
    /// Uses synchronous hipMemcpy because callers often pass temporary host
    /// buffers (std::vector locals in SoA::upload) that are destroyed before
    /// any stream synchronisation point.  Async copy from freed memory was
    /// the root cause of the "double free or corruption" crash.
    void upload(const T* host_data, int count, hipStream_t /*stream*/ = nullptr) {
        if (count <= 0) return;
        if (n < count) alloc(count);
        hipError_t err = hipMemcpy(ptr, host_data, count * sizeof(T),
                                     hipMemcpyHostToDevice);
        if (err != hipSuccess)
            throw std::runtime_error(std::string("hipMemcpy H2D failed: ") + hipGetErrorString(err));
    }

    /// Upload from std::vector
    void upload(const std::vector<T>& v, hipStream_t stream = nullptr) {
        upload(v.data(), static_cast<int>(v.size()), stream);
    }

    /// Download count elements from GPU to host
    void download(T* host_data, int count, hipStream_t stream = nullptr) const {
        if (count <= 0 || !ptr) return;
        hipError_t err = hipMemcpyAsync(host_data, ptr, count * sizeof(T),
                                          hipMemcpyDeviceToHost, stream);
        if (err != hipSuccess)
            throw std::runtime_error(std::string("hipMemcpyAsync D2H failed: ") + hipGetErrorString(err));
    }

    /// Zero-fill count elements on GPU
    void zero(int count, hipStream_t stream = nullptr) {
        if (count <= 0 || !ptr) return;
        hipMemsetAsync(ptr, 0, count * sizeof(T), stream);
    }

    void free() {
        if (ptr) { hipFree(ptr); ptr = nullptr; n = 0; }
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
    CudaBuffer<double> dc6dcn_ij;  ///< [n] dC6(i,j)/dCN(i) per pair (updated per step)
    CudaBuffer<double> dc6dcn_ji;  ///< [n] dC6(i,j)/dCN(j) per pair (updated per step)
    int n = 0;

    void upload(const std::vector<GFNFFDispersion>& v, hipStream_t stream = nullptr);
};

// ============================================================================
// RepulsionSoA: GFNFFRepulsion AoS → SoA (bonded and nonbonded share struct)
// ============================================================================
struct RepulsionSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> alpha, repab, r_cut;
    int n = 0;

    void upload(const std::vector<GFNFFRepulsion>& v, hipStream_t stream = nullptr);
};

// ============================================================================
// CoulombSoA: GFNFFCoulomb AoS → SoA (charges come from per-step d_charges)
// ============================================================================
struct CoulombSoA {
    CudaBuffer<int>    idx_i, idx_j;
    CudaBuffer<double> gamma_ij, r_cut;
    int n = 0;

    // ROCm gather path (Claude Generated June 2026): per-atom CSR adjacency built from the
    // pair list. The pair k_coulomb does 6 FP64 atomicAdd per pair (~6M on polymer); on RDNA
    // (slow/emulated FP64 atomics) that is the dominant MD cost (185 ms). k_coulomb_gather
    // sums per atom (no in-loop atomics, only 3 atomicAdd/atom) — bit-identical gamma/rcut.
    CudaBuffer<int>    csr_off;      ///< [N+1] row offsets
    CudaBuffer<int>    csr_partner;  ///< [2*n] partner atom per directed edge
    CudaBuffer<double> csr_gamma;    ///< [2*n] gamma_ij per directed edge
    CudaBuffer<double> csr_rcut;     ///< [2*n] r_cut per directed edge
    int csr_N = 0;

    void upload(const std::vector<GFNFFCoulomb>& v, int natoms, hipStream_t stream = nullptr);
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
    CudaBuffer<int>    nr_hb;      ///< HB count for bond (>= 1 triggers alpha modification)
    CudaBuffer<double> hb_cn_H;    ///< HB coordination number of H atom
    CudaBuffer<int>    hb_H_atom;  ///< Index of H atom in bond (-1 if no HB)
    int n = 0;

    void upload(const std::vector<Bond>& v, const std::vector<int>& atom_types,
                hipStream_t stream = nullptr);
};

// ============================================================================
// HBAlphaSoA: HB neighbor pairs (H, B) for alpha-modulation CN chain-rule gradient
// Claude Generated (March 2026): Flattened from BondHBEntry B_atoms
// ============================================================================
struct HBAlphaSoA {
    CudaBuffer<int>    idx_H, idx_B;    ///< H and B atom indices per pair
    CudaBuffer<double> rcov_sum;        ///< Scaled covalent radius sum per pair
    int n = 0;

    void upload(const std::vector<int>& h_idx, const std::vector<int>& b_idx,
                const std::vector<double>& rcov, hipStream_t stream = nullptr);
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
                hipStream_t stream = nullptr);
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
                hipStream_t stream = nullptr);
};

// ============================================================================
// InversionSoA: Inversion AoS → SoA
// ============================================================================
struct InversionSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k, idx_l;
    CudaBuffer<int>    ati, atj, atk, atl;  ///< Atomic numbers for rcov damping lookup
    CudaBuffer<double> fc, omega0, C0, C1, C2;
    CudaBuffer<int>    potential_type;
    int n = 0;

    void upload(const std::vector<Inversion>& v,
                const std::vector<int>& atom_types,
                hipStream_t stream = nullptr);
};

// ============================================================================
// STorsionSoA: Triple bond torsions (simple cos(2φ))
// Claude Generated (March 2026): GPU port of calcSTorsions
// ============================================================================
struct STorsionSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k, idx_l;
    CudaBuffer<double> erefhalf;
    int n = 0;

    void upload(const std::vector<GFNFFSTorsion>& v, hipStream_t stream = nullptr);
};

// ============================================================================
// BATMSoA: Bonded ATM 3-body with topology charge scaling
// Claude Generated (March 2026): GPU port of calcBATM
// ============================================================================
struct BATMSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k;
    CudaBuffer<double> zb3atm_i, zb3atm_j, zb3atm_k;
    int n = 0;

    void upload(const std::vector<GFNFFBatmTriple>& v, hipStream_t stream = nullptr);
};

// ============================================================================
// ATMSoA: Axilrod-Teller-Muto 3-body dispersion
// Claude Generated (March 2026): GPU port of calcATM + calcATMGradient
// ============================================================================
struct ATMSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k;
    CudaBuffer<int>    ati, atj, atk;  ///< Atomic numbers for rcov lookup
    CudaBuffer<double> C6_ij, C6_ik, C6_jk;
    CudaBuffer<double> s9, a1, a2, alp, triple_scale;
    int n = 0;

    void upload(const std::vector<ATMTriple>& v,
                const std::vector<int>& atom_types,
                hipStream_t stream = nullptr);
};

// ============================================================================
// XBondSoA: Halogen bonds (3-body A-X...B)
// Claude Generated (March 2026): GPU port of calcHalogenBonds
// ============================================================================
struct XBondSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k;   ///< A, X, B
    CudaBuffer<int>    elem_A, elem_B;         ///< Atomic numbers for vdW radii
    CudaBuffer<double> q_X, q_B, acidity_X;
    CudaBuffer<double> r_cut;
    int n = 0;

    void upload(const std::vector<GFNFFHalogenBond>& v,
                const std::vector<int>& atom_types,
                hipStream_t stream = nullptr);
};

// ============================================================================
// HBondSoA: Hydrogen bonds (3-body A-H...B) with packed neighbor lists
// Claude Generated (March 2026): GPU port of calcHydrogenBonds
// Variable-length neighbors packed into flat arrays with per-HB offsets.
// ============================================================================
struct HBondSoA {
    CudaBuffer<int>    idx_i, idx_j, idx_k;   ///< A, H, B
    CudaBuffer<int>    elem_A, elem_B;         ///< Atomic numbers for vdW radii
    CudaBuffer<double> q_H, q_A, q_B;
    CudaBuffer<double> basicity_A, basicity_B;
    CudaBuffer<double> acidity_A, acidity_B;
    CudaBuffer<double> r_cut;
    CudaBuffer<int>    case_type;

    // Packed neighbor indices for case 2+ (neighbors_B)
    CudaBuffer<int>    nb_B_offset, nb_B_count; ///< [n] offset and count into nb_B_flat
    CudaBuffer<int>    nb_B_flat;               ///< Flat array of all neighbor_B indices

    // Case 3: acceptor parent index
    CudaBuffer<int>    acceptor_parent;         ///< [n] acceptor_parent_index (-1 if N/A)

    // Case 3: packed neighbors_C for torsion calculation
    CudaBuffer<int>    nb_C_offset, nb_C_count; ///< [n] offset and count into nb_C_flat
    CudaBuffer<int>    nb_C_flat;               ///< Flat array of all neighbor_C indices

    // Case 4: repz values for LP distance
    CudaBuffer<double> repz_B;                  ///< [n] GFNFFParameters::repz[Z_B-1]

    int n = 0;
    int total_nb_B = 0;  ///< Total packed neighbor_B count
    int total_nb_C = 0;  ///< Total packed neighbor_C count

    void upload(const std::vector<GFNFFHydrogenBond>& v,
                const std::vector<int>& atom_types,
                hipStream_t stream = nullptr);
};

// ============================================================================
// CoordSoA: Per-step atomic coordinates in Structure-of-Arrays layout
// Claude Generated (March 2026): Replaces d_coords[3*N] AoS buffer.
// SoA layout gives fully coalesced warp reads: x[0..31] fits in one cache line,
// giving ~3× better bandwidth vs stride-3 AoS for N > 100 atoms.
// Gradients remain AoS (add_grad uses atomicAdd scatter — no SoA benefit).
// ============================================================================
struct CoordSoA {
    CudaBuffer<double> d_x;  ///< [N] x-coordinates (Bohr)
    CudaBuffer<double> d_y;  ///< [N] y-coordinates (Bohr)
    CudaBuffer<double> d_z;  ///< [N] z-coordinates (Bohr)
    int N = 0;

    void alloc(int n) { d_x.alloc(n); d_y.alloc(n); d_z.alloc(n); N = n; }
    bool empty() const { return N == 0 || d_x.empty(); }
};

// ============================================================================
// RefCoordSoA: Reference geometry for topology displacement check
// Claude Generated (March 2026): Parallel SoA layout to CoordSoA.
// Updated only on topology rebuild; compared against CoordSoA each step.
// ============================================================================
struct RefCoordSoA {
    CudaBuffer<double> d_rx;  ///< [N] reference x-coordinates (Bohr)
    CudaBuffer<double> d_ry;  ///< [N] reference y-coordinates (Bohr)
    CudaBuffer<double> d_rz;  ///< [N] reference z-coordinates (Bohr)
    int N = 0;

    void alloc(int n) { d_rx.alloc(n); d_ry.alloc(n); d_rz.alloc(n); N = n; }
    bool empty() const { return N == 0 || d_rx.empty(); }
};

#endif // USE_ROCM_GFNFF
