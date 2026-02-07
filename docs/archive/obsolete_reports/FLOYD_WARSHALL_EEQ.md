# Floyd-Warshall Topological Distances for EEQ Phase 1

**Created**: 2025-12-28
**Status**: TODO - High Priority
**Impact**: EEQ topology charges currently 4-5× too large

## Problem

Curcuma's EEQ Phase 1 (topology charges) uses **geometric distances** between atoms, but XTB uses **topological distances** computed via Floyd-Warshall shortest-path algorithm.

**Result**: Topology charges are ~4-5× too large
- Curcuma: q(C) = -0.35, q(H) = +0.088 (CH₄)
- XTB Target: q(C) = -0.074, q(H) = +0.018

## XTB Implementation

**Reference**: `external/gfnff/src/gfnff_ini.f90:431-461`

### Algorithm Steps

1. **Initialize Bond Distance Matrix** (lines 431-442)
   ```fortran
   allocate(rabd(nat,nat), source=0.0e0_sp)
   rabd = rabd_cutoff  ! Large default value

   do i = 1,nat
     rabd(i,i) = 0.0
     do j = 1,topo%nb(20,i)  ! Loop over bonded neighbors
       k = topo%nb(j,i)
       rabd(k,i) = param%rad(at(i)) + param%rad(at(k))  ! Sum of covalent radii
       rabd(i,k) = rabd(k,i)
     end do
   end do
   ```

2. **Floyd-Warshall Shortest Path** (lines 443-453)
   ```fortran
   do k = 1,nat
     do i = 1,nat
       if (rabd(i,k) > gen%tdist_thr) cycle
       do j = 1,nat
         if (rabd(k,j) > gen%tdist_thr) cycle
         if (rabd(i,j) > (rabd(i,k)+rabd(k,j))) then
           rabd(i,j) = rabd(i,k)+rabd(k,j)  ! Update to shorter path
         end if
       end do
     end do
   end do
   ```

3. **Pack into Triangular Array** (lines 455-461)
   ```fortran
   do i = 1,nat
     do j = 1,i-1
       ij = lin(j,i)
       if (rabd(j,i) .gt. gen%tdist_thr) rabd(j,i) = rabd_cutoff
       rtmp(ij) = gen%rfgoed1 * rabd(j,i) / 0.52917726d0
     end do
   end do
   ```

4. **Use in EEQ Matrix** (gfnff_engrad.F90:1516-1527)
   ```fortran
   do j = 1,i-1
     ij = k+j
     gammij = 1./sqrt(topo%alpeeq(i)+topo%alpeeq(j))
     tmp = erf(gammij*r(ij))  ! r(ij) is topological distance!
     A(j,i) = tmp/r(ij)
     A(i,j) = A(j,i)
   end do
   ```

## Key Parameters

```fortran
! From gfnff_glob.f90 or similar
real(wp) :: rabd_cutoff = 1.0d8  ! Large value for unconnected atoms
real(wp) :: gen%tdist_thr        ! Threshold for topological distance
real(wp) :: gen%rfgoed1 = 1.0    ! Scaling factor for topological distances
```

## Implementation Plan for Curcuma

### 1. Add Topology Input to EEQSolver

EEQSolver needs access to molecular topology (bonded neighbor lists):

```cpp
struct TopologyInfo {
    std::vector<std::vector<int>> neighbor_lists;  // neighbor_lists[i] = atoms bonded to i
    std::vector<double> covalent_radii;            // From parameters
};

Vector calculateTopologyCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector& cn,
    const TopologyInfo& topology  // NEW parameter
);
```

### 2. Implement Floyd-Warshall Function

```cpp
// eeq_solver.cpp - private helper function
Matrix EEQSolver::computeTopologicalDistances(
    const std::vector<int>& atoms,
    const TopologyInfo& topology
) const
{
    const int natoms = atoms.size();
    const double RABD_CUTOFF = 1.0e8;
    const double TDIST_THR = 1.0e6;   // From XTB gen parameters
    const double RFGOED1 = 1.0;       // Scaling factor
    const double BOHR_TO_ANGSTROM = 0.52917726;

    // 1. Initialize with large values
    Matrix rabd = Matrix::Constant(natoms, natoms, RABD_CUTOFF);

    // 2. Set diagonal to zero
    for (int i = 0; i < natoms; ++i) {
        rabd(i, i) = 0.0;
    }

    // 3. Set bonded distances (sum of covalent radii)
    for (int i = 0; i < natoms; ++i) {
        double rad_i = topology.covalent_radii[i];
        for (int j : topology.neighbor_lists[i]) {
            double rad_j = topology.covalent_radii[j];
            rabd(i, j) = rad_i + rad_j;
            rabd(j, i) = rad_i + rad_j;  // Symmetric
        }
    }

    // 4. Floyd-Warshall shortest path
    for (int k = 0; k < natoms; ++k) {
        for (int i = 0; i < natoms; ++i) {
            if (rabd(i, k) > TDIST_THR) continue;
            for (int j = 0; j < natoms; ++j) {
                if (rabd(k, j) > TDIST_THR) continue;
                if (rabd(i, j) > rabd(i, k) + rabd(k, j)) {
                    rabd(i, j) = rabd(i, k) + rabd(k, j);
                }
            }
        }
    }

    // 5. Apply cutoff and scaling
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            if (rabd(i, j) > TDIST_THR) {
                rabd(i, j) = RABD_CUTOFF;
            }
            rabd(i, j) = RFGOED1 * rabd(i, j) / BOHR_TO_ANGSTROM;
        }
    }

    return rabd;
}
```

### 3. Use Topological Distances in EEQ Matrix

Replace current geometric distance calculation (eeq_solver.cpp:270-292):

```cpp
// Get topological distance matrix
Matrix topo_dist = computeTopologicalDistances(atoms, topology);

// Setup off-diagonal Coulomb matrix with topological distances
for (int i = 0; i < natoms; ++i) {
    for (int j = 0; j < i; ++j) {
        double r = topo_dist(i, j);  // Topological distance

        if (r > 1e6) {
            // Unconnected atoms - skip or use large distance
            A(i, j) = 0.0;
            A(j, i) = 0.0;
            continue;
        }

        double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
        double erf_gamma = std::erf(gammij * r);
        double coulomb = erf_gamma / r;

        A(i, j) = coulomb;
        A(j, i) = coulomb;
    }
}
```

### 4. Get Topology from GFNFF

GFNFF class already has topology information - pass it to EEQSolver:

```cpp
// gfnff_method.cpp
bool GFNFF::calculateTopologyCharges(TopologyInfo& topo_info) const
{
    // Build TopologyInfo structure
    TopologyInfo topology;
    topology.neighbor_lists.resize(m_atomcount);
    topology.covalent_radii.resize(m_atomcount);

    for (int i = 0; i < m_atomcount; ++i) {
        topology.neighbor_lists[i] = /* get from m_topology or bond list */;
        topology.covalent_radii[i] = /* get from parameter tables */;
    }

    // Pass topology to EEQSolver
    topo_info.topology_charges = m_eeq_solver->calculateTopologyCharges(
        m_atoms,
        geometry_bohr,
        total_charge,
        topo_info.cn,
        topology  // NEW
    );

    return true;
}
```

## Expected Impact

After Floyd-Warshall implementation:
- Topology charges should match XTB: q(C) ≈ -0.08, q(H) ≈ +0.02
- Dispersion error should drop to ~0% (currently ~2.6% for CH₄)
- All downstream calculations (Phase 2, final energies) will be correct

## Testing Plan

1. **Unit Test**: Floyd-Warshall with known topology
   - Linear chain: CH₃-CH₂-CH₃
   - Verify distances: 1-2 bonded, 1-3 through-bond = 2×bond length

2. **Integration Test**: CH₄ topology charges
   - Compare qa(C), qa(H) with XTB reference
   - Target: < 1% deviation

3. **Dispersion Test**: Full GFN-FF energy
   - CH₄ dispersion energy should match XTB within 0.1%

## References

- **XTB Source**: `external/gfnff/src/gfnff_ini.f90:431-461`
- **Floyd-Warshall Algorithm**: Standard graph shortest-path algorithm (O(N³))
- **Covalent Radii**: `gfnff_param.f90` parameter tables

---

**Priority**: HIGH - Blocking accurate GFN-FF dispersion calculation
