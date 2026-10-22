# Per-triple GFN-FF hydrogen-bond energies

GFN-FF evaluates its hydrogen-bond term over a list of triples A-H...B (`GFNFF::getLastHBonds()`, the
`GFNFFHydrogenBond` entries). The workspace used to add every triple's energy straight into the sum.
It can now keep the single values, for example to count the triples below an energy limit as a hydrogen-bond
criterion.

```cpp
GFNFF* gfnff = energyCalculator.Interface()->gfnffInstance();   // CPU and GPU methods alike
gfnff->setRecordHBondEnergies(true);          // from the next evaluation on
energyCalculator.CalculateEnergy(false);
std::vector<GFNFFHBondEnergy> e = gfnff->hydrogenBondEnergies();   // donor, hydrogen, acceptor, case_type, energy (Eh)
```

- The entries are in the order of `getLastHBonds()` and carry its atom indices (0-based). A triple outside its
  cut-off (`r_AB > r_cut`) has energy 0. The sum of the entries is the `HBond` term of `getEnergyDecomposition()`
  (CPU; the GPU decomposition reports only the total).
- Off by default: recording costs one double per triple and a write per triple and step. The partitions of
  `FFWorkspace` own disjoint index ranges, so no lock is needed; the vector is sized when recording is switched on
  and whenever the list is replaced (`updateHBonds`, `setInteractionLists`).
- The values are the terms of the list at the last evaluation. After a re-detection of the list the entries
  of the next evaluation apply to the new list.
- The flag lives on the `GFNFF` object. When `EnergyCalculator` replaces a GPU method by the CPU one (counted GPU
  fallback, e.g. a charged system with several fragments), take `gfnffInstance()` again and switch it on there.

## GPU (CUDA, ROCm)

- The device evaluates the term in `k_hbonds` (`cuda/gfnff_kernels.cu`, `rocm/gfnff_rocm.hip`). With recording on,
  each thread writes its triple's energy into a device array of one double per triple (`hb_triple_energy`, allocated
  when recording is first switched on and grown with the list); the value is the one the block reduction adds to the
  device HB term. With recording off the kernel gets `nullptr` and makes no extra memory access.
- The GPU wrapper (`GFNFFGpuMethodImpl`, `qm_methods/gfnff_gpu_method_impl.h`) registers its device workspace with
  `GFNFF::setHBondEnergySource()`, so `setRecordHBondEnergies()` and `hydrogenBondEnergies()` on `gfnffInstance()`
  reach the device. `hydrogenBondEnergies()` then synchronises the device and downloads the energies together with the
  atom indices and case types of the device list (`FFWorkspaceGPU::downloadHBondEnergies`); nothing is copied unless
  it is called. The wrapper unregisters itself in its destructor.
- Order: that of `getLastHBonds()`, as on the CPU. ROCm stores the triples grouped by case type
  (`hbCasePermutation`); the download puts them back into list order.
- Switching the recording, or a list that outgrows the device array, drops a captured Phase-1 graph, because the
  array is a kernel argument.
- Vulkan has no GFN-FF port (`curcuma_vulkan_create_gfnff` returns nullptr): `-gpu vulkan` runs GFN-FF on the CPU,
  which records as above. The HB term has no single-precision kernel (only the halogen bond has `k_xbonds_mixed`).

## Status

🤖 AI-generated. CPU path ⚙️ machine-tested; GPU path compiled only.

Checked (`test_hbond_energies`, 12 checks, constructed water dimers): nothing recorded unless asked, one entry per
triple in list order with the list's atoms, sum equals the `HBond` term to 1e-12 Eh, recording leaves the energy
unchanged, the bonded dimer's strongest triple is donor O to acceptor O, its energy shrinks to a tenth when
the dimer is stretched to 6 A, six stacked dimers give the same entries with 1 and 4 threads.

GPU: `test_hbond_energies_gpu cuda|rocm` (ctest `hbond_energies_cuda`, `hbond_energies_rocm`) compares, for the
registry molecules `H2O_dimer` and `monosaccharide`, the GPU entries with the CPU ones (same triples, per-triple
energy within 1e-8 Eh, sum within 1e-6 Eh of the CPU `HBond` term), and checks that recording leaves the GPU energy
unchanged and that nothing is returned before recording is switched on. Without a device it prints SKIP and exits 0.
The CUDA plugin (CUDA 13.3, sm_86) and the ROCm plugin (hipcc, gfx1150) were compiled, the test was built in the CUDA
tree and ran only its SKIP path (no device on that machine); neither the comparison nor the kernels have run on GPU
hardware, so no GPU value is verified. With recording off the arithmetic of `k_hbonds` is
unchanged as far as the machine code shows: the sm_86 SASS and the gfx1150 assembly of the kernel before and after
the change have the same number of each FP64 instruction (SASS DFMA/DMUL/DADD 1137/531/284, gfx1150
v_fma/v_mul/v_add_f64 1057/544/502); the differences are instruction order, registers and the guarded store.

Not checked: that the 4-thread run really split the triples over several partitions (the systems are small),
re-detection during an MD run, react topology mode, systems with the case 2 to 4 triples (amides, phenols),
anything on GPU hardware (values, ROCm permutation, graph capture with recording, multi-GPU). Human production
testing pending.
