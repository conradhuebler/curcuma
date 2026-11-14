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

## Logging during MD (`hb_log_frequency`)

Opt-in. `-hb_log_frequency N` (default 0 = off) writes every N steps one row to `<basename>.hbonds.csv`:
`step, time_fs, hbonds_total, a_donates, b_donates, energy_sum_kj_mol`. A triple counts when its energy is below
`-hb_log_energy_limit` (kJ/mol, negative = attractive, default 0 = every attractive triple). `-hb_log_selection`
and `-hb_log_selection_b` (selection grammar, `!F1` = all but the first fragment) choose the groups exactly as in
`-trajprops` ([TRAJPROPS.md](TRAJPROPS.md)); empty counts every triple.

```
curcuma -md system.xyz -method gfnff -hb_log_frequency 100 -hb_log_energy_limit -8 -hb_log_selection F1 -hb_log_selection_b '!F1'
```

- Off means off: the only code on the step path is a pointer test. With it on, the force field records the triple
  energies (one write per triple and step) and nothing else changes; the trajectory is identical to the run
  without the option (`test_md_hbond_log`, separate processes, NVE and a fixed seed).
- GFN-FF only. With another method, or a selection that matches nothing, the run says so and continues without
  the log; it does not stop.
- The row reflects the force evaluation of the step it is written in. The logging code is one function
  (`SimpleMD::LogHydrogenBonds`) over `curcuma::classifyHBondEnergies()`, the same counting as `-trajprops`.
- GPU runs: the energies are downloaded from the device at each logged step (see the GPU section); not run on
  hardware.
- The limit has no standard value. Choose it from the histogram of `-trajprops -hb_criterion energy` for a
  trajectory of the same system.

Checked (`test_md_hbond_log`, 16 checks, three stacked water dimers, 60 fs): with the option off no log file; the
settings reach `SimpleMD`; two runs without the option agree (the control); with the log on the trajectory equals
the reference run, also with limit 0; six rows for 60 steps every 10th; total = a_donates + b_donates; limit 0
never counts fewer bonds than limit -5 on the same trajectory; with `uff` the log is refused and the run goes on.
On the CLI a run of the same dimers gave three bonds at limit -5 with a triple-energy sum of about -35 kJ/mol.

Reference comparison (`scripts/md_reference_compare.py`): the reference binary is curcuma master at the commit this work
started from (8edb8d15, built from a `git archive` of it with the same CMake options); the new binary has all
branches of this work (per-triple recording, GPU code, `trajprops`, the log). Three stacked water dimers, 100 fs NVE,
and a 653-atom lattice of waters, 50 fs with CSVR, each with one thread and a fixed seed: the reference agrees with
itself, the new binary without options equals the reference, and with `-hb_log_frequency 10` it equals the reference
as well (trajectory md5, comment lines compared without the version string). Both binaries are CPU builds
(`-DUSE_CUDA=OFF -DUSE_ROCM=OFF`): the GPU code was compared only in what it compiles to, not in a run.

Not checked: a solvated polymer, runs of days, restart/rescue (the log file is opened afresh at each start),
`topology react`, the GPU path, a restart appending to the old log. Two `SimpleMD` runs inside one process do not
start from the same random state (the second run's initial velocities differ, with the same seed), which is why
the test starts one process per run; that is independent of this option and not investigated.
