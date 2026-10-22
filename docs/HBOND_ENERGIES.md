# Per-triple GFN-FF hydrogen-bond energies

GFN-FF evaluates its hydrogen-bond term over a list of triples A-H...B (`GFNFF::getLastHBonds()`, the
`GFNFFHydrogenBond` entries). The workspace used to add every triple's energy straight into the sum.
It can now keep the single values, for example to count the triples below an energy limit as a hydrogen-bond
criterion.

```cpp
GFNFF* gfnff = energyCalculator.Interface()->gfnffInstance();
gfnff->setRecordHBondEnergies(true);          // from the next evaluation on
energyCalculator.CalculateEnergy(false);
std::vector<GFNFFHBondEnergy> e = gfnff->hydrogenBondEnergies();   // donor, hydrogen, acceptor, case_type, energy (Eh)
```

- The entries are in the order of `getLastHBonds()` and carry its atom indices (0-based). A triple outside its
  cut-off (`r_AB > r_cut`) has energy 0. The sum of the entries is the `HBond` term of `getEnergyDecomposition()`.
- Off by default: recording costs one double per triple and a write per triple and step. The partitions of
  `FFWorkspace` own disjoint index ranges, so no lock is needed; the vector is sized when recording is switched on
  and whenever the list is replaced (`updateHBonds`, `setInteractionLists`).
- CPU workspace only. When a GPU method owns the hydrogen-bond list (`setSkipHostHBWorkspaceUpdate`) the CPU list
  is not the evaluated one, and `hydrogenBondEnergies()` returns an empty vector. The GPU kernel `k_hbonds`
  computes the same value per thread; writing it out is the next step.
- The values are the terms of the list at the last evaluation. After a re-detection of the list the entries
  of the next evaluation apply to the new list.

## Status

🤖 AI-generated, ⚙️ machine-tested.

Checked (`test_hbond_energies`, 13 checks, constructed water dimers): nothing recorded unless asked, one entry per
triple in list order with the list's atoms, sum equals the `HBond` term to 1e-12 Eh, recording leaves the energy
unchanged, the bonded dimer's strongest triple is donor O to acceptor O, its energy shrinks to a tenth when
the dimer is stretched to 6 A, six stacked dimers give the same entries with 1 and 4 threads.

Not checked: that the 4-thread run really split the triples over several partitions (the systems are small),
re-detection during an MD run, react topology mode, systems with the case 2 to 4 triples (amides, phenols),
and the GPU path. Human production testing pending.
