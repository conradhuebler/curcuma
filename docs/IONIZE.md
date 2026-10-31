# Counter-ion placement (`-ionize`)

Neutralises a solvated system: solvent molecules are replaced by monovalent ions until the total charge
equals `target_charge`. The atom count per volume stays what the solvent packing gave (the approach of
GROMACS `genion`). Code: `src/capabilities/ionize.{h,cpp}`, test: `test_cases/test_ionize.cpp`.

```
curcuma -ionize solvated.xyz -charge 15              # 15 Cl-, result solvated_ionized.xyz in the BMT directory
curcuma -ionize solvated.xyz -charge -4 -cation K    # 4 K+
```

## Behaviour

- `charge` is the total charge of the input structure. It is not derived from the geometry; the caller knows
  the protonation state. `|charge - target_charge|` ions are added: anions (`anion`, default Cl) for a positive
  system, cations (`cation`, default Na) for a negative one. All ions are monovalent.
- Solvent molecules are the connected fragments (`Molecule::GetFragments`) whose composition equals the
  `solvent` formula (default `H2O`). Everything else, ions already present included, counts as solute.
- A solvent molecule is a candidate if its heaviest atom is at least `min_solute_distance` (6 A) from every
  solute atom. Candidates are tried in an order fixed by `seed`; one is taken if it is at least
  `min_ion_distance` (7 A) from the ions chosen so far. The ion sits at the position of the heaviest atom.
- With a unit cell (`Molecule::hasPBC`) distances use the minimum image convention.
- Output: the molecule without the replaced solvent, ions appended last, cell and charge carried over.
  Fewer fitting ions than requested is an error; no partial structure is written.
- Results (`Results()`): `output`, `ion`, `ions_added`, `solvent_molecules`, `candidates`, `final_charge`,
  `atoms`, `replaced_atoms`.

## Status

🤖 AI-generated, ⚙️ machine-tested.

Checked: `test_ionize` (20 checks, a constructed box of 999 water molecules around a methane) and
`curcuma -ionize` on a 653-atom box of the same kind (expected atom count, exit code 1 for an impossible
request). The unit test covers counts of
atoms and elements, ions last, minimum distances to the solute and between ions, solvent molecules left
whole, same seed same result, cation/anion choice, neutral input unchanged, error for too many ions and
for a missing solvent, minimum-image distances in a cubic cell.

Not checked: real solvent boxes from a packing tool or an MD run, large systems (timing), solvents
other than water, ions of charge other than +-1, triclinic cells, ions with a non-water solvent that has
internal fragments, how the result behaves in a subsequent simulation. Human production testing pending.

Not implemented: placement by electrostatic potential, a salt concentration (ions beyond neutralisation),
re-equilibration around the new ions.
