# ConfSearch Restart / Checkpointing (`-restart`)

> 🤖 AI-generated, machine-tested (Jun 2026). Human production testing pending.

ConfSearch can be interrupted (crash, time limit, kill) and resumed from a
self-contained checkpoint, so a long temperature-schedule search does not have
to start over.

## Usage

```bash
# first run (writes checkpoints as it goes)
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2 -restart

# ... interrupted (Ctrl-C / kill / time limit) ...

# resume: same command, the checkpoint is picked up automatically
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2 -restart
```

`-restart` both **enables checkpoint writing** and **resumes** if a checkpoint
is present. Without `-restart`, ConfSearch behaves exactly as before (no
checkpoint files, no resume).

## The checkpoint file

`<basename>.confsearch.restart.json` is written **into the BMT output dir and
copied back to the start directory (CWD)** under that stable name. Because BMT
creates a fresh timestamped directory every run, the CWD copy is what a resume
run reads — it survives the timestamp change.

It is **self-contained** (everything needed to resume is inside the JSON):

- the full **shared bias pool** — every structure's geometry plus its
  `counter`/`index`/`persistent`/`energy` metadata, restored exactly;
- the accumulated **cumulative conformer pool** (completed cycles);
- the next-cycle **seeds** (`m_in_stack`);
- **energy progress**: exploration `global_min`, `best`, `initial`;
- the **symmetry-permutation cache** and the **topology reference**;
- the **md/opt methods** and the schedule position (`next_T`, cycle, phase).

All frames share one atomic-number list, so each structure is stored as a flat
`x|y|z|...` geometry string + energy (no XYZ parser needed on load).

## When checkpoints are written

- **`post_md`** — after every MD phase (within the cycle). The expensive MD
  exploration and the grown bias pool are now safe.
- **`post_cycle`** — after every completed temperature cycle. Cumulative pool,
  seeds and energies for that cycle are final.

## What resume does

1. Reads the CWD checkpoint, validates that the md/opt methods match the current
   run (otherwise it warns and starts fresh).
2. **Skips the initial pre-optimisation.**
3. Restores the bias pool (`SharedBiasPool::restoreStructures`), rebuilds the
   cumulative file in the new BMT dir, restores seeds/energies/permutations/topo.
4. Continues the temperature schedule from `next_T`. The interrupted cycle
   re-enters at `post_md`: its MD is **not** re-run (the restored bias pool is
   re-exported to `*.bias.xyz`), and Phase 2..4 proceed normally.

## Known limitations (not yet human-tested)

- **Granularity:** a resume re-enters the interrupted cycle at `post_md`, so that
  single cycle's optimisation phases (including an expensive `opt_method`
  Phase 3b, e.g. r2scan) are re-run. Finer mid-cycle resume (skip a completed
  `post_filter`/`post_refine`) is already supported by the serializer but is not
  yet used by the resume loop.
- **`stop` file:** the graceful-halt `stop` file is consumed by the inner
  SimpleMD during the MD phase, so ConfSearch's own `CheckStop` rarely sees it.
  The reliable way to interrupt is to **kill the process** (Ctrl-C / SIGTERM /
  time limit) — a valid checkpoint is on disk after every phase.
- **Method change:** resuming with different `md_method`/`opt_method` than the
  checkpoint is refused (fresh start) to avoid mixing incompatible state.

## Flag routing note

Since July 2026 ConfSearch registers its own parameters (including `restart`) in
the ParameterRegistry, so the flat `-restart` flag stays in
`controller["confsearch"]` and is read directly. `-confsearch.restart` works too.

Before that migration ConfSearch used a static JSON parameter set and owned no
names, so `-restart` was silently auto-routed to `controller["confscan"]` (the
only module registering the name) and had to be read back from there.
