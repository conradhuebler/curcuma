# Wall energy in the MD total energy - STATUS (Sep 13, 2026)

Decision (operator, 2026-09-12): **Etot = Epot + Ekin + Wall**, Epot stays the PURE
force-field energy. State: IMPLEMENTED, machine-tested, **uncommitted**.
Branch reactff2-llm, base 48aea7c4. Baseline binary (HEAD content) and changed binary
were both built in `build/` (build/curcuma was STALE before: its simplemd.cpp.o was from
Sep 12 14:36, commit 4601be27 touched the file at 17:39).

## Code: src/capabilities/simplemd.cpp, 5 assembly sites (all `= m_Epot + m_Ekin + m_wall_potential`)
- 2037 prepareRun(): wall never evaluated at t=0, so its term is 0 (= printed Wall column).
- 2731 step() rescue branch: value was STALE (previous geometry) -> `m_wall_potential = 0.0`
  first, because Energy() rebuilt a PURE-FF gradient with no wall force. m_rescue default off.
- 2777 step() periodic print: CURRENT - WallPotential() ran inside Verlet()/Rattle() this step.
- 2832 finalizeRun(): CURRENT - same last integrated step as Epot/Ekin (loop exits before the
  next integrator call). Sep 12 fix (stale final row) preserved.
- 4921 Results(): same convention; `wall_energy` key added only when m_wall_type != 0.
No site uses `+=`, so the wall enters exactly once per row; the Wall column and
m_average_wall_potential are untouched. Epot/m_Etot consumers: PrintStatus header rows,
a_aver_Etot (AverageQuantities), restart `average_Etot` - all consistent by construction.
SimpleMD::Results() has no in-tree caller (only the RMSD driver's Results() is called).

## Wall-free bit-identity (water, gfnff, max_time 1000, dt 0.5, wall_type none, 1 thread)
17-digit doubles from `input.final.json` (before == after, byte-identical file):
  average_Epot -0.32589396996387304   average_Ekin 0.0015976767466255731
  average_Etot -0.3255599651102289    average_T    336.3376693726637
md5 identical: input.trj.xyz, input_step_0.json, input_step_500.json.
stdout differs only in the echoed binary path + wall-clock stamps; input.topo.json only in
its `timestamp` field. => physics bit-identical.

## Falsifier: NVE, wall does work - 4 H (cli_simplemd_14 input), 6000 K, spheric r=3.0,
## dt 0.125 fs, 10 ps, print 0.1 ps, -md.rm_COM 0 -md.rmrottrans 0, seed 42 (both binaries)
Wall work is real: |Wall|max = 0.1782 Eh, nonzero in 99% of printed rows, T_max = 10063 K.
  before (Epot+Ekin)  : t=0 0.194211 -> t=5ps 0.135416 -> 10ps 0.156528
                        min 0.015185 max 0.194211 (swing 0.179 Eh = 92% of the total)
                        fitted slope -1.169e-04 +/- 1.6e-03 Eh/ps (fit meaningless, the
                        series swings over its whole range)
  after (Epot+Ekin+Wall): t=0 0.194211 -> t=5ps 0.194186 -> 10ps 0.194187
                        min 0.194163 max 0.194211 (spread 4.8e-05 Eh = 0.03 kcal/mol)
                        fitted slope -3.50e-10 +/- 9.3e-08 Eh/ps
  5 ps snapshot: before reported 0.135416 while 0.058770 Eh sat in the wall (Wall column);
  after reports 0.194186 = Epot(-0.000003)+Ekin(0.135419)+Wall(0.058770).
  Invariant max|Etot-(Epot+Ekin+Wall)| over every printed row = 1.0e-06 Eh (print rounding).
  finalizeRun row confirmed current: max_time 1000 / print_frequency 300, so the last row is
  emitted by finalizeRun and not by the periodic print. Its before value 0.024715 vs after
  0.194186 = Epot(-0.000007)+Ekin(0.024722)+Wall(0.169471) - same last step as Epot/Ekin.

## Secondary finding (pre-existing, NOT caused by this change; the trajectory is bit-identical)
The same 4 H + wall run with the DEFAULT `-md.rm_COM 100 -md.rmrottrans 1` loses ~2e-2 Eh/ps:
dt 0.125, wall on, rm on  -2.089e-02 +/- 4.9e-04 Eh/ps (maxdev 5.5e-02)
dt 0.125, wall on, rm off -2.534e-07 +/- 3.4e-07 Eh/ps (maxdev 2.4e-05)
dt 0.125, no wall, rm on  -3.167e-07 +/- 1.8e-07 Eh/ps
dt 0.125, no wall, rm off -5.656e-07 +/- 3.3e-07 Eh/ps
=> the periodic COM/rotation projection is only a sink while a wall re-creates that motion.
   cli_simplemd_14's stated setup (defaults + dt 0.25) therefore cannot show conservation
   before OR after: dt 0.25 with the wall leaves 2.9e-04 maxdev (-1.5e-05 +/- 3.9e-06 Eh/ps).
   At dt 0.125 it is 2.4e-05.

## Tests (build/, after binary)
ctest -R "cli_simplemd_": 22/22 pass, 93.30 s (slowest cli_simplemd_10 69.4 s).
ctest -R "gfnff": 65/65 pass, 94.88 s, 0 failures.
No simplemd test reads Etot with a wall: 18 uses r[5] but runs wall_type none; 16/17 read
T only (r[7]); 14 greps REACT lines; 19 compares mean Epot (no wall). No test needs a decision.

## Reproduce
binaries: .../scratchpad/wall/curcuma_head (HEAD) and curcuma_wall (changed);
run_measure.sh / probe.sh / probe2.sh / analyse.py sit next to them.
`cd build && make -j4`, check the EXIT STATUS (not a grep for "error:"). `cp` is aliased to -i.
