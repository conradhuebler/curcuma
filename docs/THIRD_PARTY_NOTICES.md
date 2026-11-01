# Third-party notices

Data from other projects that is compiled into the curcuma binary, with its license. Code libraries under `external/`
carry their own license files and are not listed here.

## Basis Set Exchange (built-in basis sets of the native QM engine)

**What is derived from it.** The GENERATED headers `src/core/energy_calculators/qm_methods/qm_basis_data_*.h` and
`qm_basis_tables.h` (written by `scripts/qm_basis_export_tables.py`, BSE repository commit
`4adaf1372c7101620ca1a9f3130be9ae97fb8f30`, master of 2026-08-07, downloaded 2026-11-01; every generated header carries
source, commit and the BSD-3 text below). Stage 1 (Nov 2026) built in five sets for H-Ne; stage 2 (Nov 2026) every def2
set of the BSE for every element it has:

| Built-in sets | Elements | BSE source (`basis_set_exchange/tests/sources/`) |
|---|---|---|
| def2-SVP, def2-SV(P), def2-TZVP, def2-TZVPP, def2-QZVP, def2-QZVPP | H-Rn (86) | `<name>.1.tm.bz2` |
| def2-SVPD, def2-TZVPD, def2-TZVPPD, def2-QZVPD, def2-QZVPPD | H-La, Hf-Rn (72) | `<name>.1.tm.bz2` |
| def2-mTZVP, def2-mTZVPP | H-Rn, Th-Lr (100) | `def2-mtzvp(p).1.tm.bz2` |
| def2-universal-JFIT, def2-universal-JKFIT, def2-SV(P)-JKFIT | H-Rn (86) | `<name>.1.tm.bz2` |
| def2-mTZVPP-RIJ | H-Ar, Fe (19) | `def2-mtzvpp-rij.1.tm.bz2` |
| def2-{SVP,SV(P),SVPD,TZVP,TZVPP,TZVPD,TZVPPD,QZVPP,QZVPPD}-RIFIT | H-La, Hf-Rn (72) | `<name>.1.tm.bz2` |
| def2-QZVP-RIFIT | H-Ca, Ga-Sr, In-Ba, Tl-Rn (42) | `def2-qzvp-rifit.1.tm.bz2` |
| def2-ECP (effective core potentials, raw text, unused by the engine) | Rb-Rn (50) | `def2-ecp.1.tm.bz2` |
| MINIX (MINIS part) | H-Ne | `basis_set_exchange/data/huzinaga/MINI-scaled.0.json` ('Scaled MINI') |

def2 data: every exponent/coefficient token unchanged except a Fortran `D` exponent marker written as `E`; shell order as
in the `.tm` source. The generator checks each set against the BSE json data composed the way BSE composes it
(`data/<name>.1.table.json` -> element json -> component files): same elements, same numbers as IEEE doubles for every
primitive; it reports where only the spelling of a number (`.80` vs `0.80`) or the shell order differs (shell order differs
for many elements, e.g. H-Ne of the diffuse sets, because BSE composes base + add files). The ECP sections inside the
orbital `.tm` files equal `def2-ecp.1.tm.bz2`; BSE's def2-mTZVP/-mTZVPP carry no ECP for any element.

**Turbomole origin.** The BSE file headers of the def2 sets state `Version: 1  (Data from Turbomole 7.3)`; the BSE notes
name Florian Weigend as original contributor (PI F. Weigend, R. Ahlrichs). The owner decided that the BSE BSD-3 notice
suffices; the Turbomole terms themselves were **not verified**. For MINI-scaled the BSE notes name David Feller as original
contributor (15 Jan 2007).

**Shipped files.** The files in `src/core/energy_calculators/qm_methods/` stay in the repository and take precedence over
the built-in copies when found in the basis search path. `def2-SV(P).dat` (86 elements) has the same numbers as BSE for
every element; for 37 elements beyond Ne (Sc-Zn, Rb, In-Xe, Ce-Lu, Tl-Rn) its shell order differs from the `.tm` source
(not a number difference; an earlier note of "33 differing elements beyond Ar" was this order difference).

**MINIX.** Composition from R. Sure, S. Grimme, J. Comput. Chem. 2013, 34, 1672-1685, Table 1:

| Elements | MINIX |
|---|---|
| H-He, B-Ne | MINIS |
| Li-Be | MINIS + 1(p) |
| Na-Mg | MINIS + 1(p) |
| Al-Ar | MINIS + 1(d) |
| K-Zn | SV |
| Ga-Kr | SVP |
| Rb-Xe | def2-SV(P) with ECP |

Built in: H-Ne only (Rb-Xe needs def2-SV(P) plus the def2 ECP in the SCF, which the engine does not have). Per element:

| Element | Exponents | Contraction coefficients |
|---|---|---|
| H, He, B-Ne | BSE MINI-scaled | BSE MINI-scaled, each contraction renormalised to a unit-norm contracted function, 10 decimals |
| Li, Be (s shells) | BSE MINI-scaled | as above |
| Li, Be (added p shell: 0.1 / 0.25, coefficient 1) | ORCA 6.1 `orca_exportbasis -b minix` (not in BSE) | ORCA 6.1 |

The renormalised BSE coefficients equal the ORCA 6.1 export and the shipped `MINIX.dat` token for token (raw BSE
coefficients differ from them by up to 1.1e-6 relative, a pure normalisation factor per contraction).

**ORCA 6.1 cross-check** (`orca_exportbasis -b <name>`, H-Ne, measured 2026-11-01 with
`scripts/qm_basis_export_tables.py --compare-orca DIR`; ORCA prints 10 decimals and renormalises contractions, so
coefficients are compared after normalising both sides; shells matched by L, length and exponents):

| Set (ORCA name) | shells | exponents, max rel. diff | coefficients (normalised), max rel. diff | max abs. diff (ORCA scale) |
|---|---|---|---|---|
| def2-SVP, def2-SV(P) | 52, 51 | 5.1e-10 | 1.5e-08 | 8.6e-11 |
| def2-SVPD | 71 | 1.6e-09 | 1.5e-08 | 8.6e-11 |
| def2-TZVP, -TZVPP, -TZVPD, -TZVPPD | 91-115 | <= 1.7e-09 | 2.7e-07 | 1.9e-10 |
| def2-QZVP, -QZVPP, -QZVPD, -QZVPPD | 149-166 | 2.0e-09 | 1.7e-06 | 2.4e-10 |
| def2-mTZVP, def2-mTZVPP | 74, 78 | 1.1e-09 | 2.7e-07 | 1.9e-10 |
| def2-universal-JFIT (def2/J), -JKFIT (def2/JK) | 130, 199 | 0 | <= 3.5e-09 | 5.5e-11 |
| def2-mTZVPP-RIJ (def2-mTZVPP/J) | 103 | 1.5e-10 | 4.7e-10 | 5.3e-11 |
| def2-<X>-RIFIT (def2-<X>/C), X = SVP, SVPD, TZVP, TZVPP, TZVPD, TZVPPD, QZVPP, QZVPPD | 129-313 | 0 | <= 1.3e-09 | 5.5e-11 |
| MINIX | 26 | 0 | 0 | 0 |

All differences are at the 1e-10 print resolution of ORCA (the larger relative values belong to small coefficients). Not
compared (ORCA 6.1 has no such name): def2-SV(P)-JKFIT, def2-SV(P)-RIFIT, def2-QZVP-RIFIT, def2-ECP. Elements beyond Ne
were not compared with ORCA.

**License of the Basis Set Exchange** (BSD-3-Clause, https://github.com/MolSSI-BSE/basis_set_exchange, file `LICENSE`),
reproduced verbatim:

```
Copyright (c) 2020 The Molecular Sciences Software Institute, Virginia Tech

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

**References.** def2: F. Weigend, R. Ahlrichs, Phys. Chem. Chem. Phys. 7, 3297-3305 (2005); A. Schaefer, H. Horn,
R. Ahlrichs, J. Chem. Phys. 97, 2571-2577 (1992); QZVP: F. Weigend, F. Furche, R. Ahlrichs, J. Chem. Phys. 119,
12753-12762 (2003); diffuse sets: D. Rappoport, F. Furche, J. Chem. Phys. 133, 134105 (2010); lanthanides: R. Gulde,
P. Pollak, F. Weigend, J. Chem. Theory Comput. 8, 4062-4068 (2012). def2-mTZVP: J. G. Brandenburg et al., J. Chem. Phys.
148, 064104 (2018); def2-mTZVPP and -RIJ: S. Grimme et al., J. Chem. Phys. 154, 064103 (2021). def2-universal-JFIT:
F. Weigend, Phys. Chem. Chem. Phys. 8, 1057-1065 (2006); def2-universal-JKFIT, def2-SV(P)-JKFIT: F. Weigend, J. Comput.
Chem. 29, 167-175 (2008). RIFIT: F. Weigend, M. Haeser, H. Patzelt, R. Ahlrichs, Chem. Phys. Lett. 294, 143-152 (1998);
C. Haettig, Phys. Chem. Chem. Phys. 7, 59-66 (2005); A. Hellweg, C. Haettig, S. Hoefener, W. Klopper, Theor. Chem. Acc.
117, 587-597 (2007); A. Hellweg, D. Rappoport, Phys. Chem. Chem. Phys. 17, 1010-1017 (2015). def2-ECP (per element, as
listed by BSE): Andrae et al. Theor. Chim. Acta 77, 123 (1990); Dolg et al. Theor. Chim. Acta 75, 173 (1989) and J. Chem.
Phys. 90, 1730 (1989); Leininger et al. Chem. Phys. Lett. 255, 274 (1996); Kaupp et al. J. Chem. Phys. 94, 1360 (1991);
Metz et al. J. Chem. Phys. 113, 2563 (2000) and Theor. Chem. Acc. 104, 22 (2000); Peterson et al. J. Chem. Phys. 119,
11113 (2003); Dolg et al. Theor. Chim. Acta 85, 441 (1993); Cao, Dolg J. Chem. Phys. 115, 7348 (2001). MINIX/HF-3c:
R. Sure, S. Grimme, J. Comput. Chem. 2013, 34, 1672-1685. Basis Set Exchange: B. P. Pritchard et al., J. Chem. Inf.
Model. 59, 4814-4820 (2019). The orbital and fitting references are cited at run time (`qm_citations.h`) whether the
file or the built-in copy is loaded; the ECP and lanthanide references are not (nothing uses those data yet).

Status: AI-generated, machine-tested (Claude Generated, Nov 2026).
