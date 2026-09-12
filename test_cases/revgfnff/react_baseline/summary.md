# React baseline with corrected forces (post Known Issue #28)

AI-generated (scripts/react_baseline.py), machine-measured. Same seed (42) for every run; MD with CSVR (10 fs), dt 0.5 fs.

| run | formed | broken | rebuilds | exchange | dE_jump n / median / min / max [kJ/mol] | NaN | max r [A] | final fragments | wall [s] |
|---|---:|---:|---:|---:|---|---|---:|---|---:|
| R1_h4_T1000 | 0 | 0 | 0 | 0 | 0 | no | 4.83 | H + H + H + H | 0.2 |
| R1_h4_T2000 | 0 | 0 | 0 | 0 | 0 | no | 5.78 | H + H + H + H | 0.2 |
| R1_h4_T3000 | 0 | 0 | 0 | 0 | 0 | no | 8.64 | H + H + H + H | 0.3 |
| R1_h4_T4000 | 0 | 0 | 0 | 0 | 0 | no | 7.21 | H + H + H + H | 0.3 |
| R1_h4_T5000 | 1 | 1 | 2 | 0 | 2 / -222.4 / -451.4 / +6.6 | no | 9.71 | H + H + H + H | 0.3 |
| R1_h4_T6000 | 3 | 3 | 6 | 0 | 6 / -209.9 / -548.0 / +33.5 | no | 8.78 | H + H + H + H | 0.3 |
| R2_2h2_break_early_1.45 | 0 | 2 | 2 | 0 | 2 / +469.6 / +463.8 / +475.3 | no | 8.87 | H + H + H + H | 0.2 |
| R2_2h2_break_default_2.6 | 0 | 0 | 0 | 0 | 0 | no | 1.69 | H + H + H2 | 0.2 |
| R3_n2_3h2_filters_on | 34 | 33 | 67 | 0 | 67 / -74.0 / -575.6 / +629.2 | no | 5.58 | H + H2 + N2H3 | 2.1 |
| R3_n2_3h2_filters_off | 147 | 150 | 288 | 0 | 288 / +1.9 / -1389.6 / +1124.2 | no | 11.93 | H + H + H + H + H + H + N2 | 2.1 |
| R4_n4h4_slack_on | 20 | 19 | 37 | 0 | 37 / +24.1 / -598.0 / +610.7 | no | 5.14 | H + N2H + N2H2 | 1.7 |
| R4_n4h4_slack_off | 189 | 190 | 366 | 0 | 366 / -3.2 / -891.6 / +1261.4 | no | 6.96 | H + H + H + N2 + N2H | 1.7 |
| R5_wall5ps_harmonic_298.15 | 2 | 3 | 5 | 0 | 5 / +33.4 / -464.0 / +447.8 | no | 11.63 | H + H + H + H + H + H2 + H2 + H2 + N2 + N2H | 0.8 |
| R5_wall5ps_harmonic_10000 | 9 | 8 | 17 | 0 | 17 / +14.3 / -461.2 / +517.5 | no | 4.74 | H + H + H + H2 + H2 + N2H2 + N2H3 | 1.0 |
| R5_wall5ps_logfermi_298.15 | 6 | 5 | 11 | 0 | 11 / +27.1 / -541.3 / +556.8 | no | 5.61 | H + H2 + H2 + H2 + H2 + N2H + N2H2 | 0.9 |
| R5_wall5ps_logfermi_10000 | 9 | 6 | 15 | 0 | 15 / +41.9 / -741.2 / +763.7 | no | 3.44 | H + H + H2 + H2 + N2H2 + N2H4 | 1.0 |
| R5_wall20ps_harmonic | 42 | 41 | 83 | 0 | 83 / -39.1 / -585.6 / +597.1 | no | 10.33 | H + H + H + H2 + N2H4 + NH + NH2 | 3.5 |
| R5_wall20ps_logfermi | 28 | 26 | 52 | 0 | 52 / -53.1 / -541.3 / +584.5 | no | 4.49 | H + H + H2 + N2H4 + NH + NH3 | 3.6 |
| R5_wall20ps_pbc | 57 | 53 | 109 | 0 | 109 / -52.5 / -571.9 / +598.7 | no | 5.41 | NH3 + NH3 + NH3 + NH3 | 4.0 |
| R6_2h2_2500K | 1 | 2 | 3 | 0 | 3 / +12.3 / -455.4 / +20.3 | no | 5.75 | H + H + H + H | 0.2 |
