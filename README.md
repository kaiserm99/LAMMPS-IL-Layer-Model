# Zu simulieren:


## ILMs


| Simulation Szenario | 3:1:2 | 5:1:4 | 7:1:6 |
| ------------------- | ----- | ----- | ----- |
| Neutral             | 1_1   | 2_1   | 3_1   |
| Oberfläche Negativ  | 1_2   | 2_2   | 3_2   |
| Oberfläche Positiv  | 1_3   | 2_3   | 3_3   |


## ILs

| Simulation Szenario | P66614:BTA | P66614:BMB | P66614:Doc | P66614:Cl                                       |
| ------------------- | ---------- | ---------- | ---------- | ----------------------------------------------- |
| Neutral             | 1-1        | 1-2        | 1-3        | 1-4                                             |
| Oberfläche Negativ  | 2-1        | 2-2        | 2-3        | 2-4 (nicht möglich weil rutscht komplett durch) |
| Oberfläche Positiv  | 3-1 (1x)   | 3-2 (2x)   | 3-3   (3x) | Das geht nicht                                  |


## Größen

- P66614: 1.357 nm --> 0.0025
- BMB: 0.707 nm --> 0.001302
- Doc: 0.943 nm  --> 0.001736
- BTA: 0.364 nm --> 0.0006707
- Cl: 0.167 nm --> 0.000308

mpirun -np 14 lmp_auto -in in.simulation