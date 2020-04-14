# ADKODA

## HepMC3

HepMC3 needs to be installed, and the following paths need to be set:

`HEPMC3_BIN`, `HEPMC3_INCLUDE` , `HEPMC3_LIB`.

One also needs to add `HEPMC3_LIB` to `LD_LIBRARY_PATH`.

Following Appendix A in the [HepMC3 manual](https://arxiv.org/pdf/1912.08005.pdf).

Status codes for particles:

| Status        | Definiton     		|
|:-------------:|:-----------------------------:|
| 1             | Decayed physical particle 	|
| 2      	| Undecayed physical particle   |
| 52 		| Decays into carbon copy      	|

Status codes for vertices:

| Status        | Definiton     		|
|:-------------:|:-----------------------------:|
| 1             | 1->2 splitting	 	|
| 2      	| 1->1 splitting (carbon copy)  |

