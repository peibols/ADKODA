<<<<<<< HEAD
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

## Pythia8

Pythia8 needs to be installed, and the following paths need to be set:

`PYTHIA8_INCLUDE` , `PYTHIA8_LIB`.

One also needs to add `PYTHIA8_LIB` to `LD_LIBRARY_PATH`.

Additionally, one needs to set the path `PYTHIA8DATA` to point to the `xmldoc` folder in the Pythia8 installation. For example, in my profile I have:

`export PYTHIA8DATA="/Users/peibols/Software/pythia8240/share/Pythia8/xmldoc"`.

### Identified unsolved issues

+ When Pythia8 is used, we also keep the beam remnants, necessary to keep track of momentum conservation and color flow. They have `status=63`. They are properly ignored in `Evolution.cc`, but they seem to cause, sometimes, trouble in test function `Test_PrintLundPlane`.

+ HepMC3 writer is not properly adapted to account for including beam remnants

+ Still need to set the daughters of incoming hard partons with `status=-21`. This is mostly necessary in order for them to be printed in the HepMC3 format.
=======
# BerGen

To compile and run:
$./compiler.sh main
$./main bergen_input
>>>>>>> medium_dani
