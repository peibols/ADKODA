# to compile: ./compiler.sh main
#!/bin/bash

g++ -g -std=c++11 -Wall \
	BerGEN.cc utils/fjcore.cc Util.cc Shower.cc Parton.cc InPartons.cc Evolution.cc Tests.cc \
	HepMC3_writer.cc ${HEPMC3_LIB}/libHepMC3-static.a -I${HEPMC3_INCLUDE} -L${HEPMC3_LIB} -lHepMC3 \
	$1.cc -o $1
