#!/bin/bash

g++ -g -std=c++11 -Wall -Wno-deprecated \
	BerGEN.cc utils/fjcore.cc Util.cc Shower.cc Parton.cc InPartons.cc Evolution.cc Tests.cc \
	Cascade.cc MedEvolution.cc Broadening.cc \
        HepMC_writer.cc ${HEPMC_LIB}/libHepMC.so.4 -I${HEPMC_INCLUDE} -L${HEPMC_LIB} -lHepMC \
	HepMC3_writer.cc ${HEPMC3_LIB}/libHepMC3.so.3 -I${HEPMC3_INCLUDE} -L${HEPMC3_LIB} -lHepMC3 \
	InPythia.cc ${PYTHIA_LIB}/libpythia8.a -I${PYTHIA_INCLUDE} -L${PYTHIA_LIB} -Wl,-rpath,${PYTHIA_LIB} -lpythia8 -ldl \
        -I${FASTJET3_INCLUDE} -L${FASTJET3_LIB} -Wl,-rpath,${FASTJET3_LIB} -lfastjet `${FASTJET}/fastjet-config --cxxflags --libs --plugins` -lRecursiveTools -lNsubjettiness -lLundPlane \
	$1.cc -o $1
