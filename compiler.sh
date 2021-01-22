# to compile: ./compiler.sh main
#!/bin/bash

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/ata053/Physics/hepmc3-install/lib:/Users/ata053/Physics/pythia8235/lib

#HEPMC3_INCLUDE=/Users/ata053/Physics/hepmc3-install/include
#HEPMC3_LIB=/Users/ata053/Physics/hepmc3-install/lib

FASTJET3_BIN=/Users/peibols/Software/fastjet-install/bin
FASTJET3_INCLUDE=/Users/peibols/Software/fastjet-install/include
FASTJET3_LIB=/Users/peibols/Software/fastjet-install/lib

g++ -g -std=c++11 -Wall \
	BerGEN.cc utils/fjcore.cc Util.cc Shower.cc Parton.cc InPartons.cc Evolution.cc Tests.cc \
	Cascade.cc MedEvolution.cc Broadening.cc \
	HepMC3_writer.cc ${HEPMC3_LIB}/libHepMC3-static.a -I${HEPMC3_INCLUDE} -L${HEPMC3_LIB} -lHepMC3 \
	$1.cc -o $1 \
        -I${FASTJET3_INCLUDE} -L${FASTJET3_LIB} -Wl,-rpath,${FASTJET3_LIB} -lfastjet `${FASTJET3_BIN}/fastjet-config --cxxflags --libs --plugins` -lRecursiveTools -lNsubjettiness
       
