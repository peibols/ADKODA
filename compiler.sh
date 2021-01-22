# to compile: ./compiler.sh main
#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/ata053/Physics/HepMC3/hepmc3-install/lib:/Users/ata053/Physics/HepMC2/hepmc-install/lib:/Users/ata053/Physics/pythia8235/lib



HEPMC_INCLUDE=/Users/ata053/Physics/HepMC2/hepmc-install/include
HEPMC_LIB=/Users/ata053/Physics/HepMC2/hepmc-install/lib

HEPMC3_INCLUDE=/Users/ata053/Physics/HepMC3/hepmc3-install/include
HEPMC3_LIB=/Users/ata053/Physics/HepMC3/hepmc3-install/lib

FASTJET=/Users/ata053/Physics/fastjet-install/bin

g++ -g -std=c++11 -Wall \
	BerGEN.cc utils/fjcore.cc Util.cc Shower.cc Parton.cc InPartons.cc Evolution.cc Tests.cc \
	HepMC_writer.cc ${HEPMC_LIB}/libHepMC.a -I${HEPMC_INCLUDE} -L${HEPMC_LIB} -lHepMC \
	HepMC3_writer.cc ${HEPMC3_LIB}/libHepMC3-static.a -I${HEPMC3_INCLUDE} -L${HEPMC3_LIB} -lHepMC3 \
	`${FASTJET}/fastjet-config --cxxflags --libs --plugins` -lRecursiveTools -lLundPlane \
	$1.cc -o $1
