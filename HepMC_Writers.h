#include "Parton.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC/IO_GenEvent.h"

namespace Adkoda {

void write_HepMC_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC::IO_GenEvent &outfile, int nEv);

void write_HepMC3_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC3::WriterAscii &outfile, int iEv);

} // end namespace Adkoda

