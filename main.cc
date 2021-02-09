#include <string>
#include <fstream>
#include "HepMC_Writers.h"

#include "BerGEN.h"
#include "Tests.h"

using namespace Adkoda;
using namespace Util;
using namespace std;

int main(int argc, char **argv) {

  //Read arguments and pass as input file
  assert(argc==2);
  string input_file = *(argv+1);

  //HepMC 2&3
  HepMC::IO_GenEvent outfile_test_HepMC("test/test_HepMC.hepmc", std::ios::out);
  HepMC3::WriterAscii outfile_test_HepMC3("test/test_HepMC3.hepmc");

  Tests tests;

  cout << "#Start Program" << endl;

  BerGEN bergen(input_file);
  int nEv = bergen.number_events();
  bergen.init();

  for (int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 1000 == 0) cout << "#Event: " <<  iEv << endl;
    bergen.next();
    //bergen.print();

    double event_weight = bergen.get_event_weight();
    double event_xsec   = bergen.get_event_xsec();
    std::vector<Parton> parton_list = bergen.get_parton_list();
  
    write_HepMC_event(parton_list, event_xsec, event_weight, outfile_test_HepMC, iEv);
    write_HepMC3_event(parton_list, event_xsec, event_weight, outfile_test_HepMC3, iEv);

    tests.run(parton_list, event_weight, event_xsec, iEv);

  }

  cout << "#End Program" << endl;

  tests.close(nEv);

  return 0;

} // End program
