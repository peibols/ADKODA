#include <string>
#include <fstream>
#include "BerGEN.h"
#include "Tests.h"
#include "HepMC3/WriterAscii.h"

using namespace Adkoda;
using namespace Util;
using namespace std;

int main(int argc, char **argv) {

  //Read arguments and pass as input file
  assert(argc==2);
  string input_file = *(argv+1);

  remove("test/test_weights.out");
  remove("test/test_FinalPartons.out");
  remove("test/test_LundPlane.out");
  remove("test/test_kinematics.out");
  remove("test/test_veto.out");
  remove("test/test_HepMC3.out");

  HepMC3::WriterAscii outfile_test_HepMC3("test/test_HepMC3.hepmc");
  ofstream outfile_test_weights;
  if (!outfile_test_weights.is_open()) outfile_test_weights.open("test/test_weights.out", ios_base::app);
  ofstream outfile_test_FinalPartons;
  if (!outfile_test_FinalPartons.is_open()) outfile_test_FinalPartons.open("test/test_FinalPartons.out", ios_base::app);
  ofstream outfile_test_LundPlane;
  if (!outfile_test_LundPlane.is_open()) outfile_test_LundPlane.open("test/test_LundPlane.out", ios_base::app);

  cout << "#Start Program" << endl;

  BerGEN bergen(input_file);
  int nEv = bergen.number_events();
  bergen.init();
  for (int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 1 == 0) cout << "\n #Event: " <<  iEv << endl;
    bergen.next();
    bergen.print();

    //Test modules
    //double event_weight = bergen.get_event_weight();
    //double event_xsec   = bergen.get_event_xsec();
    std::vector<Parton> parton_list = bergen.get_parton_list();
    //Test_Weights(parton_list, event_xsec, event_weight, outfile_test_weights);
    //Test_PrintFinalPartons(parton_list, event_xsec, event_weight, outfile_test_FinalPartons);
    //Test_EnergyMomentumConservation(parton_list);
    //Test_PrintLundPlane_history(parton_list, event_xsec, event_weight, outfile_test_LundPlane);

    //HepMC3 ASCII writer
    //write_HepMC3_event(parton_list, event_xsec, event_weight, outfile_test_HepMC3, iEv);

    // Freezing time hist
    for (unsigned int i=0; i<parton_list.size(); i++) {
      double time = parton_list[i].x().t();
      double en = parton_list[i].e();
      //std::cout << " en= " << en << " time= " << time << std::endl;
    }

  }
  cout << "#End Program" << endl;

  //Closing outfiles
  outfile_test_weights.close();
  outfile_test_FinalPartons.close();
  outfile_test_LundPlane.close();

  return 0;
} // End program
