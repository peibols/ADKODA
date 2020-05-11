#include <string>
#include <fstream>
#include "BerGEN.h"
#include "Tests.h"
#include "HepMC/IO_GenEvent.h"
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
  remove("test/test_LundPlane_history.out");
  remove("test/test_LundPlane_history_weight.out");
  remove("test/test_LundPlane_FJ.out");
  remove("test/test_LundPlane_FJ_weight.out");
  remove("test/test_JetMass_FJ.out");
  remove("test/test_kinematics.out");
  remove("test/test_veto.out");
  remove("test/test_HepMC.hepmc");
  remove("test/test_HepMC3.hepmc");

  HepMC::IO_GenEvent outfile_test_HepMC("test/test_HepMC.hepmc", std::ios::out);
  HepMC3::WriterAscii outfile_test_HepMC3("test/test_HepMC3.hepmc");
  ofstream outfile_test_weights;
  if (!outfile_test_weights.is_open()) outfile_test_weights.open("test/test_weights.out", ios_base::app);
  ofstream outfile_test_FinalPartons;
  if (!outfile_test_FinalPartons.is_open()) outfile_test_FinalPartons.open("test/test_FinalPartons.out", ios_base::app);
  ofstream outfile_test_LundPlane_history;
  if (!outfile_test_LundPlane_history.is_open()) outfile_test_LundPlane_history.open("test/test_LundPlane_history.out", ios_base::app);
  ofstream outfile_test_LundPlane_history_weight;
  if (!outfile_test_LundPlane_history_weight.is_open()) outfile_test_LundPlane_history_weight.open("test/test_LundPlane_history_weight.out", ios_base::app);
  ofstream outfile_test_LundPlane_FJ;
  if (!outfile_test_LundPlane_FJ.is_open()) outfile_test_LundPlane_FJ.open("test/test_LundPlane_FJ.out", ios_base::app);
  ofstream outfile_test_LundPlane_FJ_weight;
  if (!outfile_test_LundPlane_FJ_weight.is_open()) outfile_test_LundPlane_FJ_weight.open("test/test_LundPlane_FJ_weight.out", ios_base::app);
  ofstream outfile_test_JetMass_FJ;
  if (!outfile_test_JetMass_FJ.is_open()) outfile_test_JetMass_FJ.open("test/test_JetMass_FJ.out", ios_base::app);

  cout << "#Start Program" << endl;

  BerGEN bergen(input_file);
  int nEv = bergen.number_events();
  bergen.init();

  for (int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 1000 == 0) cout << "#Event: " <<  iEv << endl;
    bergen.next();
    //bergen.print();

    //Test modules
    double event_weight = bergen.get_event_weight();
    double event_xsec   = bergen.get_event_xsec();
    std::vector<Parton> parton_list = bergen.get_parton_list();
    Test_Weights(parton_list, event_xsec, event_weight, outfile_test_weights);
    //Test_PrintFinalPartons(parton_list, event_xsec, event_weight, outfile_test_FinalPartons);
    //Test_EnergyMomentumConservation(parton_list);
    Test_PrintLundPlane_history(parton_list, event_xsec, event_weight, outfile_test_LundPlane_history, outfile_test_LundPlane_history_weight);
    Test_PrintLundPlane_FJ(parton_list, event_xsec, event_weight, outfile_test_LundPlane_FJ, outfile_test_LundPlane_FJ_weight);
    //Test_JetMass_FJ(parton_list, event_xsec, event_weight, outfile_test_JetMass_FJ);

    //write_HepMC_event(parton_list, event_xsec, event_weight, outfile_test_HepMC, iEv);
    //write_HepMC3_event(parton_list, event_xsec, event_weight, outfile_test_HepMC3, iEv);

  }
  cout << "#End Program" << endl;

  //Closing outfiles
  outfile_test_weights.close();
  outfile_test_FinalPartons.close();
  outfile_test_LundPlane_history.close();
  outfile_test_LundPlane_history_weight.close();
  outfile_test_LundPlane_FJ.close();
  outfile_test_LundPlane_FJ_weight.close();
  outfile_test_JetMass_FJ.close();

  return 0;
} // End program
