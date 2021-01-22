#include "Shower.h"
#include "BerGEN.h"
#include "InPartons.h"
#include "data.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC/IO_GenEvent.h"

namespace Adkoda {

void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

void Test_EnergyMomentumConservation(std::vector<Parton> parton_list);

void Test_FirstSplit(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight);

void Test_PrintLundPlane_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight);

void Test_JetMass_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

void Test_IteratedSoftDrop_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

std::vector<int> ListFinalDescendants(std::vector<Parton> parton_list, int ip);

std::vector<int> ListSplits(std::vector<Parton> parton_list, int ip);

std::vector<int> CountRealSplits(std::vector<std::vector<int>> split_arr);

int CountSplitsInMedium(std::vector<Parton> parton_list, std::vector<int> splits);

void write_HepMC_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC::IO_GenEvent &outfile, int nEv);

void write_HepMC3_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC3::WriterAscii &outfile, int iEv);

} // end namespace Adkoda
