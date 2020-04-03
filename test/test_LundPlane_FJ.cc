// g++ -O3 -Wall test_LundPlane_FJ.cc -o test_LundPlane_FJ `/Users/ata053/Physics/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lRecursiveTools -lLundPlane
// ./test_LundPlane_FJ input_file

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/contrib/SoftDrop.hh"
using namespace std;
using namespace fastjet;

void read_event(vector<PseudoJet> &event, double &event_weight, int &multiplicity, std::ifstream &infile);

int main(int argc, char **argv){

  //Read input
  assert(argc==2);
  string input_file = *(argv+1);
  std::ifstream infile(input_file);
  if(!infile.is_open()) { perror("Error: opening file."); exit(EXIT_FAILURE); }
  //Prepare output
  remove("test_LundPlane_FJ.dat");
  ofstream outfile("test_LundPlane_FJ.dat");

  int nEv  = 10000;
  int nJet = 0;
  for (int iEv = 0; iEv < nEv; iEv++){
    vector<PseudoJet> event;
    double event_weight = 0.;
    int multiplicity = 0;
    read_event(event, event_weight, multiplicity, infile);

    double R = 3., ptmin = 40., ptmax = 60.;
    JetDefinition jet_def(antikt_algorithm, R);
    ClusterSequence cs(event, jet_def);
    vector<PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

    double z_cut = 0.0;
    double beta  = 0.0;
    contrib::SoftDrop sd(beta, z_cut);

    contrib::LundGenerator lund;
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      PseudoJet sd_jet = sd(jets[ijet]);
      if (sd_jet.pt() > ptmax || abs(sd_jet.rap()) > 5.) continue;
      vector<contrib::LundDeclustering> declusts = lund(sd_jet);
      for (unsigned idecl = 0; idecl < declusts.size(); idecl++) {
        outfile << declusts[idecl].pair().m2()/declusts[idecl].z()/(1.-declusts[idecl].z()) << " " <<
                   declusts[idecl].z() << " " << declusts[idecl].pair().E() << " " <<
                   declusts[idecl].kt()*declusts[idecl].kt() << " " << declusts[idecl].Delta() << " " <<
                   declusts[idecl].pair().m2() << " " <<
                   2.*declusts[idecl].pair().E()/declusts[idecl].pair().m2() << " " << 1. << endl;
      }
      outfile << "#EndJet pt: " << jets[ijet].pt() << endl;
      nJet++;
    }
    outfile << "#EndEvent: " << iEv << " weight: " << event_weight << endl;

  }
  outfile << "#Total_events: " << nEv << " Total_jets: " << nJet << endl;
  infile.close();
  outfile.close();

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event, double &event_weight, int &multiplicity, std::ifstream &infile){
  string line, scheck, trash1, trash2;
  while (std::getline(infile, line)) {
    std::istringstream ifi(line);
    ifi >> scheck;
    if (scheck == "#event") {
      ifi >> trash1 >> event_weight >> trash2 >> multiplicity;
      return;
    }
    double px, py, pz, E;
    px = strtod(scheck.c_str(), NULL);
    ifi >> py >> pz >> E;
    PseudoJet particle(px, py, pz, E);
    event.push_back(particle);
  }
}
