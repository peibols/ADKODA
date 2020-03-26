//----------------------------------------------------------------------
// g++ -O3 -Wall -g -I/Users/ata053/Physics/fastjet-install/include   -c -o test_lundplane.o test_lundplane.cc
// g++ -O3 -Wall -g -I/Users/ata053/Physics/fastjet-install/include   -c -o LundGenerator.o LundGenerator.cc
// g++ -O3 -Wall -g -I/Users/ata053/Physics/fastjet-install/include -o test_lundplane test_lundplane.o -L. -lLundPlane -lm -Wl,-rpath,/Users/ata053/Physics/fastjet-install/lib -lm -L/Users/ata053/Physics/fastjet-install/lib -lfastjettools -lfastjet


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/contrib/SoftDrop.hh"
using namespace std;
using namespace fastjet;

void read_event(vector<PseudoJet> &event);

int main(){

  int nEv = 5000;
  int nJet = 0;
  remove("test_JetRates.dat");
  ofstream out_LP;
  out_LP.open("JetRates.dat");
  for (int iEv = 0; iEv < nEv; iEv++){

    vector<PseudoJet> event;
    read_event(event);

    double R = 3., ptmin = 40., ptmax = 60.;
    JetDefinition jet_def(ee_kt_algorithm, R); //Durham algorithm
    ClusterSequence cs(event, jet_def);
    vector<PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

    //double z_cut = 0.0;
    //double beta  = 0.0;
    //contrib::SoftDrop sd(beta, z_cut);

    //----------------------------------------------------------
    // create an instance of LundGenerator, with default options
    contrib::LundGenerator lund;
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
      //PseudoJet sd_jet = sd(jets[ijet]);
      //if (sd_jet.pt() > ptmax || abs(sd_jet.rap()) > 3.) continue;
      //if (jets[ijet].pt() > ptmax || abs(jets[ijet].rap()) > 3.) continue;
      //vector<contrib::LundDeclustering> declusts = lund(sd_jet);
      vector<contrib::LundDeclustering> declusts = lund(jets[ijet]);
      for (unsigned idecl = 0; idecl < declusts.size(); idecl++) {
        out_LP << (declusts[idecl].m()*declusts[idecl].m())/declusts[idecl].z()/(1.-declusts[idecl].z()) << " " <<
                   declusts[idecl].z() << " " << declusts[idecl].pair().E() << " " <<
                   declusts[idecl].kt()*declusts[idecl].kt() << " " << declusts[idecl].Delta() << " " <<
                   declusts[idecl].m()*declusts[idecl].m() << " " <<
                   2.*declusts[idecl].pair().E()/declusts[idecl].m()/declusts[idecl].m() << " " << 1. << endl;
      }
      out_LP << "#EndJet pt: " << jets[ijet].pt() << endl;
      nJet++;
    }
    out_LP << "#EndEvent: " << iEv << endl;

  }
  out_LP << "# Total number of analyzed events: " << nEv << endl;
  out_LP << "# Total number of analyzed jets: " << nJet << endl;
  out_LP.close();

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){
  string line;
  string scheck;
  while (getline(cin, line)) {
    //istringstream linestream(line);
    // take substrings to avoid problems when there is extra "pollution"
    // characters (e.g. line-feed).
    istringstream ifi(line);
    ifi >> scheck;
    if (scheck=="#") return;
    double px,py,pz,E;
    px=strtod(scheck.c_str(), NULL);
    //cout << " px= " << px << endl;
    ifi >> py >> pz >> E;
    //if (line.substr(0,1) == "#") {return;}
    //if (line.substr(0,1) == "#") {continue;}

    //linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}
