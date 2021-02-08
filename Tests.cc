#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Recluster.hh"
#include "fastjet/ClusterSequence.hh"

#include "Tests.h"
#include "BerGEN.h"
#include "Util.h"

using namespace fastjet;

namespace Adkoda {

void Test_FirstMass(std::vector<Parton> parton_list, std::ofstream &outfile, double mass_hist[2][20], double &njets) {

  std::vector<PseudoJet> particles;
  for (int i=0; i<parton_list.size(); i++) {
    if (parton_list[i].stat()<0) continue;
    fastjet::PseudoJet p;
    p.reset_momentum(parton_list[i].px(),parton_list[i].py(),parton_list[i].pz(),parton_list[i].e());
    p.set_user_index(i);
    particles.push_back(p);
  }

  double R=0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(particles, jet_def);
  double jet_ptmin=50.;
  std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(jet_ptmin));
  for (unsigned int ijet=0; ijet<jets.size(); ijet++) {
    //cout << " jetpt= " << jets[ijet].pt() << endl;
    double jetm=jets[ijet].m();
    //cout << " jet mass= " << jetm << endl;
    std::vector<PseudoJet> constituents = sorted_by_pt(jets[ijet].constituents());
    int ihard = constituents[0].user_index();
    int anc=ihard;
    int orig=-1000;
    double highest_splitmass=0.;
    while (true) {
      if (parton_list[anc].stat()==23) {
        cout << " hard never splitted! " << endl;
        break;
      }
      if (parton_list[anc].stat()==-23) {
        orig = anc;
        //cout << " found the orig = " << orig << endl;
        break;
      }
      int mom=parton_list[anc].mom1();
      int d1=parton_list[mom].d1();
      int d2=parton_list[mom].d2();
      if (d1!=0 && d2!=0 && d1!=d2) { 
        double split_mass = std::pow(parton_list[d1].e()+parton_list[d2].e(),2.)
				-std::pow(parton_list[d1].px()+parton_list[d2].px(),2.)
				-std::pow(parton_list[d1].py()+parton_list[d2].py(),2.)
				-std::pow(parton_list[d1].pz()+parton_list[d2].pz(),2.);
        //cout << " split mass= " << std::sqrt(split_mass) << endl;
        PseudoJet p1(parton_list[d1].px(),parton_list[d1].py(),parton_list[d1].pz(),parton_list[d1].e());
        PseudoJet p2(parton_list[d2].px(),parton_list[d2].py(),parton_list[d2].pz(),parton_list[d2].e());
        double deltaR=p1.delta_R(p2);
        double delR1=p1.delta_R(jets[ijet]);
        double delR2=p2.delta_R(jets[ijet]);
        //if (deltaR > R+0.1) cout << " split beyond R! " << endl;
        if (delR1>R || delR2>R) {
          //cout << " split beyond R! " << endl;
        }
        else highest_splitmass = std::sqrt(split_mass);
      }
      anc=parton_list[anc].mom1();
    }
    //cout << " HIGHEST= " << highest_splitmass << endl;
    int nm=int(jetm/2.);
    mass_hist[0][nm]+=1.;
    nm=int(highest_splitmass/2.);
    mass_hist[1][nm]+=1.;
    njets+=1.;
  }

}

//TEST MODULE: Save the initiator e-e+ kinematics and event_weights
void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  outfile << parton_list[4].e()  << " " << parton_list[4].px() << " " <<
             parton_list[4].py() << " " << parton_list[4].pz() << " " <<
             event_xsec << event_weight << std::endl;
}

//TEST MODULE: Save final particles (px, py, pz, E)
void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int multi = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    //if (parton_list[ip].stat()<0 && abs(parton_list[ip].id())!=11) continue; //MCnet python durham algorithm needs the incoming  particles.
    if (parton_list[ip].stat() < 0) continue;
     outfile << parton_list[ip].px() << " " << parton_list[ip].py() << " "
             << parton_list[ip].pz() << " " << parton_list[ip].e() << std::endl;
     multi++;
  }
  outfile << "#event xsec: " << event_xsec << ", weight: " << event_weight << ", multiplicity: " << multi << std::endl;
}


void Test_EnergyMomentumConservation(std::vector<Parton> parton_list) {
  double total_p[4] = {0., 0., 0., 0.};
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].stat()<0) continue; //MCnet python durham algorithm needs the incoming  particles.
    total_p[0] += parton_list[ip].e();
    total_p[1] += parton_list[ip].px();
    total_p[2] += parton_list[ip].py();
    total_p[3] += parton_list[ip].pz();
  }
  for (unsigned a = 0; a < 4; a++) std::cout << "Total momentum: p_" << a << " = " << total_p[a] << std::endl;
}

//TEST MODULE: Save the primary Lund plane using the history and FastJet definition of the variables
void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && ip != m1 && m2 == -1) { m2 = ip; break; }
  }
  if (parton_list[m1].stat() < 0) { //It decays
    int mom = m1;
    int d1  = parton_list[mom].d1();
    int d2  = parton_list[mom].d2();
    if (d1==d2) mom = m2; //Find the initial splitter
    //Soft Drop parameters
    double z_cut = 0.;
    double beta = 0.;
    double R = 3.;
    while (parton_list[mom].stat() < 0) {
      d1 = parton_list[mom].d1();
      d2 = parton_list[mom].d2();
      if (abs(parton_list[d1].stat()) == 51) { //A splitter
        if (parton_list[d1].pt() < parton_list[d2].pt()) { //Follow the primary
          int temp = d1;
          d1 = d2;
          d2 = temp;
        }
        double Delta    = Util::Delta(parton_list[d1].p(), parton_list[d2].p());
        double z        = parton_list[d2].pt() / (parton_list[d1].pt() + parton_list[d2].pt());
        if (z > z_cut * pow(Delta/R, beta) && Delta < R) { // Soft Drop condition.
          double kt     = parton_list[d2].pt() * Delta;
          double m2     = Util::m2(parton_list[d1].p(), parton_list[d2].p()); //Mass of the splitter
          double qt2    = m2 / z / (1.-z);
          double energy = parton_list[d1].e() + parton_list[d2].e(); //Energy of the splitter
          double tf     = 2. * energy / m2;
          outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " "
                  << Delta << " " << m2 << " " << tf << " " << "1." << std::endl;
        }
        mom = d1;
      }
      else if (parton_list[d1].stat()>0) break;
      else mom = d1;
    }
  }
  outfile << "#xsec: " << event_xsec << "weight: " << event_weight << std::endl;
}

/*
//TEST MODULE: Quench particle history by counting the splitting inside the medium
void Simple_Quenching(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int inSplit = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (parton_list[ip].x().t()*0.1973 < Rmed_in_fm) inSplit++;
  }
  Test_PrintFinalPartons(parton_list, event_xsec, std::pow(Quench, inSplit), outfile);
}
*/

} // end namespace Adkoda
