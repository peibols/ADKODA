#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "Tests.h"
#include "BerGEN.h"
#include "Util.h"

namespace Adkoda {

//TEST MODULE: Save the initiator e-e+ kinematics and event_weights
void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  double LinGeV = 1./0.193;
  double qhatinGeV = 1.;
  double nSplitIn = CountSplitsInMedium(parton_list, LinGeV, qhatinGeV);
  outfile << parton_list[4].e()  << " " << parton_list[4].px() << " " <<
             parton_list[4].py() << " " << parton_list[4].pz() << " " <<
             event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
}

//TEST MODULE: Save final particles (px, py, pz, E)
void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int multi = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].stat()<0 && abs(parton_list[ip].id())!=11) continue; //MCnet python durham algorithm needs the incoming  particles.
    //if (parton_list[ip].stat() < 0) continue;
     outfile << parton_list[ip].px() << " " << parton_list[ip].py() << " "
             << parton_list[ip].pz() << " " << parton_list[ip].e() << " "
             << parton_list[ip].id() << " "
             << parton_list[ip].col() << " " << parton_list[ip].acol() << std::endl;
     multi++;
  }
  outfile << "#event xsec: " << event_xsec << " weight: " << event_weight << " multiplicity: " << multi << std::endl;
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
void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight) {

  double LinGeV = 1./0.193;
  double qhatinGeV = 1.;
  double nSplitIn = CountSplitsInMedium(parton_list, LinGeV, qhatinGeV);
  outfile_weight << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;

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
    double R = 0.4;
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
          FourVector xSplit = parton_list[d1].x();
          FourVector xSpect = parton_list[parton_list[parton_list[mom].dippart()].d1()].x();
          outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " " <<
                  Delta << " " << m2 << " " << tf << " " <<
                  event_xsec << " " << event_weight << " " << nSplitIn << " " <<
                  (xSplit-xSpect).p3abs() << " " << xSplit.t() << std::endl;
        }
        mom = d1;
      }
      else if (parton_list[d1].stat()>0) break;
      else mom = d1;
    }
  }
}


//TEST MODULE: Save the primary Lund plane FastJet alhgorithm
void Test_PrintLundPlane_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight) {
  std::vector<fastjet::PseudoJet> event;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (parton_list[ip].stat() < 0) continue;
    fastjet::PseudoJet particle(parton_list[ip].p().x(), parton_list[ip].p().y(), parton_list[ip].p().z(), parton_list[ip].p().t());
    event.push_back(particle);
  }

  double R = 0.4, ptmin = 90., ptmax = 110.;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  fastjet::ClusterSequence cs(event, jet_def);
  std::vector<fastjet::PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

  double LinGeV = 1./0.193;
  double qhatinGeV = 1.;
  double nSplitIn = CountSplitsInMedium(parton_list, LinGeV, qhatinGeV);

  double z_cut = 0.0;
  double beta  = 0.0;
  fastjet::contrib::SoftDrop sd(beta, z_cut);

  fastjet::contrib::LundGenerator lund;
  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    fastjet::PseudoJet sd_jet = sd(jets[ijet]);
    if (sd_jet.pt() > ptmax) continue;
    std::vector<fastjet::contrib::LundDeclustering> declusts = lund(sd_jet);
    for (unsigned idecl = 0; idecl < declusts.size(); idecl++) {
      double z   = declusts[idecl].z();
      double m2  = declusts[idecl].pair().m2();
      double qt2 = m2/z/(1.-z);
      double en  = declusts[idecl].pair().E();
      double kt2 = declusts[idecl].kt()*declusts[idecl].kt();
      outfile << qt2 << " " << z << " " << en << " " << kt2 << " " <<
                 declusts[idecl].Delta() << " " << m2 << " " << 2.*en/m2 << " " <<
                 event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
    }
    if (declusts.size()>0) outfile_weight << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
  }
}

//TEST MODULE: Save the jet mass using FastJet alhgorithm
void Test_JetMass_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  std::vector<fastjet::PseudoJet> event;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (parton_list[ip].stat() < 0) continue;
    fastjet::PseudoJet particle(parton_list[ip].p().x(), parton_list[ip].p().y(), parton_list[ip].p().z(), parton_list[ip].p().t());
    event.push_back(particle);
  }

  double LinGeV = 1./0.193;
  double qhatinGeV = 1.;
  double nSplitIn = CountSplitsInMedium(parton_list, LinGeV, qhatinGeV);

  double R = 0.4, ptmin = 90.;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  fastjet::ClusterSequence cs(event, jet_def);
  std::vector<fastjet::PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

  double z_cut = 0.0;
  double beta  = 0.0;
  fastjet::contrib::SoftDrop sd(beta, z_cut);

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    fastjet::PseudoJet sd_jet = sd(jets[ijet]);
    //outfile << jets[ijet].m() << " " << jets[ijet].perp() << " " << event_xsec << " " << event_weight << std::endl;
    if (sd_jet.perp()<110) {
      outfile << sd_jet.m() << " " << sd_jet.perp() << " " << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
    }
  }
}


//Count the number of splittings resolved by the medium
double CountSplitsInMedium(std::vector<Parton> parton_list, double LinGeV, double qhatinGeV) {
  LinGeV=1./0.193;
  qhatinGeV=5.;
  int nSplitIn = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].stat() != -23 && parton_list[ip].stat() != -51 && parton_list[ip].stat() != -52) continue; //Skip all non shower particles
    if (parton_list[ip].d1() != parton_list[ip].d2() && parton_list[ip].dippart() != 0) { //Splitter
      int rec = parton_list[ip].dippart();
      FourVector xSplit = parton_list[parton_list[ip].d1()].x();
      FourVector xSpect = parton_list[parton_list[rec].d1()].x();
      double rSplit = xSplit.p3abs();
      double lDip = (xSplit-xSpect).p3abs();
      double lRes = std::sqrt(12./qhatinGeV/xSplit.t());
      //std::cout << rSplit << " ?< \t" << LinGeV << "\t " << lDip << " ?> \t" << lRes << std::endl;
      //if (rSplit < LinGeV && lDip > lRes) nSplitIn++;
      if (lDip > lRes) nSplitIn++;
    }
  }
  return nSplitIn;
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
