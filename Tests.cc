#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <algorithm>

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/LundGenerator.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/IteratedSoftDrop.hh"

#include "Tests.h"
#include "BerGEN.h"
#include "Util.h"

namespace Adkoda {

//TEST MODULE: Save the initiator e-e+ kinematics and event_weights
void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  outfile << parton_list[4].e()  << " " << parton_list[4].px() << " " <<
             parton_list[4].py() << " " << parton_list[4].pz() << " " <<
             event_xsec << " " << event_weight << std::endl;
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
             << parton_list[ip].col() << " " << parton_list[ip].acol() << " "
             << parton_list[ip].x().x() << " " << parton_list[ip].x().y() << " "
             << parton_list[ip].x().z() << " " << parton_list[ip].x().t() << std::endl;
     multi++;
  }
  outfile << "#event xsec: " << event_xsec << " weight: " << event_weight << " multiplicity: " << multi << std::endl;
}

//TEST MODULE: Sum up the energy momentum of the final-state particles
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

//TEST MODULE: Save the first splitting using the history
void Test_FirstSplit(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  //Find primaries
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hard-scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && ip != m1 && m2 == -1) { m2 = ip; break; }
  }
  if (parton_list[m1].stat() < 0) { //It decays
    int d1  = parton_list[m1].d1();
    int d2  = parton_list[m1].d2();
    if (d1==d2) m1 = m2; //Find the initial splitter
  } else return;

  //List all final particles of the splitter
  std::vector<int> FinalDesc_arr = ListFinalDescendants(parton_list, m1);
  std::vector<std::vector<int>> ConstSplit_arr;
  for (unsigned int j=0; j<FinalDesc_arr.size(); j++){
    std::vector<int> ConstSplit = ListSplits(parton_list, FinalDesc_arr[j]);
    if (ConstSplit.size()>0) ConstSplit_arr.push_back(ConstSplit);
  }
  std::vector<int> ConstSplitReal = CountRealSplits(ConstSplit_arr);
  int nSplitIn = CountSplitsInMedium(parton_list, ConstSplitReal);

  int d1 = parton_list[m1].d1();
  int d2 = parton_list[m1].d2();
  double Delta  = Util::Delta(parton_list[d1].p(), parton_list[d2].p());
  double kt     = parton_list[d2].p().pt() * Delta;
  double pt_ev  = parton_list[d1].scale();
  double mass2     = Util::m2(parton_list[d1].p(), parton_list[d2].p()); //Mass of the splitter
  //double z      = parton_list[d2].p().pt() / (parton_list[d1].p().pt() + parton_list[d2].p().pt());
  double z      = 0.5*(1.-std::sqrt(1.-4.*pt_ev*pt_ev/mass2));
  double qt2    = mass2 / z / (1.-z);
  double energy = parton_list[d1].e() + parton_list[d2].e(); //Energy of the splitter after the shift.
  FourVector xSplit_beg = parton_list[m1].x();
  FourVector xSplit_end = parton_list[d1].x();
  FourVector xSpect_end = parton_list[parton_list[parton_list[m1].dippart()].d1()].x();
  double tf = xSplit_end.t()-xSplit_beg.t();
  outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " " <<
             Delta << " " << mass2 << " " << tf << " " <<
             event_xsec << " " << event_weight << " " << nSplitIn << " " <<
             pt_ev << std::endl;
}

//TEST MODULE: Save the primary Lund plane using the history and FastJet definition of the variables
void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight) {
  //Find primary splitter
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hard-scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && ip != m1 && m2 == -1) { m2 = ip; break; }
  }
  if (parton_list[m1].stat() < 0) { //It decays
    int d1  = parton_list[m1].d1();
    int d2  = parton_list[m1].d2();
    if (d1==d2) m1 = m2; //Find the initial splitter
  }
  else return;

  //Soft Drop parameters
  double z_cut = 0.;
  double beta  = 0.;
  double R     = 10.4;

  int nSplitIn = 0;
  int nSplit = 0;

  int mom = m1;
  while (parton_list[mom].stat() < 0) {
    int d1 = parton_list[mom].d1();
    int d2 = parton_list[mom].d2();
    if (abs(parton_list[d1].stat()) == 51) { //A splitter
      if (parton_list[d1].pt() < parton_list[d2].pt()) { //Follow the primary
        int temp = d1;
        d1 = d2;
        d2 = temp;
      }
      double Delta    = Util::Delta(parton_list[d1].p(), parton_list[d2].p());
      double z        = parton_list[d2].p().pt() / (parton_list[d1].p().pt() + parton_list[d2].p().pt());
      if (z > z_cut * pow(Delta/R, beta) && Delta < R) { // Soft Drop condition.
        double kt     = parton_list[d2].p().pt() * Delta;
        double pt_ev  = parton_list[d1].scale();
        double m2     = Util::m2(parton_list[d1].p(), parton_list[d2].p()); //Mass of the splitter
        double qt2    = m2 / z / (1.-z);
        double energy = parton_list[d1].e() + parton_list[d2].e(); //Energy of the splitter after the shift.
        FourVector xSplit_beg = parton_list[mom].x();
        FourVector xSplit_end = parton_list[d1].x();
        FourVector xSpect_end = parton_list[parton_list[parton_list[mom].dippart()].d1()].x();
        double tf = xSplit_end.t()-xSplit_beg.t();
        outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " " <<
                   Delta << " " << m2 << " " << tf << " " <<
                   event_xsec << " " << event_weight << " " << nSplitIn << " " <<
                   pt_ev << std::endl;
        nSplit++;
        if (nSplit==2) break;
      }
      mom = d1;
    }
    else if (parton_list[d1].stat()>0) break;
    else mom = d1;
  }
  outfile << "#" << std::endl;
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

  int nSplitIn=0;

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
    particle.set_user_index(ip); //Label particle
    event.push_back(particle);
  }

  double R = 0.4, ptmin = 90.;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  //fastjet::JetDefinition jet_def(fastjet::kt_algorithm, R);
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, R);
  fastjet::ClusterSequence cs(event, jet_def);
  std::vector<fastjet::PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

  double z_cut = 0.0;
  double beta  = 0.0;
  fastjet::contrib::SoftDrop sd(beta, z_cut);

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    if (jets[ijet].perp() < 100.) {
      std::vector<fastjet::PseudoJet> constituents = jets[ijet].constituents();
      std::vector<std::vector<int>> ConstSplit_arr;
      for (unsigned int j=0; j<constituents.size(); j++){
        std::vector<int> ConstSplit = ListSplits(parton_list, constituents[j].user_index());
        if (ConstSplit.size()>0) ConstSplit_arr.push_back(ConstSplit);
      }
      std::vector<int> ConstSplitReal = CountRealSplits(ConstSplit_arr);
      int nSplitIn = CountSplitsInMedium(parton_list, ConstSplitReal);
      fastjet::PseudoJet sd_jet = sd(jets[ijet]);
      outfile << sd_jet.m() << " " << sd_jet.perp() << " " << sd_jet.e() << " " << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
    }
  }
}

//TEST MODULE: Save the SoftDrop using FastJet alhgorithm
void Test_IteratedSoftDrop_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  std::vector<fastjet::PseudoJet> event;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (parton_list[ip].stat() < 0) continue;
    fastjet::PseudoJet particle(parton_list[ip].p().x(), parton_list[ip].p().y(), parton_list[ip].p().z(), parton_list[ip].p().t());
    particle.set_user_index(ip); //Label particle
    event.push_back(particle);
  }

  double R = 0.4, ptmin = 90., ptmax = 100.;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  fastjet::ClusterSequence cs(event, jet_def);
  std::vector<fastjet::PseudoJet> jets = sorted_by_rapidity(cs.inclusive_jets(ptmin));

  double z_cut = 0.1;
  double beta  = 0.0;
  double theta_cut = 0.0;
  fastjet::contrib::IteratedSoftDrop isd(beta, z_cut, theta_cut, R);
  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    if (jets[ijet].perp() < ptmax) {
      std::vector<fastjet::PseudoJet> constituents = jets[ijet].constituents();
      std::vector<std::vector<int>> ConstSplit_arr;
      for (unsigned int j=0; j<constituents.size(); j++){
        std::vector<int> ListSplit = ListSplits(parton_list, constituents[j].user_index());
        if (ListSplit.size()>0) ConstSplit_arr.push_back(ListSplit);
      }
      std::vector<int> ConstSplitReal = CountRealSplits(ConstSplit_arr);
      int nSplitIn = CountSplitsInMedium(parton_list, ConstSplitReal);
      fastjet::contrib::IteratedSoftDropInfo isd_info = isd(jets[ijet]);
      if (isd_info.multiplicity()>0) {
        outfile << isd_info[0].first << " " << isd_info[0].second << " " << isd_info.multiplicity() << " " << jets[ijet].m() << " " << jets[ijet].perp() << " " << jets[ijet].e() << " " << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
      } else outfile << 0.0 << " " << 0.0 << " " << isd_info.multiplicity() << " " << jets[ijet].m() << " " << jets[ijet].perp() << " " << jets[ijet].e() << " " << event_xsec << " " << event_weight << " " << nSplitIn << std::endl;
    }
  }
}

//List the all final state descendants for a given parton
std::vector<int> ListFinalDescendants(std::vector<Parton> parton_list, int ip) {
  std::vector<int> FinalDescendants;
  for (unsigned int i = 0; i < parton_list.size(); i++) {
    if (parton_list[i].stat() < 0) continue; //Skip non-finals
    std::vector<int> iSplits = ListSplits(parton_list, i); //List all ancestors of a final
    if(std::find(iSplits.begin(), iSplits.end(), ip) != iSplits.end()) {
      FinalDescendants.push_back(i);
    }
  }
  return FinalDescendants;
}

//List the splittings for which both daughters are exist
std::vector<int> CountRealSplits(std::vector<std::vector<int>> split_arr) {
  std::vector<int> RealSplit;
  while (split_arr.size() > 1) {
    int max_val = 0, max_i = 0;
    for (int j=0; j<split_arr.size(); j++) {
      if (split_arr[j].size()>0) {
        if (split_arr[j][0] > max_val) { max_val = split_arr[j][0]; max_i = j; }
      }
      else split_arr.erase(split_arr.begin()+j);
    }
    if (split_arr.size() < 2) break;
    bool Pair = false;
    for (int j=0; j<split_arr.size(); j++) {
      if (split_arr[j].size()>0) {
        if (max_val==split_arr[j][0] && max_i!=j) {
          Pair=true;
          RealSplit.push_back(split_arr[j][0]);
          split_arr[max_i].erase(split_arr[max_i].begin());
          split_arr.erase(split_arr.begin()+j);
          break;
        }
      } else std::cout << "Something is wrong with CountRealSplits!" << std::endl;
    }
    if (Pair==false) split_arr[max_i].erase(split_arr[max_i].begin());
  }
  return RealSplit;
}

//List the ancestor splittings of a given parton
std::vector<int> ListSplits(std::vector<Parton> parton_list, int ip) {
  std::vector<int> Splits;
  int i = ip;
  if (abs(parton_list[i].stat())==23) return Splits;
  do {
    int mom = parton_list[i].mom1();
    if (parton_list[mom].d1() != parton_list[mom].d2()) {
      Splits.push_back(mom);
      i = mom;
    } else i = mom;
  } while (abs(parton_list[i].stat())!=23);
  return Splits;
}


//Count the number of resolved splittings by the medium
int CountSplitsInMedium(std::vector<Parton> parton_list, std::vector<int> splits) {
  double LinGeV=4./0.193;
  double qhatinGeV=0.3;
  int nSplitIn = 0;
  for (unsigned int ip = 0; ip < splits.size(); ip++) {
    int a = splits[ip];
    int r = parton_list[a].dippart();
    FourVector xSplit_beg = parton_list[a].x();
    FourVector xSplit_end = parton_list[parton_list[a].d1()].x();
    FourVector xSpect_end = parton_list[parton_list[r].d1()].x();
    double rSplit = xSplit_end.p3abs();
    double lDip = (xSplit_end-xSpect_end).p3abs();
    double tDip = xSplit_end.t()-xSplit_beg.t(); //FIXME use the true time which the dipole spent in the medium
    double lRes = std::sqrt(12./qhatinGeV/tDip);
    if (rSplit < LinGeV && lDip < lRes && lRes < LinGeV) nSplitIn++;
  }
  return nSplitIn + 1; //We are interested in not the splitting but the number of resulted partons.
}



} // end namespace Adkoda
