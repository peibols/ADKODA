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

using namespace Util;

namespace Adkoda {

Tests::Tests() {
  
  remove("test/test_weights.out");
  remove("test/test_veto.out");
  remove("test/test_FinalPartons.out");
  remove("test/test_LundPlane_history.out");
  remove("test/test_LundPlane_history_weight.out");
  remove("test/test_LundPlane_FJ.out");
  remove("test/test_LundPlane_FJ_weight.out");
  remove("test/test_JetMass_FJ.out");
  remove("test/test_IteratedSoftDrop_FJ.out");
  remove("test/test_kinematics.out");

  if (!outfile_test_weights.is_open()) outfile_test_weights.open("test/test_weights.out", std::ios_base::app);
  if (!outfile_test_FinalPartons.is_open()) outfile_test_FinalPartons.open("test/test_FinalPartons.out", std::ios_base::app);
  if (!outfile_test_FirstSplit.is_open()) outfile_test_FirstSplit.open("test/test_FirstSplit.out", std::ios_base::app);
  if (!outfile_test_LundPlane_history.is_open()) outfile_test_LundPlane_history.open("test/test_LundPlane_history.out", std::ios_base::app);
  if (!outfile_test_LundPlane_history_weight.is_open()) outfile_test_LundPlane_history_weight.open("test/test_LundPlane_history_weight.out", std::ios_base::app);
  if (!outfile_test_LundPlane_FJ.is_open()) outfile_test_LundPlane_FJ.open("test/test_LundPlane_FJ.out", std::ios_base::app);
  if (!outfile_test_LundPlane_FJ_weight.is_open()) outfile_test_LundPlane_FJ_weight.open("test/test_LundPlane_FJ_weight.out", std::ios_base::app);
  if (!outfile_test_JetMass_FJ.is_open()) outfile_test_JetMass_FJ.open("test/test_JetMass_FJ.out", std::ios_base::app);
  if (!outfile_test_IteratedSoftDrop_FJ.is_open()) outfile_test_IteratedSoftDrop_FJ.open("test/test_IteratedSoftDrop_FJ.out", std::ios_base::app);

  max_nx = 100;
  x_bin = 1./double(max_nx);
  for (unsigned a=0; a<100; a++) {
    for (unsigned b=0; b<6; b++) {
      dx_hist[a][b]=0.;
    }
  }

  max_kx = 100;
  kx_bin = 20./double(max_kx);
  for (unsigned a=0; a<101; a++) {
    for (unsigned b=0; b<6; b++) {
      kx_hist[a][b]=0.;
    }
  }

  ntau_cuts = 6;
  tau_cuts[0] = 0.01;
  tau_cuts[1] = 0.02;
  tau_cuts[2] = 0.1;
  tau_cuts[3] = 0.2;
  tau_cuts[4] = 1.;
  tau_cuts[5] = 2.;

}

void Tests::run(std::vector<Parton> parton_list, double event_weight, double event_xsec, int iEv) {

  //Test_Weights(parton_list, event_xsec, event_weight, outfile_test_weights);
  Test_PrintFinalPartons(parton_list, event_xsec, event_weight, outfile_test_FinalPartons);
  //Test_EnergyMomentumConservation(parton_list);
  //Test_FirstSplit(parton_list, event_xsec, event_weight, outfile_test_FirstSplit);
  //Test_PrintLundPlane_history(parton_list, event_xsec, event_weight, outfile_test_LundPlane_history, outfile_test_LundPlane_history_weight);
  //Test_PrintLundPlane_FJ(parton_list, event_xsec, event_weight, outfile_test_LundPlane_FJ, outfile_test_LundPlane_FJ_weight);
  //Test_JetMass_FJ(parton_list, event_xsec, event_weight, outfile_test_JetMass_FJ);
  //Test_IteratedSoftDrop_FJ(parton_list, event_xsec, event_weight, outfile_test_IteratedSoftDrop_FJ);

  CascadeDist(parton_list);

}

void Tests::close(int nEv) {
  
  outfile_test_weights.close();
  outfile_test_FinalPartons.close();
  outfile_test_FirstSplit.close();
  outfile_test_LundPlane_history.close();
  outfile_test_LundPlane_history_weight.close();
  outfile_test_LundPlane_FJ.close();
  outfile_test_LundPlane_FJ_weight.close();
  outfile_test_JetMass_FJ.close();
  outfile_test_IteratedSoftDrop_FJ.close();

  // Print Dx hist
  std::ofstream dxfile("Dx_hist.dat");
  for (int a=0; a<max_nx; a++) {
    double xcen = x_bin*double(a)+x_bin/2.;
    dxfile << xcen << " ";
    for (int b=0; b<ntau_cuts; b++) {
      dxfile << dx_hist[a][b]/x_bin*xcen*sqrt(xcen)/nEv << " " << sqrt(dx_hist[a][b])/x_bin*xcen*sqrt(xcen)/nEv << " ";
    }
    dxfile << endl;
  }
  dxfile.close();

  // Print Kx hist
  std::ofstream kxfile("Kx_hist.dat");
  for (int a=0; a<=max_kx; a++) {
    double xcen = -10.1+kx_bin*double(a)+kx_bin/2.;
    kxfile << xcen << " ";
    for (int b=0; b<ntau_cuts; b++) {
      kxfile << kx_hist[a][b]/kx_bin/nEv << " " << sqrt(kx_hist[a][b])/kx_bin/nEv << " ";
    }
    kxfile << endl;
  }
  kxfile.close();

}

void Tests::CascadeDist(std::vector<Parton> parton_list) {

  double pplus;
  double alphas_med = 0.3;
  double qhat = 1.5 * 0.1973;
  double angle, k[3];
  for (unsigned int i=0; i<parton_list.size(); i++) {

    Parton p = parton_list[i];
    if (p.stat()==-23) { 
      FourVector pIni = p.p();
      AlignWithZ(pIni, angle, k);
      pplus=1./2.*(pIni.t()+pIni.z());
      continue;
    }
    if (p.x().t()==0.) continue;
    
    double ti = p.x().t();
    double taui = alphas_med * std::sqrt(qhat/pplus) * ti / 0.1973;

    double tf = p.xf().t();
    if (tf==-1000.) tf = 10000000.; 
    double tauf = alphas_med * std::sqrt(qhat/pplus) * tf / 0.1973;

    //cout << " pplus= " << pplus << " taui= " << taui << " tauf= " << tauf << endl;

    for (int c=0; c<ntau_cuts; c++) {
      if (taui<=tau_cuts[c] && tauf>tau_cuts[c]) {
        FourVector pVec = p.p();
        Rotation(pVec, angle, k);
        // x hist
        double t_pplus = 1./2.*(pVec.t()+pVec.z());
        double x = t_pplus/pplus;
        int nx = int(x/x_bin);
        //cout << " x= " << x << " nx= " << nx << endl;
        if (nx>=0 && nx<max_nx) {
          dx_hist[nx][c]+=1.;
        }
        // kx hist
        double kx = pVec.x();
        double sign=0.;
        if (kx>0.) sign=1.;
        else sign=-1.;
        int nkx = max_kx/2+int((kx+sign*kx_bin/2.)/kx_bin);
        if (nkx>=0 && nkx<=max_kx) {
          kx_hist[nkx][c]+=1.;
        }
      }
    }      
    
  }

}

/*
//TEST MODULE: First mass of jet according to some ad-hoc definition	
void Tests::Test_FirstMass(std::vector<Parton> parton_list, std::ofstream &outfile, double mass_hist[2][20], double &njets) {

  std::vector<PseudoJet> particles;
  for (unsigned int i=0; i<parton_list.size(); i++) {
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
    //int orig=-1000;
    double highest_splitmass=0.;
    while (true) {
      if (parton_list[anc].stat()==23) {
        cout << " hard never splitted! " << endl;
        break;
      }
      if (parton_list[anc].stat()==-23) {
        //orig = anc;
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
        //double deltaR=p1.delta_R(p2);
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
*/

//TEST MODULE: Save the initiator e-e+ kinematics and event_weights
void Tests::Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  outfile << parton_list[4].e()  << " " << parton_list[4].px() << " " <<
             parton_list[4].py() << " " << parton_list[4].pz() << " " <<
             event_xsec << " " << event_weight << std::endl;
}

//TEST MODULE: Save final particles (px, py, pz, E)
void Tests::Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
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
void Tests::Test_EnergyMomentumConservation(std::vector<Parton> parton_list) {
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
void Tests::Test_FirstSplit(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  //Find primaries
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hard-scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && int(ip) != m1 && m2 == -1) { m2 = ip; break; }
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
  //FourVector xSpect_end = parton_list[parton_list[parton_list[m1].dippart()].d1()].x();
  double tf = xSplit_end.t()-xSplit_beg.t();
  outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " " <<
             Delta << " " << mass2 << " " << tf << " " <<
             event_xsec << " " << event_weight << " " << nSplitIn << " " <<
             pt_ev << std::endl;
}

//TEST MODULE: Save the primary Lund plane using the history and FastJet definition of the variables
void Tests::Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight) {
  //Find primary splitter
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hard-scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && int(ip) != m1 && m2 == -1) { m2 = ip; break; }
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
        //FourVector xSpect_end = parton_list[parton_list[parton_list[mom].dippart()].d1()].x();
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
void Tests::Test_PrintLundPlane_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight) {
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
void Tests::Test_JetMass_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
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
void Tests::Test_IteratedSoftDrop_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
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
std::vector<int> Tests::ListFinalDescendants(std::vector<Parton> parton_list, int ip) {
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
std::vector<int> Tests::CountRealSplits(std::vector<std::vector<int>> split_arr) {
  std::vector<int> RealSplit;
  while (split_arr.size() > 1) {
    int max_val = 0, max_i = 0;
    for (unsigned int j=0; j<split_arr.size(); j++) {
      if (split_arr[j].size()>0) {
        if (split_arr[j][0] > max_val) { max_val = split_arr[j][0]; max_i = j; }
      }
      else split_arr.erase(split_arr.begin()+j);
    }
    if (split_arr.size() < 2) break;
    bool Pair = false;
    for (unsigned int j=0; j<split_arr.size(); j++) {
      if (split_arr[j].size()>0) {
        if (max_val==split_arr[j][0] && max_i!=int(j)) {
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
std::vector<int> Tests::ListSplits(std::vector<Parton> parton_list, int ip) {
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
int Tests::CountSplitsInMedium(std::vector<Parton> parton_list, std::vector<int> splits) {
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
