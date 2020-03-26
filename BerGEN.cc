#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include "BerGEN.h"
#include "Util.h"

namespace Adkoda {

BerGEN::BerGEN(std::string input_file) {
  DATA = read_in_parameters(input_file);
}

void BerGEN::init() {

  inpartons = new InPartons(DATA);
  shower    = new Shower(DATA);

  remove("test/FinalPartons.out");
  remove("test/test.out");
  remove("test/test_LundPlane.out");
  return;
}

void BerGEN::next() {

  shower->init(*inpartons);
  shower->run();
  shower->print();


  //TEST MODULE: Save final particles (px, py, pz, E) for FastJet
  event_weight = shower->get_event_weight();
  parton_list = shower->get_parton_list();
  std::ofstream outfile;
  if (!outfile.is_open()) outfile.open("test/FinalPartons.out", std::ios_base::app);
  double total_p[4]={0.};
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    if (parton_list[ip].stat()<0 && abs(parton_list[ip].id())!=11) continue;
    total_p[0]+=parton_list[ip].px();
    total_p[1]+=parton_list[ip].py();
    total_p[2]+=parton_list[ip].pz();
    total_p[3]+=parton_list[ip].e();
    outfile << parton_list[ip].px() << " "
		<< parton_list[ip].py() << " "
		<< parton_list[ip].pz() << " "
		<< parton_list[ip].e() << " "
    << parton_list[ip].id() << " "
    << parton_list[ip].col() << " " << parton_list[ip].acol() << " "
    << event_weight
    << std::endl;
  }
  outfile << "# end" << std::endl;
  //Test energy-momentum conservation.
  //for (unsigned a=0; a<4; a++) std::cout << " total comp " << a << " = " << total_p[a] << std::endl;

/*
  //TEST MODULE: Save Lund Plane
  parton_list = shower->get_parton_list();
  std::ofstream OutputFile;
  if (!OutputFile.is_open()) OutputFile.open("test/test_LundPlane.out", std::ios_base::app);
  int m1 = -1;
  int m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].stat() == -23 && m1 == -1) m1 = ip;
    if (parton_list[ip].stat() == -23 && ip != m1 && m2 == -1) {m2 = ip; break;}
  }
  if (parton_list[m1].stat()<0){ //It decays
    int mom = m1;
    int d1 = parton_list[mom].d1();
    int d2 = parton_list[mom].d2();
    if (d1==d2) mom = m2; //Find the initial splitter
    double z_cut = 0.;
    double beta = 0.;
    double R = M_PI;
    while (parton_list[mom].stat()<0) {
      d1 = parton_list[mom].d1();
      d2 = parton_list[mom].d2();
      if (abs(parton_list[d1].stat()) == 51) { //A splitter
        if (parton_list[d1].pt() < parton_list[d2].pt()) { //Follow the primary
          int temp = d1;
          d1 = d2;
          d2 = temp;
        }
        double Delta = std::sqrt(std::pow(parton_list[d1].rapidity() - parton_list[d2].rapidity(), 2.) +
                                 std::pow(parton_list[d1].phi() - parton_list[d2].phi(), 2.));
        double z     = parton_list[d2].pt() / (parton_list[d1].pt() + parton_list[d2].pt());
        if (z > z_cut * pow(Delta/R, beta) && Delta < R) { // Soft Drop condition.
          double kt      = parton_list[d2].pt() * Delta;
          double mass    = Util::m(parton_list[d1].p(), parton_list[d2].p());
          double qtilde2 = mass*mass / z / (1.-z);
          double energy  = parton_list[d1].e() + parton_list[d2].e();
          double tf      = 2. * energy / mass / mass;
          std::cout << "mom: " << mom << "\t d1: " << d1 << "\t d2: " << d2 << "\t z: " << z << "\t Delta:" << Delta << "\t kt2: " << kt*kt << "\t pt2:" << z*(1.-z)*mass*mass << "\t m2: " << mass*mass << endl;
          OutputFile << qtilde2 << " " << z << " " << energy << " " << kt*kt << " "
                     << Delta << " " << mass*mass << " " << tf << " " << "1." << std::endl;
        }
        mom = d1;
      }
      else if (parton_list[d1].stat()>0) break;
      else mom = d1;
    }
  }

*/

}



InitData BerGEN::read_in_parameters(std::string input_file) {

  InitData parameter_list;

  std::string tempinput;

  // # of events
  int temp_number_events = 1;
  tempinput = Util::StringFind4(input_file, "number_events");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_number_events;
  parameter_list.number_events = temp_number_events;

  // Choice of parton shower evolution variable
  int temp_evol_scale = 0;
  tempinput = Util::StringFind4(input_file, "evol_scale");
  if (tempinput != "empty") std::stringstream(tempinput) >> temp_evol_scale;
  parameter_list.evol_scale = temp_evol_scale;

  // Switch for parton gun
  bool temp_parton_gun = true;
  tempinput = Util::StringFind4(input_file, "parton_gun");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_parton_gun;
  parameter_list.parton_gun = temp_parton_gun;

  // Type of partons in parton gun
  int temp_hard_partons = 0;
  tempinput = Util::StringFind4(input_file, "hard_partons");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_hard_partons;
  parameter_list.hard_partons = temp_hard_partons;

  // Initial COM of parton gun
  double temp_pt_max = 100.0;
  tempinput = Util::StringFind4(input_file, "pt_max");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_pt_max;
  parameter_list.pt_max = temp_pt_max;

  // Pt cutoff of the shower
  double temp_pt_min = 1.0;
  tempinput = Util::StringFind4(input_file, "pt_min");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_pt_min;
  parameter_list.pt_min = temp_pt_min;

  return parameter_list;

}

} // end namespace Adkoda
