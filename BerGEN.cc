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

  remove("FinalPartons.out");
  remove("test.out");

  return;
}

void BerGEN::next() {

  shower->init(*inpartons);
  shower->run();

  // Output
  shower->print();
  parton_list = shower->get_parton_list();
  std::ofstream outfile;
  if (!outfile.is_open()) outfile.open("FinalPartons.out", std::ios_base::app);
  double total_p[4]={0.};
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    if (parton_list[ip].stat()<0) continue;
    total_p[0]+=parton_list[ip].px();
    total_p[1]+=parton_list[ip].py();
    total_p[2]+=parton_list[ip].pz();
    total_p[3]+=parton_list[ip].e();
    outfile << parton_list[ip].px() << " "
		<< parton_list[ip].py() << " "
		<< parton_list[ip].pz() << " "
		<< parton_list[ip].e() << " "
    << std::endl;
  }
  outfile << "# end" << std::endl;
  //for (unsigned a=0; a<4; a++) std::cout << " total comp " << a << " = " << total_p[a] << std::endl;
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
