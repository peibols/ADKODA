#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include "BerGEN.h"
#include "Util.h"

namespace Adkoda {

BerGEN::BerGEN(std::string input_file) { DATA = read_in_parameters(input_file); }

void BerGEN::init() {

  inpartons = new InPartons(DATA);
  shower    = new Shower(DATA);

  return;
}

void BerGEN::next() {

  shower->init(*inpartons);
  event_weight = shower->get_event_weight();
  event_xsec = shower->get_event_xsec();
  shower->run();
  parton_list = shower->get_parton_list();
  //shower->print();

}

void BerGEN::print() { //FIXME call this from Shower as a friend
  std::cout << "#Event weight: " << event_weight << std::endl;
  std::cout << "#ip\t ID\t Stat\t m1\t m2\t d1\t d2\t c\t ac\t px\t py\t pz\t E\t m\t pt2\t x\t y\t z\t t" << std::endl;
  for (int ip=0; ip<parton_list.size(); ip++) {
    std::cout << ip << "\t ";
    parton_list[ip].display();
  }
  return;
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

  // Splitting kernel of the shower
  int temp_shower_kernel = 0;
  tempinput = Util::StringFind4(input_file, "shower_kernel");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_shower_kernel;
  parameter_list.shower_kernel = temp_shower_kernel;

  return parameter_list;

}

} // end namespace Adkoda
