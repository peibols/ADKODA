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
  cascade   = new Cascade(DATA);

  return;
}

void BerGEN::next() {

  shower->init(inpartons);
  event_weight = shower->get_event_weight();
  event_xsec = shower->get_event_xsec();
  shower->run();
  parton_list = shower->get_parton_list();
  if (DATA.do_quenching) {
    cascade->init(parton_list);
    cascade->run();
    parton_list = cascade->get_parton_list();
  }
  //shower->print();
  if (DATA.do_third_stage) {
    shower->third_stage_init();
    shower->run();
  }
  parton_list = shower->get_parton_list();

}

void BerGEN::print() { //FIXME call this from Shower as a friend
  std::cout << "#Event weight: " << event_weight << std::endl;
  std::cout << "#ip\t ID\t Stat\t m1\t m2\t d1\t d2\t c\t ac\t px\t py\t pz\t E\t m\t pt2\t x\t y\t z\t t" << std::endl;
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
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
  int temp_parton_gun = true;
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

  // Medium is applied or not (to chop vacuum shower)
  int temp_medium = 0;
  tempinput = Util::StringFind4(input_file, "medium");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_medium;
  parameter_list.medium = temp_medium;

  // Jet transport parameter in GeV3
  double temp_qhat0 = 0.;
  tempinput = Util::StringFind4(input_file, "qhat0");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_qhat0;
  parameter_list.qhat0 = temp_qhat0;

  // Medium fixed size in GeV-1
  double temp_L0 = 0.;
  tempinput = Util::StringFind4(input_file, "L0");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_L0;
  parameter_list.L0 = temp_L0;

  // Medium fixed temperature in GeV
  double temp_T0 = 0.;
  tempinput = Util::StringFind4(input_file, "T0");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_T0;
  parameter_list.T0 = temp_T0;

  // alphas med
  double temp_alphas_med = 0.3;
  tempinput = Util::StringFind4(input_file, "alphas_med");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_alphas_med;
  parameter_list.alphas_med = temp_alphas_med;
  
  // eps med
  double temp_eps_med = 0.01;
  tempinput = Util::StringFind4(input_file, "eps_med");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_eps_med;
  parameter_list.eps_med = temp_eps_med;

  // xmin med
  double temp_xmin_med = 0.01;
  tempinput = Util::StringFind4(input_file, "xmin_med");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_xmin_med;
  parameter_list.xmin_med = temp_xmin_med;
  
  // Do medium cascade or not
  bool temp_do_quenching = false;
  tempinput = Util::StringFind4(input_file, "do_quenching");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_do_quenching;
  parameter_list.do_quenching = temp_do_quenching;

  // Do medium cascade or not
  bool temp_do_third_stage = false;
  tempinput = Util::StringFind4(input_file, "do_third_stage");
  if (tempinput != "empty") std::istringstream(tempinput) >> temp_do_third_stage;
  parameter_list.do_third_stage = temp_do_third_stage;

  return parameter_list;

}

} // end namespace Adkoda
