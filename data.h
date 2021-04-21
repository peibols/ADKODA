#ifndef DATA_H
#define DATA_H

typedef struct init_data {

  int number_events;

  int evol_scale;
  int shower_kernel;

  int parton_gun;
  int hard_partons;	// 0 for gg, 1 for qq, else for qg

  double pt_max; //COM energy
  double pt_min; //Minimal scale in shower

  int medium;
  bool do_quenching;
  bool do_third_stage;
  double xmin_med;
  double eps_med;
  double L0;  // in GeV-1
  double qhat0; // in GeV^3
  double T0; // in GeV
  double alphas_med;

} InitData;

#endif  //  DATA_H
