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

  bool do_quenching;
  double xmin_med;
  double eps_med;
  double L_med;  // in fm
  double qhat; // in GeV^2/fm
  double alphas_med;

} InitData;

#endif  //  DATA_H
