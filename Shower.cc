#include <iostream>

#include "Shower.h"

namespace Adkoda {

Shower::Shower(const InitData &DATA_in) : DATA(DATA_in) {

  // Introduce kernels
  Pgg* pgg = new Pgg (21, 21, 21);
  kernels.push_back(pgg);
  //Pqg* pqg = new Pqg (1, 1, 21);
  //kernels.push_back(pqg);
  //Pqq* pqq = new Pqq (21, 1, -1);
  //kernels.push_back(pqq);
  //FIXME by including more flavors, will q/g ratio change?
/*   for (int fl = 1; fl <= 6; fl++) {
    Pqg* pqg = new Pqg (fl, fl, 21);
    kernels.push_back(pqg);
    Pqg* pqg = new Pqg (-fl, -fl, 21);
    kernels.push_back(pqg);
    Pqq* pqq = new Pqq (21, fl, -fl);
    kernels.push_back(pqq);
  }
*/

  // Alpha_s and scale cutoff
  pt_min = DATA.pt_min;
  std::cout << " pt_min= " << pt_min << std::endl;
  max_alpha_s = alpha_s( pt_min * pt_min );
  std::cout << " max_alpha_s= " << max_alpha_s << std::endl;

  // Random number generator
  dis = std::uniform_real_distribution<double> { 0.0, 1.0 };
  std::random_device rd;
  std::mt19937 gen(rd());

  std::cout << "Shower CONSTRUCTED" << std::endl;

}

void Shower::init ( InPartons inpartons) {

  parton_list = inpartons.PartonList();

  // Fix minimum and maximum scale
  if      (DATA.evol_scale==0) {              // pt ordering
    t_min = std::pow(DATA.pt_min, 2.0);
    t_max = std::pow(DATA.pt_max, 2.0)/4.0;
  }
  else if (DATA.evol_scale==1) {              // m ordering
    t_min = std::pow(DATA.pt_min, 2.0) * 4.0;
    t_max = std::pow(DATA.pt_max, 2.0);
  }
  else if (DATA.evol_scale==2) {              // tf ordering
    t_min = 2.0 * DATA.pt_min*DATA.pt_min / (DATA.pt_max/2.0);
    t_max = DATA.pt_max;
  }
  else if (DATA.evol_scale==3) {              // qt ordering
    t_min = DATA.pt_min * DATA.pt_min * 16.0;
    t_max = std::pow(DATA.pt_max, 2.0);
  }


  // Update max_color index
  max_colour = 101;
  for (unsigned int ip=0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].col()  > max_colour) max_colour = parton_list[ip].col();
    if (parton_list[ip].acol() > max_colour) max_colour = parton_list[ip].acol();
  }

  std::cout << "Shower INITIALIZED" << std::endl;
  return;

}

void Shower::run () {

  std::cout << "Initial Parton List size = " << parton_list.size() << endl;
  std::cout << "Shower RUNNING" << std::endl;
  bool do_evolve = 1;
  while (do_evolve) do_evolve = evolve(); // FIXME how does this evolve() called?!

  std::cout << "Shower FINISHED" << std::endl;
  std::cout << "Final Parton List size = " << parton_list.size() << endl;

}


/*
double Shower::beta0(double nf) { return 11.0/6.0*CA - 2./3.*TR*nf; }

double Shower::alpha_s( double t ) { // FIXME use a particle list for the mass, also in NLO
  double tref, asref, b0;
  double mb2 = std::pow(4.75, 2.);
  double mc2 = std::pow(1.3, 2.);
  double asmz = 0.118;
  if (t >= mb2) {
    tref = MZ2;
    asref = asmz;
    b0 = beta0(5.)/(2.0*M_PI);
    return 1.0 / ( 1.0/asref + b0*std::log(t/tref) );
  }
  else if (mb2 > t && t >= mc2) {
    tref = mb2;
    asref = 1.0 / ( 1.0/asmz + beta0(5.)/(2.0*M_PI)*std::log(t/MZ2) );
    b0 = beta0(4.)/(2.0*M_PI);
    return 1.0 / ( 1.0/asref + b0*std::log(t/tref) );
  }
  else {
    tref = mc2;
    asref = 1.0 / ( 1.0/asmz + beta0(5.)/(2.0*M_PI)*std::log(t/MZ2) );
    b0 = beta0(3.)/(2.0*M_PI);
  }

}
*/

double Shower::alpha_s( double t ) {
  double alpha = 12.0 * M_PI / ( 33.0 - 2.0 * Nf ) / std::log( t / ( Lambda * Lambda ) );
  return alpha;
}



void Shower::print() {

  for (int ip=0; ip<parton_list.size(); ip++) {
    std::cout << ip << " ";
    parton_list[ip].display();
  }

  return;

}

} // end namespace Adkoda
