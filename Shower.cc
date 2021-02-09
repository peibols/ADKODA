#include <iostream>

#include "Shower.h"
#include "Util.h"

namespace Adkoda {

Shower::Shower(const InitData &DATA_in) : DATA(DATA_in) {
  // Introduce kernels
  //FIXME shouldnt this sum depend on the max scale?
  if (DATA.shower_kernel == 0) { //Alterelli-Parisi splitting kernels
    Pgg* pgg = new Pgg (21, 21, 21);
    kernels.push_back(pgg);

/*
    Pgq* pgq;
    Pqq* pqq;
     for (int fl = 1; fl <= 5; fl++) {
      pqq = new Pqq (fl, fl, 21);
      kernels.push_back(pqq);
      pqq = new Pqq (-fl, -fl, 21);
      kernels.push_back(pqq);
      pgq = new Pgq (21, fl, -fl);
      kernels.push_back(pgq);
    }
*/

  } else if (DATA.shower_kernel == 1) { //Catani-Seymour splitting kernels
    Pgg_CS* pgg = new Pgg_CS (21, 21, 21);
    kernels.push_back(pgg);
    Pgq_CS* pgq;
    Pqq_CS* pqq;
     for (int fl = 1; fl <= 5; fl++) {
      pqq = new Pqq_CS (fl, fl, 21);
      kernels.push_back(pqq);
      pqq = new Pqq_CS (-fl, -fl, 21);
      kernels.push_back(pqq);
      pgq = new Pgq_CS (21, fl, -fl);
      kernels.push_back(pgq);
    }
  } else std::cout << "ERROR in bergen_input: shower_kernel = 0, 1." << std::endl;

  // Alpha_s and scale cutoff
  pt_min = DATA.pt_min;
  //std::cout << " pt_min= " << pt_min << std::endl;
  max_alpha_s = alpha_s( pt_min * pt_min );
  //std::cout << " max_alpha_s= " << max_alpha_s << std::endl;

  // Medium parameters
  alphas_med = DATA.alphas_med;
  L0 = DATA.L0; // in GeV^-1
  qhat0= DATA.qhat0; // in GeV^3

  //std::cout << "Shower CONSTRUCTED" << std::endl;

}

void Shower::init (InPartons *inpartons) {

  if (DATA.parton_gun==0 || DATA.parton_gun==1) parton_list = inpartons->PartonList();
  else {
    parton_list = inpartons->PythiaPartonList();
  }
  event_weight = inpartons->event_weight;
  event_xsec = inpartons->event_xsec;
  hard_pt_max = inpartons->hard_pt_max;

  std::cout << " Weight= " << event_weight << " Xsec= " << event_xsec << " Pt Max= " << hard_pt_max << std::endl;

  // Fix minimum and maximum scale
  //double ecms = Util::m2(parton_list[parton_list.size()-2].p(), parton_list[parton_list.size()-1].p());;
  if      (DATA.evol_scale==0) {              // pt ordering
    t_min = std::pow(DATA.pt_min, 2.);
    t_max = std::pow(hard_pt_max,2.);         // Right way for PYTHIA8 input
    //t_max = std::pow(hard_pt_max/2.,2.);
    //t_max = std::min(std::pow(DATA.pt_max, 2.), ecms);
  }
  else if (DATA.evol_scale==1) {              // m ordering
    t_min = std::pow(DATA.pt_min, 2.) * 4.;
    //t_max = std::min(std::pow(DATA.pt_max, 2.) * 4., ecms);
    t_max = std::pow(hard_pt_max/2., 2.) * 4.;
  }
  else if (DATA.evol_scale==2) {              // tf ordering
    t_min = std::pow(DATA.pt_min, 2.) * 2. / (hard_pt_max/2.);
    //t_max = std::min(std::pow(DATA.pt_max, 2.) * 2. / (DATA.pt_max/2.), ecms/2./(DATA.pt_max/2.));
    t_max = std::pow(hard_pt_max/2., 2.) * 2. / DATA.pt_min;
  }
  else if (DATA.evol_scale==3) {              // qt ordering
    t_min = std::pow(DATA.pt_min, 2.) * 16.;
    //t_max = std::min(std::pow(DATA.pt_max, 2.) * 16., 4.*ecms);
    t_max = std::pow(hard_pt_max/2., 2.) * 16.;
  }

  // Update max_color index
  max_colour = 101;
  for (unsigned int ip=0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].col()  > max_colour) max_colour = parton_list[ip].col();
    if (parton_list[ip].acol() > max_colour) max_colour = parton_list[ip].acol();
  }

  //std::cout << "Shower INITIALIZED" << std::endl;
  return;

}

void Shower::run () {

  std::cout << "Initial Parton List size = " << parton_list.size() << endl;
  print();
  std::cout << "Shower RUNNING" << std::endl;

  // Vacuum Shower
  bool do_evolve = 1;
  while (do_evolve) do_evolve = evolve(); // FIXME how does this evolve() called?!

  std::cout << "Shower FINISHED" << std::endl;
  std::cout << "Final Parton List size = " << parton_list.size() << endl;
  print();

}

double Shower::beta0(int nf) { return 11./6.*CA - 2./3.*TR*nf; }

double Shower::beta1(int nf) { return 17./6.*CA*CA - (5./3.*CA + CF)*TR*nf; }

double Shower::alpha_s0(double t) { // FIXME use a particle list for the mass
  double tref, asref, b0;
  double mb2  = std::pow(4.75, 2.);
  double mc2  = std::pow(1.3, 2.);
  double asmz = 0.118;
  double tmin = pt_min*pt_min;
  if (t >= mb2) {
    tref  = MZ2;
    asref = asmz;
    b0    = beta0(5)/(2.*M_PI);
  }
  else if (mb2 > t && t >= mc2) {
    tref  = mb2;
    asref = 1. / ( 1./asmz + beta0(5)/(2.*M_PI)*std::log(mb2/MZ2) );
    b0    = beta0(4)/(2.*M_PI);
  }
  else if (mc2 > t && t > tmin) {
    tref  = mc2;
    asref = 1. / ( 1./asmz + beta0(5)/(2.*M_PI)*std::log(mb2/MZ2) + beta0(4)/(2.*M_PI)*std::log(mc2/mb2));
    b0    = beta0(3)/(2.0*M_PI);
  } else {
    tref  = mc2;
    asref = 1. / ( 1./asmz + beta0(5)/(2.*M_PI)*std::log(mb2/MZ2) + beta0(4)/(2.*M_PI)*std::log(mc2/mb2));
    b0    = beta0(3)/(2.0*M_PI);
    return 1. / ( 1./asref + b0*std::log(tmin/tref) );
  }
  return 1. / ( 1./asref + b0*std::log(t/tref) );
}

double Shower::alpha_s(double t) { // FIXME use a particle list for the mass
  double tref, asref, b0, b1, wr, w;
  double mb2  = std::pow(4.75, 2.);
  double mc2  = std::pow(1.3, 2.);
  double asmz = 0.118;
  double tmin = pt_min*pt_min;
  if (t >= mb2) {
    tref  = MZ2;
    asref = asmz;
    b0    = beta0(5)/(2.*M_PI);
    b1    = beta1(5)/std::pow(2.*M_PI, 2.);
  }
  else if (mb2 > t && t >= mc2) {
    tref  = mb2;
    wr    = 1. + beta0(5) / (2.*M_PI) * asmz * std::log(mb2/MZ2);
    asref = asmz / wr * ( 1. - beta1(5)/std::pow(2.*M_PI, 2.)/(beta0(5)/(2.*M_PI))*asmz*std::log(wr)/wr );
    b0    = beta0(4) / (2.*M_PI);
    b1    = beta1(4) / std::pow(2.*M_PI, 2.);
  }
  else if (mc2 > t && t > tmin) {
    tref  = mc2;
    double wr0    = 1. + beta0(5) / (2.*M_PI) * asmz * std::log(mb2/MZ2);
    double asref0 = asmz / wr0 * ( 1. - beta1(5)/std::pow(2.*M_PI, 2.)/(beta0(5)/(2.*M_PI))*asmz*std::log(wr0)/wr0 );
    wr    = 1. + beta0(4) / (2.*M_PI) * asref0 * std::log(mc2/mb2);
    asref = asref0 / wr * ( 1. - beta1(4)/std::pow(2.*M_PI, 2.)/(beta0(4)/(2.*M_PI))*asref0*std::log(wr)/wr );
    b0    = beta0(3) / (2.0*M_PI);
    b1    = beta1(3) / std::pow(2.*M_PI, 2.);
  } else {
    tref  = mc2;
    double wr0    = 1. + beta0(5) / (2.*M_PI) * asmz * std::log(mb2/MZ2);
    double asref0 = asmz / wr0 * ( 1. - beta1(5)/std::pow(2.*M_PI, 2.)/(beta0(5)/(2.*M_PI))*asmz*std::log(wr0)/wr0 );
    wr    = 1. + beta0(4) / (2.*M_PI) * asref0 * std::log(mc2/mb2);
    asref = asref0 / wr * ( 1. - beta1(4)/std::pow(2.*M_PI, 2.)/(beta0(4)/(2.*M_PI))*asref0*std::log(wr)/wr );
    b0    = beta0(3) / (2.0*M_PI);
    b1    = beta1(3) / std::pow(2.*M_PI, 2.);
    w     = 1. + b0 * asref * std::log(tmin/tref);
    return asref / w * ( 1. - b1 / b0 * asref * std::log(w)/w );
  }
  w = 1. + b0*asref*std::log(t/tref);
  return asref / w * ( 1. - b1 / b0 * asref * std::log(w)/w );
}

/*
double Shower::alpha_s( double t ) {
  double Nf = 5.;
  double alpha;
  if (t > 0.95) alpha = 12.0 * M_PI / ( 33.0 - 2.0 * Nf ) / std::log( t / ( Lambda * Lambda ) );
  else          alpha = 12.0 * M_PI / ( 33.0 - 2.0 * Nf ) / std::log( 0.95 / ( Lambda * Lambda ) );
  return alpha;
}
*/

void Shower::print() {
  std::cout << "#Event weight: " << event_weight << std::endl;
  std::cout << "#ip\t ID\t Stat\t m1\t m2\t d1\t d2\t c\t ac\t px\t\t py\t\t pz\t\t E\t\t m\t\t"
  //<< "pt2\t\t"
  //<< "x\t\t y\t\t z\t\t t"
  << std::endl;
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    std::cout << ip << "\t ";
    parton_list[ip].display();
  }
  return;
}

} // end namespace Adkoda
