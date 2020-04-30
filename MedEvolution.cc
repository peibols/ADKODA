
#include <iostream>
#include <fstream>

#include "Cascade.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

bool Cascade::evolve(double &start_time, Parton p, std::vector<Parton> &cascade_list, int mom)
{
  
  double t_ini = start_time;
  double t = 10000000.;

  double xi=p.xfrac();

  int wKernel;

  for (int iKernel=0; iKernel<kernels.size(); iKernel++) {
    if (kernels[iKernel]->flav(0) != p.id()) continue; // Skip if kernel not applies

    // Generate next t
    double f = kernels[iKernel]->Integral(eps_med, 1.-eps_med);
    double tt =  t_ini - std::log(dis(gen)) * std::sqrt(xi) / f;
    //std::cout << " f= " << f << std::endl;

    if ( tt < t ) {
      t 	= tt;
      wKernel	= iKernel;
    }
  }

  double time = tkin(t);

  double form_time = time - p.x().t();
  //std::cout << " form time= " << form_time << std::endl << std::endl;

  FourVector pLab = p.p();
  Rotation(pLab, -angle, k);
  FourVector x;
  x.Set(pLab.x()/pLab.p3abs()*form_time,
	pLab.y()/pLab.p3abs()*form_time,
	pLab.z()/pLab.p3abs()*form_time,
	form_time);
 
  FourVector xSplit;
  xSplit = x + p.x();
  double split_pos = xSplit.p3abs() / 0.1973; // in GeV^-1
  
  //std::cout << " time = " << time << " split_pos = " << split_pos * 0.1973 << std::endl; 

  if (split_pos > L_med) { // Finish attempts if outside of medium
    std::cout << "Outside medium." << std::endl;
    return 0;  
  }
  
  // Generate final z value
  double z = kernels[wKernel]->GenerateZ(eps_med, 1.-eps_med, dis(gen));
  
  // Accept / Reject veto procedure (assume fixed medium coupling FTM)
  double f = kernels[wKernel]->Value(z,0.);
  double g = kernels[wKernel]->Estimate(z);
  if (f/g > 1. ) {
    std::cout << " f/g= " << f/g << std::endl;
    exit(1);
  }
  if (f/g < dis(gen) ) {
    start_time = t;
    //std::cout << " Failed attempt. Going to next t." << std::endl;
    return 1;
  }

  // Add daughters
  //std::cout << " Adding cascade daughters at t_kin= " << time << " and tau= " << t << " and z= " << z << std::endl;
  //std::cout << " with z = " << z << std::endl;
  //std::cout << " xi= " << xi << std::endl;

  FourVector pDau1, pDau2;
  pDau1.Set(p.px()*z,p.py()*z,p.pz()*z,p.e()*z);
  pDau2.Set(p.px()*(1.-z),p.py()*(1.-z),p.pz()*(1.-z),p.e()*(1.-z));

  int dau1_id=21;
  int dau2_id=21;

  int dau1_stat=101;
  int dau2_stat=101;
  
  if (pDau1.t() < med_cutoff ) {
    dau1_stat=-102;
    //std::cout << " Thermalised d1 e= " << pDau1.t() << std::endl;
  }
  if (pDau2.t() < med_cutoff ) {
    dau2_stat=-102;
    //std::cout << " Thermalised d2 e= " << pDau2.t() << std::endl;
  }

  Parton daughter1(Parton(dau1_id, dau1_stat, pDau1, x+p.x())); //Status: 101, active cascade particles
  Parton daughter2(Parton(dau2_id, dau2_stat, pDau2, x+p.x()));
  
  daughter1.set_mom1(mom), daughter1.set_mom2(0);
  daughter2.set_mom1(mom), daughter2.set_mom2(0);
  
  double d1_xfrac = z * p.xfrac();
  double d2_xfrac = (1.-z) * p.xfrac();
  daughter1.set_xfrac(d1_xfrac);
  daughter2.set_xfrac(d2_xfrac);
  
  cascade_list.push_back(daughter1);
  int Dau1 = int(cascade_list.size()) - 1;
  cascade_list.push_back(daughter2);
  int Dau2 = int(cascade_list.size()) - 1;

  // Update splitter
  cascade_list[mom].set_stat(-p.stat()); //It decayed: state -> decayed
  cascade_list[mom].set_d1(Dau1);
  cascade_list[mom].set_d2(Dau2);

  return 0;
}

} // end namespace Adkoda

