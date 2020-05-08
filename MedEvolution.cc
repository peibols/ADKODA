
#include <iostream>
#include <fstream>

#include "Cascade.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

bool Cascade::evolve(double &start_time, std::vector<Parton> &cascade_list)
{
  
  double t_ini = start_time;
  double tlarge=1000000000.;
  double t = tlarge;

  int wKernel;
  int wSplit;

  for (int iSplit=0; iSplit<cascade_list.size(); iSplit++) {
  
    Parton p = cascade_list[iSplit];
 
    if (p.stat()<0) continue;

    double xi=p.xfrac();

    // Splitting kernels
    for (int iKernel=0; iKernel<kernels.size(); iKernel++) {
      if (kernels[iKernel]->flav(0) != p.id()) continue; // Skip if kernel not applies

      // Generate next t
      double f = kernels[iKernel]->Integral(eps_med, 1.-eps_med);
      double tt =  t_ini - std::log(dis(gen)) * std::sqrt(xi) / f;
      //std::cout << " f= " << f << std::endl;

      if ( tt < t ) {
        t 	= tt;
        wKernel	= iKernel;
        wSplit  = iSplit;
      }
    } // end splitting kernels loop

    // Broadening
    double f = BroadInt();
    double tt = t_ini - std::log(dis(gen)) / f;
    if ( tt < t ) {
      t 	= tt;
      wKernel	= -1;	// broadening flag
      wSplit  = iSplit;
    }

  } // end cascade loop

  if (t == tlarge) {
    return 0;
  }

  if (t>0.2) return 0;

  start_time = t; // Update time

  Parton p = cascade_list[wSplit];
  
  double time = tkin(t); // in fm

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
    cascade_list[wSplit].set_stat(-103); // Outside medium stat
    std::cout << "Outside medium." << std::endl;
    return 1;  
  }
 
  if (wKernel!=-1) { // Splitting 
  
    // Generate final z value
    double z = kernels[wKernel]->GenerateZ(eps_med, 1.-eps_med, dis(gen));
    //cout << " z= " << z << endl;
    if (z<0. || z>1.) {
      cout << " wtf z= " << z << endl;
    }  
    // Accept / Reject veto procedure (assume fixed medium coupling FTM)
    double f = kernels[wKernel]->Value(z,0.);
    double g = kernels[wKernel]->Estimate(z);
    if (f/g > 1. ) {
      std::cout << " f/g= " << f/g << std::endl;
      exit(1);
    }
//  if (f/g < dis(gen) ) {
    //std::cout << " Failed attempt. Going to next t." << std::endl;
//    return 1;
//  }

    FourVector pDau1, pDau2;
    pDau1.Set(p.px()*z,p.py()*z,p.pz()*z,p.e()*z);
    pDau2.Set(p.px()*(1.-z),p.py()*(1.-z),p.pz()*(1.-z),p.e()*(1.-z));

    int dau1_id=21;
    int dau2_id=21;

    int dau1_stat=101;
    int dau2_stat=101;
  
    double d1_xfrac = z * p.xfrac();
    double d2_xfrac = (1.-z) * p.xfrac();
    if (d1_xfrac < xmin_med) {
      dau1_stat=-102;
      //std::cout << " Thermalised d1 e= " << pDau1.t() << std::endl;
    }
    if (d2_xfrac < xmin_med) {
      dau2_stat=-102;
      //std::cout << " Thermalised d2 e= " << pDau2.t() << std::endl;
    }

    Parton daughter1(Parton(dau1_id, dau1_stat, pDau1, xSplit)); //Status: 101, active cascade particles
    Parton daughter2(Parton(dau2_id, dau2_stat, pDau2, xSplit));
  
    daughter1.set_mom1(wSplit), daughter1.set_mom2(0);
    daughter2.set_mom1(wSplit), daughter2.set_mom2(0);
  
    daughter1.set_xfrac(d1_xfrac);
    daughter2.set_xfrac(d2_xfrac);
  
    cascade_list.push_back(daughter1);
    int Dau1 = int(cascade_list.size()) - 1;
    cascade_list.push_back(daughter2);
    int Dau2 = int(cascade_list.size()) - 1;

    // Update splitter
    cascade_list[wSplit].set_stat(-p.stat()); //It decayed: state -> decayed
    cascade_list[wSplit].set_xf(xSplit);
    cascade_list[wSplit].set_d1(Dau1);
    cascade_list[wSplit].set_d2(Dau2);

  }
  else {  // Broadening

    FourVector pNow = p.p();
    
    double q;
    double newpzsq;
    int negpz=0;
    while (true) {
      q = GenerateQ(dis(gen));
      newpzsq = pNow.t()*pNow.t() - q*q;
      if (newpzsq > 0) break;
      else negpz++;
      if (negpz>1000) {
        cout << " too soft for kicks! e= " << pNow.t() << endl;
        cascade_list[wSplit].set_stat(-104);  // too soft for kicks stat
        return 1; 
      }
    }
    
    //cout << " q= " << q << endl;

    double bro_angle;
    double bro_k[3];
    AlignWithZ(pNow, bro_angle, bro_k);
    double phi = 2.*M_PI*dis(gen);
    double newpz = std::sqrt(newpzsq);
    pNow.Set(q*cos(phi),q*sin(phi),newpz,pNow.t());

    Rotation(pNow, -bro_angle, bro_k);
    
    cascade_list[wSplit].set_xf(xSplit);
    cascade_list[wSplit].reset_momentum(pNow);
  }

  return 1;
}

} // end namespace Adkoda

