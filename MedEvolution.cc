
#include <iostream>
#include <fstream>

#include "Cascade.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

bool Cascade::evolve(double &start_time, std::vector<Parton> &cascade_list, std::vector<Parton> &active_list, std::vector<int> &active_map)
{
  
  double t_ini = start_time;
  double tlarge=1000000000.;
  double t = tlarge;

  int wKernel;
  int wSplit;

  for (int iSplit=0; iSplit<active_list.size(); iSplit++) {
  
    Parton p = active_list[iSplit];
 
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

  //cout << " t= " << t << endl;
  if (t>1) {
    //cout << " after t= " << t << endl;
  }

  if (t>2) return 0;

  start_time = t; // Update time

  Parton p = active_list[wSplit];
  int map = active_map[wSplit]; 
 
  double time = tkin(t); // in fm
  //cout << " time= " << time << endl;

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
    cascade_list[map].set_stat(103); // Outside medium stat
    active_list.erase(active_list.begin()+wSplit);
    active_map.erase(active_map.begin()+wSplit);
    //std::cout << "Outside medium." << std::endl;
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

    FourVector pNow = p.p();
    double split_angle;
    double split_k[3];
    AlignWithZ(pNow, split_angle, split_k);

    FourVector pDau1, pDau2;
    pDau1.Set(pNow.x()*z,pNow.y()*z,pNow.z()*z,pNow.t()*z);
    pDau2.Set(pNow.x()*(1.-z),pNow.y()*(1.-z),pNow.z()*(1.-z),pNow.t()*(1.-z));

    Rotation(pDau1, -split_angle, split_k);
    Rotation(pDau2, -split_angle, split_k);
     
    int dau1_id=21;
    int dau2_id=21;

    int dau1_stat=101;
    int dau2_stat=101;
  
    double d1_xfrac = z * p.xfrac();
    double d2_xfrac = (1.-z) * p.xfrac();
    if (d1_xfrac < xmin_med) {
      dau1_stat=-102;
      //std::cout << " Thermalised d1 e= " << pDau1.t() << " at t= " << t << std::endl;
    }
    if (d2_xfrac < xmin_med) {
      dau2_stat=-102;
      //std::cout << " Thermalised d2 e= " << pDau2.t() << " at t= " << t << std::endl;
    }

    Parton daughter1(Parton(dau1_id, dau1_stat, pDau1, xSplit)); //Status: 101, active cascade particles
    Parton daughter2(Parton(dau2_id, dau2_stat, pDau2, xSplit));
  
    daughter1.set_mom1(map), daughter1.set_mom2(0);
    daughter2.set_mom1(map), daughter2.set_mom2(0);
  
    daughter1.set_xfrac(d1_xfrac);
    daughter2.set_xfrac(d2_xfrac);
  
    cascade_list.push_back(daughter1);
    int Dau1 = int(cascade_list.size()) - 1;
    if (dau1_stat>0) {
      active_list.push_back(daughter1);
      active_map.push_back(Dau1);
    }

    cascade_list.push_back(daughter2);
    int Dau2 = int(cascade_list.size()) - 1;
    if (dau2_stat>0) {
      active_list.push_back(daughter2);
      active_map.push_back(Dau2);
    }

    // Update splitter
    cascade_list[map].set_stat(-p.stat()); //It decayed: state -> decayed
    cascade_list[map].set_xf(xSplit);
    cascade_list[map].set_d1(Dau1);
    cascade_list[map].set_d2(Dau2);
    active_list.erase(active_list.begin()+wSplit);
    active_map.erase(active_map.begin()+wSplit);

  }
  else {  // Broadening

    double q;
    while (true) {
      q = GenerateQ(dis(gen));
      break;
    }
    
    //q=0.;

    FourVector pNow = p.p();
    
    double bro_angle;
    double bro_k[3];
    AlignWithZ(pNow, bro_angle, bro_k);
    double phi = 2.*M_PI*dis(gen);
    double t_pplus =1./2.*(pNow.t()+pNow.z());
    double newen = t_pplus + q*q/4./t_pplus; 
    double newpz = t_pplus - q*q/4./t_pplus;
    pNow.Set(q*cos(phi),q*sin(phi),newpz,newen);
    if (newen!=newen) { cout << " newen nan! " << endl; exit(1); }
    Rotation(pNow, -bro_angle, bro_k);
  
    int dau_stat = 101; 
    int dau_id = p.id();
    Parton daughter(Parton(dau_id, dau_stat, pNow, xSplit)); //Status: 101, active cascade particles
  
    daughter.set_mom1(wSplit), daughter.set_mom2(0);
  
    daughter.set_xfrac(p.xfrac());
  
    cascade_list.push_back(daughter);
    int Dau = int(cascade_list.size()) - 1;
    active_list.push_back(daughter);
    active_map.push_back(Dau);

    // Update splitter
    cascade_list[map].set_stat(-103); // kicked particle
    cascade_list[map].set_xf(xSplit);
    cascade_list[map].set_d1(Dau);
    cascade_list[map].set_d2(0);
    cascade_list[map].reset_momentum(pNow);
    if (dau_stat>0) {
      active_list.erase(active_list.begin()+wSplit);
      active_map.erase(active_map.begin()+wSplit);
    }

  }

  return 1;
}

} // end namespace Adkoda

