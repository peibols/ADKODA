#include <iostream>

#include "Cascade.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

Cascade::Cascade(const InitData &DATA_in) : DATA(DATA_in) {
  // Introduce kernels
  //FIXME shouldnt this sum depend on the max scale?
  Kgg* kgg = new Kgg (21, 21, 21);
  kernels.push_back(kgg);
/*
  Kgq* kgq;
  Kqq* kqq;
  for (int fl = 1; fl <= 5; fl++) {
    kqq = new Kqq (fl, fl, 21);
    kernels.push_back(kqq);
    kqq = new Kqq (-fl, -fl, 21);
    kernels.push_back(kqq);
    kgq = new Kgq (21, fl, -fl);
    kernels.push_back(kgq);
  }
*/

  // Medium parameters
  alphas_med = DATA.alphas_med;
  L0 = DATA.L0; // in GeV^-1
  qhat0= DATA.qhat0; // in GeV^3
  eps_med = DATA.eps_med;
  xmin_med = DATA.xmin_med;

  qmin = 0.2;

  //std::cout << "Cascade CONSTRUCTED" << std::endl;

}

void Cascade::init (std::vector<Parton> ini_list) {

  parton_list  = ini_list;
  
  return;

}

void Cascade::run () {

  // Need to move to Pz frame of initiator
  // Do whole cascade there
  // Apply broadening afterwards

  std::cout << "Initial Cascade List size = " << parton_list.size() << endl;
  std::cout << "Cascade RUNNING" << std::endl;

  std::vector<Parton> ini_list = parton_list;

  double initial_p[4]={0.};
  for (unsigned int i=0; i<ini_list.size(); i++) {
    Parton p = ini_list[i];
    //cout << i << " "; p.display();
    if (p.stat()<0) continue;
    initial_p[0]+=p.px();
    initial_p[1]+=p.py();
    initial_p[2]+=p.pz();
    initial_p[3]+=p.e();
  }
  cout << " Initial p= " << initial_p[0] << " " << initial_p[1] << " " << initial_p[2] << " " << initial_p[3] << endl;  

  for (unsigned int i=0; i<ini_list.size(); i++) {

    Parton p = ini_list[i];
    
    if (p.stat()<0) continue;

    if (abs(p.id())<=6) continue; // Ignore quarks for now

    //std::cout << "Ini E= " << pIni.t() << " id= " << p.id() << std::endl;
    
    double curr_pos = p.x().p3abs() / 0.1973; // in GeV^-1
    //std::cout << " curr_pos= " << curr_pos << " L_med= " << L_med << std::endl;
    if (curr_pos > L0) {
      std::cout << "Out of medium. " << std::endl;
      continue; // Static QGP Sphere
    }	
    
    int last_part = parton_list.size()-1;

    // Align with z
    FourVector pIni = p.p();
    AlignWithZ(pIni, angle, k);
    
    std::vector<Parton> cascade_list;
    cascade_list.push_back(p);
    cascade_list[cascade_list.size()-1].reset_momentum(pIni);

    pplus = pIni.t();
    xmin_med = DATA.T0 / pplus;
    if (xmin_med < eps_med )
    {
      cout << " PROBLEM xmin_med = " << xmin_med << " eps_med= " << eps_med << std::endl;
      cout << std::endl;
    }

    std::vector<Parton> active_list = cascade_list;
    std::vector<int> active_map;
    active_map.push_back(cascade_list.size()-1);

    double curr_time = tau(p.x().t());
    double start_time = curr_time;
    //std::cout << " Start tau= " << start_time << " t_kin= " << q.x().t() << std::endl;
    bool do_evolve=1;
    while (do_evolve) do_evolve = evolve(start_time, cascade_list, active_list, active_map);

    // Update full parton list
    //std::cout << " i= " << i << " last_part= " << last_part << std::endl;
    for (unsigned int j=0; j<cascade_list.size(); j++) {
      
      // Rotate Back
      FourVector ptemp = cascade_list[j].p();
      Rotation(ptemp, -angle, k);
      cascade_list[j].reset_momentum(ptemp);

      //std::cout << " j= " << j << " mom = " << cascade_list[j].mom1() << std::endl;

      // Update daughters if it has any
      if (cascade_list[j].d1()!=0) {
        cascade_list[j].set_d1(cascade_list[j].d1() + last_part);
        cascade_list[j].set_d2(cascade_list[j].d2() + last_part);
      }

      // Update mothers
      if (j!=0) {
        if (cascade_list[j].mom1()==0) cascade_list[j].set_mom1(i); // Maps to initial list
        else cascade_list[j].set_mom1(cascade_list[j].mom1() + last_part); // Maps to total list
      }   

      if (j==0) parton_list[i] = cascade_list[0]; // Replace initiator
      else {
        parton_list.push_back(cascade_list[j]);
      }
    }

    //break;

  }

  std::cout << "Cascade FINISHED" << std::endl;
  std::cout << "Final Cascade List size = " << parton_list.size() << endl;

}

double Cascade::tau(double t /*in fm*/)
{
  return alphas_med * std::sqrt(qhat0/pplus) * t / 0.1973; // dim less
}

double Cascade::tkin(double tau /*dim less*/)
{
  return tau / alphas_med * std::sqrt(pplus/qhat0) * 0.1973; // in fm
}


} // end namespace Adkoda
