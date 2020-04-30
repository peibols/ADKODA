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
  alphas_med = 0.3;
  L_med = 4./0.1973; // in GeV^-1
  qhat= 1.5 * 0.1973; // in GeV^3
  eps_med = 0.01;
  med_cutoff = 1.25; // in GeV

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

  for (unsigned int i=0; i<ini_list.size(); i++) {

    Parton p = ini_list[i];
    
    if (p.stat()<0) continue;

    if (abs(p.id())<=6) continue; // Ignore quarks for now

    int last_part = parton_list.size()-1;

    // Align with z
    FourVector pIni = p.p();
    AlignWithZ(pIni, angle, k);
    
    std::vector<Parton> cascade_list;
    cascade_list.push_back(p);
    cascade_list[cascade_list.size()-1].reset_momentum(pIni);

    pplus = pIni.t();

    //std::cout << "Ini E= " << pIni.t() << " id= " << p.id() << std::endl;

    for (unsigned int j=0; j<cascade_list.size(); j++) {

      Parton q = cascade_list[j];
 
      if (q.stat()<0) continue;

      double curr_pos = q.x().p3abs() / 0.1973; // in GeV^-1
      //std::cout << " curr_pos= " << curr_pos << " L_med= " << L_med << std::endl;
      if (curr_pos > L_med) {
        std::cout << "Out of medium. " << std::endl;
        continue; // Static QGP Sphere
      }	

      double curr_time = tau(q.x().t());
      double start_time = curr_time;
      //std::cout << " Start tau= " << start_time << " t_kin= " << q.x().t() << std::endl;
      bool do_evolve=1;
      while (do_evolve) do_evolve = evolve(start_time, q, cascade_list, j);
    
    }

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
  return alphas_med * std::sqrt(qhat/pplus) * t / 0.1973; // dim less
}

double Cascade::tkin(double tau /*dim less*/)
{
  return tau / alphas_med * std::sqrt(pplus/qhat) * 0.1973; // in fm
}


} // end namespace Adkoda
