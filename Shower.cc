#include <iostream>

#include "Shower.h"
#include "Util.h"

using namespace Util;

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
    
 /*   
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
*/

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

  in_third = 0;

  //std::cout << " Weight= " << event_weight << " Xsec= " << event_xsec << " Pt Max= " << hard_pt_max << std::endl;

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
    t_max = std::pow(hard_pt_max, 2.) * 16.;
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

void Shower::third_stage_init() {

  //std::cout << "Initializing 3rd Stage" << endl;
      
  in_third = 1;
      
  double Qs2 = DATA.qhat0 * DATA.L0;
  //std::cout << "Qs2= " << Qs2 << std::endl;
  t_min = std::pow(DATA.pt_min, 2.);
  //t_max = Qs2;
  //t_max = std::pow(hard_pt_max,2.);
  t_max = stop_scale;
  //t_max = std::pow(hard_pt_max, 2.) * 16.;

  //for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    //parton_list[ip].set_is_frozen(0);
  //}
  
  /*
  double lambda=dis(gen);

  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
     FourVector p = parton_list[ip].p();
     FourVector pnew;
     pnew.Set(p.x()*lambda,p.y()*lambda,p.z()*lambda,p.t()*lambda);
     parton_list[ip].reset_momentum(pnew);
  }
  */
/*
  // Find hardest, count active
  double ptmax= 0;
  int imax = -1000;
  int nactive = 0;
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    Parton p = parton_list[ip];
    if (p.stat()<0) continue;
    if (p.stat()==63) {
      parton_list[ip].set_stat(64);
      continue;
    }
    nactive++;
    if (p.pt() > ptmax) ptmax = p.pt(), imax=ip;
  }
  std::vector<bool> aredone (int(parton_list.size()),0);
  
  max_colour = 101;
  int inow = imax;
  int col[2] = {101,102};
  max_colour++;
  parton_list[inow].set_cols(col);
  aredone[inow]=1;
  nactive--;
  while (true) {

    Parton p = parton_list[inow];
    double min_delR = 1000000;
    int iclose = -1000;
    for (unsigned int ip=0; ip<parton_list.size(); ip++) {
      Parton k = parton_list[ip];
      if (aredone[ip] == 1) continue;
      if (k.stat()<0 || k.stat()==64) continue;
      double delR = Util::Delta(p.p(),k.p());
      //double phi1 = p.phi();
      //double phi2 = k.phi();
      //double deltaphi = fabs(phi1-phi2);
      //std::cout << " phi1= " << phi1 << " phi2= " << phi2 << " deltaphi= " << deltaphi << std::endl;
      if (delR < min_delR) min_delR = delR, iclose = ip;
    }
   
    if (iclose == -1000) {
      std::cout << " Could not find close!" << std::endl;
      exit(1);
    }

    // Make link
    int pcol[2] = {max_colour,max_colour+1};
    max_colour++;
    parton_list[iclose].set_cols(pcol);
    
    aredone[iclose]=1;
    nactive--;

    if (nactive==0) break;
    inow = iclose;

  }
*/
/*
  std::vector<Parton> thermal_list;
  int therm_id = 21;

  max_colour = 101;

  for (unsigned int ip=0; ip<parton_list.size(); ip++) {

    Parton p = parton_list[ip];
    if (p.stat()<0) continue;
    if (p.stat()==63) {
      parton_list[ip].set_stat(64); // Remnant frozen in third stage
      continue;
    }

    double En = p.p().t();
    double Th=0.00001;
    double T = DATA.T0;
    while (true) {
      Th = 8.*T*dis(gen);
      double func = Th * exp(-Th/T) / (T/exp(1));
      double rand = dis(gen);
      if (func>rand) break;
    }
    
    double costheta = cos(2.*M_PI*dis(gen));
    double kz = Th * costheta;
    double kT = Th * sqrt(1-costheta*costheta);
    double phi = 2.*M_PI*dis(gen);
    double kx = kT * cos(phi);
    double ky = kT * sin(phi);
    FourVector pTherm;
    pTherm.Set(kx,ky,kz,Th);
    //double angle = -1000;
    //double k[3] = {0.};
    //FourVector pShower = p.p();
    //AlignWithZ(pShower, angle, k);
    //Rotation(pTherm, -angle, k);


    FourVector xa;
    Parton thermal(Parton(therm_id, 80, pTherm, xa));
    thermal.set_mom1(0), thermal.set_mom2(0);

    int col1[2]={max_colour, max_colour+1};
    int col2[2]={max_colour+1, max_colour};

    parton_list[ip].set_cols(col1);
    thermal.set_cols(col2);
    
    max_colour+=2;

    thermal_list.push_back(thermal);
  
  }

  for (unsigned int ip=0; ip<thermal_list.size(); ip++) {
    parton_list.push_back(thermal_list[ip]);
  }

  max_colour-=1; //last color actually

*/

/*
  std::vector<Parton> new_list;
  max_colour = 101;

  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    Parton p = parton_list[ip];
    if (p.stat()<0) continue;
    if (p.stat()==63) continue;
    Parton pnew = p;
    int col[2]={max_colour,max_colour+1};
    pnew.set_cols(col);
    new_list.push_back(pnew);
    int npartners = 0;
    bool corrected = 0;
    for (unsigned int jp=0; jp<parton_list.size(); jp++) {
      if (ip==jp) continue;
      Parton k = parton_list[jp];
      if (k.stat()<0) continue;
      if (p.ColourConnected(k)) {
        double dipmass = m2(p.p(),k.p());
	double Qsvirt = std::sqrt(Qs2) * p.p().t();
	//Qsvirt=0.;
	Qsvirt*=2;
	//double thetacscale = p.stop_scale();
	//double minvirt = std::max(Qsvirt,thetacscale);
	double minvirt = Qsvirt;
	//std::cout << " dip mass= " << dipmass << " Qsvirt= " << Qsvirt << " ThetacScale= " << thetacscale << " ei= " << p.p().t() << " ej= " << k.p().t() << std::endl;
     
     	int kcol[2];
	if (npartners==0) {
          kcol[0]=max_colour+1;
	  kcol[1]=max_colour+10000;
	  npartners++;
	}
	else {
          kcol[0]=max_colour+10001;
	  kcol[1]=max_colour;
	}
     
     	if (dipmass >= minvirt || corrected == 1) {
          k.set_cols(kcol);
          k.set_stat(80);
          new_list.push_back(k);
	}
*/
/*	
        else {
          double factor = minvirt/dipmass;
	  FourVector knew;
	  knew.Set(k.p().x()*factor,k.p().y()*factor, k.p().z()*factor, k.p().t()*factor);
	  double checkmass = m2(knew,p.p());
	  //std::cout << " checkmass = " << checkmass << " minvirt= " << minvirt << " factor= " << factor << std::endl;
	  k.reset_momentum(knew);
          k.set_cols(kcol);
          k.set_stat(80);
	  new_list.push_back(k);
	  corrected = 1;
	}
*/
	/*
	else {
          double En=p.p().t();
          double Th=0.00001;
          double T = DATA.T0;
	  int ntry=0;
	  bool couldnot=0;
          while (true) {
	    if (ntry>30) {
	      couldnot=1;
	      break;
	    }
	    while (true) {
              Th = 8.*T*dis(gen);
              double func = Th * exp(-Th/T) / (T/exp(1));
              double rand = dis(gen);
              if (func>rand) break;
            }
	    if (minvirt/2/En/Th<2) break;
            else ntry++;
	  }
	  if (couldnot) {
            std::cout << "Could not find!" << std::endl;
	  }
	  else {
            double costheta = 1-minvirt/2/En/Th;
            double kz = Th * costheta;
            double kT = Th * sqrt(1-costheta*costheta);
            double phi = 2.*M_PI*dis(gen);
            double kx = kT * cos(phi);
            double ky = kT * sin(phi);
            FourVector pTherm;
            pTherm.Set(kx,ky,kz,Th);
            double angle = -1000;
            double k[3] = {0.};
            FourVector pShower = p.p();
            AlignWithZ(pShower, angle, k);
            Rotation(pTherm, -angle, k);

            FourVector xa;
            Parton thermal(Parton(21, 80, pTherm, xa));
            thermal.set_mom1(0), thermal.set_mom2(0);

	    thermal.set_cols(kcol);
	    new_list.push_back(thermal);

	   // std::cout << " Fixed virt, now= " << m2(p.p(),thermal.p()) << std::endl;
	  } 

	}

      }
    }
    max_colour+=2;
  }

  max_colour-=1;
*/
/*
  for (unsigned int ip=0; ip<new_list.size(); ip++) {
    Parton p = new_list[ip];
    std::cout << "en= " << p.p().t() << " stat= " << p.stat() << " col= " << p.col() << " acol= " << p.acol() << std::endl;
  }
*/

  //parton_list = new_list;
  
  return;
}

void Shower::run () {

  //std::cout << "Initial Parton List size = " << parton_list.size() << endl;
  //print();
  //std::cout << "Shower RUNNING" << std::endl;

  // Vacuum Shower
  bool do_evolve = 1;
  stop_scale = DATA.pt_min*DATA.pt_min;
  first_stop = 0;
  while (do_evolve) do_evolve = evolve(); // FIXME how does this evolve() called?!

  if (!in_third) std::cout << " Stop Scale= " << stop_scale << std::endl;

  //std::cout << "Shower FINISHED" << std::endl;
  //std::cout << "Final Parton List size = " << parton_list.size() << endl;
  //print();

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
  std::cout << "#ip\t ID\t Stat\t m1\t m2\t d1\t d2\t c\t ac\t px\t py\t pz\t E\t\t m\t"
  //<< "pt2\t\t"
  << "x\t y\t z\t t"
  << std::endl;
  for (unsigned int ip=0; ip<parton_list.size(); ip++) {
    std::cout << ip << "\t ";
    parton_list[ip].display();
  }
  return;
}

} // end namespace Adkoda
