#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include "global.h"
#include "Parton.h"
  
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.,1.);

double envel(double t, double z);
double integrand(double t, double z);
double prim_envel(double t, double cutoff);
double inv_prim_envel(double y, double cutoff);
double zprim_envel(double z);
double inv_zprim_envel(double y);
double int_Pgg_kernel(double z);
double alpha_s(double t);
double P_gg(double z);

void EvolveSisters(Parton& parton_b, Parton& parton_c, double tmax_b, double t_maxc, bool evolve_b, bool evolve_c)
{
  double m_min=sqrt(2.);

  double Eb0=parton_b.p().t(); 
  double Ec0=parton_c.p().t();

  double cutoff_b, cutoff_c;
  if (tmax_b>m_min*m_min) cutoff_b= 1./2.*(1.-sqrt(Eb0*Eb0-m_min*m_min)/Eb0);
  if (tmax_c>m_min*m_min) cutoff_c= 1./2.*(1.-sqrt(Ec0*Ec0-m_min*m_min)/Ec0);
 
  double t_temp_b=tmax_b; 
  double t_temp_c=tmax_c;

  double z_temp_b=1.;
  double z_temp_c=1.;

  bool evolve_b=1;
  bool evolve_c=1;
  do {
   
    if (t_max_b<m_min*m_min) {
      evolve_b=0;
      t_temp_b=0.;
    }

    if (t_max_c<m_min*m_min) {
      evolve_c=0;
      t_temp_c=0.;
    }
     
    if (evolve_b==1) {
      //Pick tb
      double Rt_b=dis(gen);
      t_temp_b=inv_prim_envel(prim_envel(tmax_b,cutoff_b)+log(Rt_b),cutoff_b);
      if (t_temp_b<m_min*m_min) {
        cout << " freezing b with t_temp= " << t_temp_b << endl;
      }
    }
    
    if (evolve_c==1) {  
      //Pick tc
      double Rt_c=dis(gen);
      t_temp_c=inv_prim_envel(prim_envel(tmax_c,cutoff_c)+log(Rt_c),cutoff_c);
      if (t_temp_c<m_min*m_min) {
        cout << " freezing c with t_temp= " << t_temp_c << endl;
      }
    }

    //Compute actual energies
    double rb=(tmom+(t_temp_c-t_temp_b)-sqrt(pow(tmom-t_temp_b-t_temp_c,2.)-4.*t_temp_b*t_temp_c))/2./tmom;
    double rc=(tmom-(t_temp_c-t_temp_b)-sqrt(pow(tmom-t_temp_b-t_temp_c,2.)-4.*t_temp_b*t_temp_c))/2./tmom;
    
    Eb=Eb0+(rc*Ec0-rb*Eb0);
    Ec=Ec0-(rc*Ec0-rb*Eb0);
    pb=sqrt(Eb*Eb-t_temp_b);
    pc=sqrt(Ec*Ec-t_temp_c);

    //Update cutoffs
    cutoff_b=1./2.*(1.-sqrt(Eb*Eb-t_temp_b)/Eb);    
    cutoff_c=1./2.*(1.-sqrt(Ec*Ec-t_temp_c)/Ec);

    if (evolve_b==1) {
      //Pick zb
      double Rz_b=dis(gen);
      z_temp_b=inv_zprim_envel(zprim_envel(cutoff_b)+Rz_b*(zprim_envel(1.-cutoff_b)-zprim_envel(cutoff_b)));     
      //Evaluate ratio of b
      double Rhm_b=dis(gen);
      double ratio_b=integrand(t_temp_b,z_temp_b)/envel(t_temp_b,z_temp_b);
      if (ratio_b>Rhm) evolve_b=0; 
    }

    if (evolve_c==1) {
      //Pick zb
      double Rz_c=dis(gen);
      z_temp_c=inv_zprim_envel(zprim_envel(cutoff_c)+Rz_c*(zprim_envel(1.-cutoff_c)-zprim_envel(cutoff_c)));     
      //Evaluate ratio of b
      double Rhm_c=dis(gen);
      double ratio_c=integrand(t_temp_c,z_temp_c)/envel(t_temp_c,z_temp_c);
      if (ratio_c>Rhm_c) evolve_c=0; 
    }

    //Update maximum virtualities
    t_max_b=t_temp_b;
    t_max_c=t_temp_c;

  } while (evolve_b==1 || evolve_c==1);

  //Put on-mass-shell if below m_min
  if (t_temp_b<m_min*m_min) t_temp_b=0., z_temp_b=1., parton_b.set_stat(-1);
  parton_b.set_virt(t_temp_b);
  parton_b.set_z(z_temp_b);
  
  if (t_temp_c<m_min*m_min) t_temp_c=0., z_temp_c=1., parton_c.set_stat(-1);
  parton_c.set_virt(t_temp_c);
  parton_c.set_z(z_temp_c);
}

void EvolveParton(Parton& parton, double tmax, bool iter)
{

  double E=parton.p().t();
  double m_min=sqrt(2.);
  double cutoff=1./2.*(1.-sqrt(E*E-m_min*m_min)/E); 

  int counter=0;
  cout << " cutoff= " << cutoff << " E= " << E << " m_min= " << m_min << " tmax= " << tmax << endl;
  double t_temp=tmax;
  double z_temp=1.; 
  while (true) {
   
    if (t_temp<m_min*m_min) break;
 
    //Pick t
    double Rt=dis(gen);
    cout << " prim_envel= " << prim_envel(tmax,cutoff) << endl;
    t_temp=inv_prim_envel(prim_envel(tmax,cutoff)+log(Rt),cutoff);
    if (t_temp<m_min*m_min) {
      cout << " freezing with t_temp= " << t_temp << endl;
      break;
    }

    //Pick z
    //Update cutoff
    cutoff=1./2.*(1.-sqrt(E*E-t_temp)/E);
    cout << " updated cutoff= " << cutoff << endl;
    double Rz=dis(gen);
    z_temp=inv_zprim_envel(zprim_envel(cutoff)+Rz*(zprim_envel(1.-cutoff)-zprim_envel(cutoff)));     

    //Evaluate ratio
    double Rhm=dis(gen);
    double ratio=integrand(t_temp,z_temp)/envel(t_temp,z_temp); 

    cout << " t_temp= " << t_temp << " z_temp= " << z_temp << endl;
    cout << " ratio= " << ratio << " Rhm= " << Rhm << endl;
    if (ratio>Rhm) break;
    counter++;
    tmax=t_temp;
    if (counter>50) exit(0);
  }

  //Put on-mass-shell if below m_min
  cout << " t_temp after= " << t_temp << endl << endl;
  if (t_temp<m_min*m_min) t_temp=0., z_temp=1., parton.set_stat(-1);

  parton.set_virt(t_temp);
  parton.set_z(z_temp);
   
}

//g(t,z)
double envel(double t, double z)
{
  double CA=3.;
  double kernel=CA*(1./z+1./(1.-z));
  double max_alphas=0.5;
  return kernel * max_alphas/2./M_PI / t;
}

//f(t,z)
double integrand(double t, double z)
{
  double kernel=P_gg(z);
  return kernel * alpha_s(t*z*(1.-z))/2./M_PI / t;
}

//G(t)
double prim_envel(double t, double cutoff)
{
  double CA=3.;
  double max_alphas=0.5;
  double int_kernel=CA*4.*atanh(1.-2.*cutoff);
  double result= int_kernel * max_alphas/2./M_PI * log(t);
  return result; 
}

//G^-1
double inv_prim_envel(double y, double cutoff)
{
  double CA=3.;
  double max_alphas=0.5;
  double int_kernel=CA*4.*atanh(1.-2.*cutoff);
  double result= exp( y / (int_kernel * max_alphas/2./M_PI) ); 
  return result;
}

//H(z)
double zprim_envel(double z)
{
  return log(z)-log(1.-z);
}

//H^-1
double inv_zprim_envel(double y)
{
  return exp(y)/(exp(y)+1.);     
}

//Integrated Pgg kernel (no alphas)
double int_Pgg_kernel(double z)
{
  double CA=3.;
  return CA*(-1./6.*z*(z*(2.*z-3.)+12.)-2.*atanh(1.-2.*z));
}
