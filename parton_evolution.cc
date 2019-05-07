#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include "global.h"
#include "Parton.h"
 
using namespace std;
 
std::random_device rd;
std::mt19937 gen(1);
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

void EvolveParton(Parton& parton, double tmax)
{

  double m_min=sqrt(t0);
  double cutoff=1./2.*(1.-sqrt(1.-4.*m_min*m_min/tmax));

  int counter=0;
  //cout << " cutoff= " << cutoff << " E= " << E << " m_min= " << m_min << " tmax= " << tmax << endl;
  double t_temp=tmax;
  double z_temp=1.; 
  while (true) {
   
    if (sqrt(t_temp)<2.*m_min) break;
    cutoff=1./2.*(1.-sqrt(1.-4.*m_min*m_min/tmax));

    //Pick t
    double Rt=dis(gen);
    //cout << " prim_envel= " << prim_envel(tmax,cutoff) << endl;
    t_temp=inv_prim_envel(prim_envel(tmax,cutoff)+log(Rt),cutoff);
    if (sqrt(t_temp)<2.*m_min) {
      //cout << " freezing with t_temp= " << t_temp << endl;
      break;
    }

    //Pick z
    double actual_cutoff=1./2.*(1.-sqrt(1.-4.*m_min*m_min/t_temp));
    double Rz=dis(gen);
    z_temp=inv_zprim_envel(zprim_envel(cutoff)+Rz*(zprim_envel(1.-cutoff)-zprim_envel(cutoff)));     
    //Evaluate ratio
    double Rhm=dis(gen);
    double ratio;
    if (z_temp<1.-actual_cutoff && z_temp>actual_cutoff) ratio=integrand(t_temp,z_temp)/envel(t_temp,z_temp); 
    else ratio=0.;
    if (ratio<0.) { cout << " RATIO NEGA IN HARD PARTON, ratio= " << ratio << endl; exit(0); }

    //cout << " t_temp= " << t_temp << " z_temp= " << z_temp << endl;
    //cout << " ratio= " << ratio << " Rhm= " << Rhm << endl;
    if (ratio>Rhm) break; 
    counter++;
    tmax=t_temp;
    //if (counter>50) exit(0);
  }

  //Put on-mass-shell if below m_min
  //cout << " t_temp after= " << t_temp << endl << endl;
  if (sqrt(t_temp)<2.*m_min) t_temp=0., z_temp=1., parton.set_stat(-1);

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


void SetKinematicsSudakovBasis(Parton pmom, Parton &pb, Parton &pc, double prog[4])
{
  double mprog2=prog[3]*prog[3]-prog[0]*prog[0]-prog[1]*prog[1]-prog[2]*prog[2];
  double modn=sqrt(prog[0]*prog[0]+prog[1]*prog[1]+prog[2]*prog[2]);
  double n[4]={-prog[0]/modn,-prog[1]/modn,-prog[2]/modn,1.};
  double ndotp=n[3]*prog[3]-n[0]*prog[0]-n[1]*prog[1]-n[2]*prog[2];   

  double mom_z=pmom.z();
  double alpha_mom=pmom.alpha();
  double ma2=pmom.virt();
  double mb2=pb.virt();
  double mc2=pc.virt();

  double alpha_b=alpha_mom*mom_z;
  double alpha_c=alpha_mom-alpha_b;
  
  double kptwo=mom_z*(1.-mom_z)*ma2-(1.-mom_z)*mb2-mom_z*mc2;
  double mom_px=pmom.p().x();
  double mom_py=pmom.p().y();

  double pcx, pcy, pbx, pby; 
  double b_pperp2, c_pperp2; 
    
  double rand_ang=2.*M_PI*dis(gen);
  double kperp_x=sqrt(kptwo)*cos(rand_ang);
  double kperp_y=sqrt(kptwo)*sin(rand_ang);
  pbx=kperp_x+mom_z*mom_px;
  pby=kperp_y+mom_z*mom_py;
  pcx=mom_px-pbx;
  pcy=mom_py-pby; 
  c_pperp2=pcx*pcx+pcy*pcy;
  b_pperp2=pbx*pbx+pby*pby;

  double beta_b=(mb2-alpha_b*alpha_b*mprog2-b_pperp2)/(2.*alpha_b*ndotp);
  double beta_c=(mc2-alpha_c*alpha_c*mprog2-c_pperp2)/(2.*alpha_c*ndotp);

  double b_pz=alpha_b*prog[2]+beta_b*n[2];  
  double c_pz=alpha_c*prog[2]+beta_c*n[2];  

  double b_en=alpha_b*prog[3]+beta_b*n[3];  
  double c_en=alpha_c*prog[3]+beta_c*n[3];  
  
  FourVector p_b(pbx,pby,b_pz,b_en);
  FourVector p_c(pcx,pcy,c_pz,c_en);
  pb.set_p(p_b);
  pc.set_p(p_c);
  pb.set_alpha(alpha_b);
  pc.set_alpha(alpha_b);
}

//Set Kinematics after z,ta,tb,tc are known
void SetKinematics(Parton pmom, Parton &pb, Parton &pc)
{
  double mom_z=pmom.z();
  double ma2=pmom.virt();
  double mb2=pb.virt();
  double mc2=pc.virt();

  double kptwo=mom_z*(1.-mom_z)*ma2-(1.-mom_z)*mb2-mom_z*mc2;
  double mom_px=pmom.p().x();
  double mom_py=pmom.p().y();

  double b_pp=pmom.pplus()*mom_z;
  double c_pp=pmom.pplus()*(1.-mom_z);

  double pcx, pcy, pbx, pby; 
  double b_pperp2, c_pperp2; 
    
  double rand_ang=2.*M_PI*dis(gen);
  double kperp_x=sqrt(kptwo)*cos(rand_ang);
  double kperp_y=sqrt(kptwo)*sin(rand_ang);
  pbx=kperp_x+mom_z*mom_px;
  pby=kperp_y+mom_z*mom_py;
  pcx=mom_px-pbx;
  pcy=mom_py-pby; 
  c_pperp2=pcx*pcx+pcy*pcy;
  b_pperp2=pbx*pbx+pby*pby;

  double b_pm=(mb2+b_pperp2)/2./b_pp;
  double c_pm=(mc2+c_pperp2)/2./c_pp;

  double b_en=1./sqrt(2.)*(b_pp+b_pm);
  double b_pz=1./sqrt(2.)*(b_pp-b_pm);

  double c_en=1./sqrt(2.)*(c_pp+c_pm);
  double c_pz=1./sqrt(2.)*(c_pp-c_pm);

  FourVector p_b(pbx,pby,b_pz,b_en);
  FourVector p_c(pcx,pcy,c_pz,c_en);
  pb.set_p(p_b);
  pc.set_p(p_c);
}

void AddDaughters(vector<Parton> &parton_list, int iP)
{
	  
  //cout << " adding daughter " << endl; 
  FourVector p1,p2,x1,x2;
  Parton d1 = Parton(21,1,p1,x1);
  Parton d2 = Parton(21,1,p2,x2);
	  
  //Follow the heir
  double z=parton_list[iP].z();
  if (parton_list[iP].stat()==2) {
    if (z>0.5) {
      d1.set_stat(2);
    }
    else { 	
      d2.set_stat(2);
    }
  }
	
  d1.set_mom(iP);
  d2.set_mom(iP);
  //cout << " iP= " << iP << endl;

  parton_list.push_back(d1);         
  parton_list[iP].set_d1(parton_list.size()-1); 
          
  parton_list.push_back(d2);         
  parton_list[iP].set_d2(parton_list.size()-1);

  parton_list[iP].set_stat(-1);

}

void rotate_prog(double prog[4]){

  //Normalise prog vector
  double mod_prog=sqrt(prog[0]*prog[0]+prog[1]*prog[1]+prog[2]*prog[2]);
  double u_prog[3]={prog[0]/mod_prog,prog[1]/mod_prog,prog[2]/mod_prog};
  
  //Define desired orientation of prog
  double u_z[3]={0.,0.,1.};

  //Compute rotation plane and normalise perp vector
  double axis[3]={u_prog[1]*u_z[2]-u_prog[2]*u_z[1],u_prog[2]*u_z[0]-u_prog[0]*u_z[2],u_prog[0]*u_z[1]-u_prog[1]*u_z[0]};
  double mod_axis=sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);
  for (unsigned a=0; a<3; a++) axis[a]=axis[a]/mod_axis;
  
  //Find rotation angle
  double sin_rot_angle=mod_axis;
  double cos_rot_angle=(u_z[0]*prog[0]+u_z[1]*prog[1]+u_z[2]*prog[2])/mod_prog;

  //Compute cross and dot products
  double kcrossv[3]={axis[1]*prog[2]-axis[2]*prog[1],axis[2]*prog[0]-axis[0]*prog[2],axis[0]*prog[1]-axis[1]*prog[0]};
  double kdotv=axis[0]*prog[0]+axis[1]*prog[1]+axis[2]*prog[2];
  
  //Perform rotation
  double rot_x=prog[0]*cos_rot_angle+kcrossv[0]*sin_rot_angle+axis[0]*kdotv*(1.-cos_rot_angle);
  double rot_y=prog[1]*cos_rot_angle+kcrossv[1]*sin_rot_angle+axis[1]*kdotv*(1.-cos_rot_angle);
  double rot_z=prog[2]*cos_rot_angle+kcrossv[2]*sin_rot_angle+axis[2]*kdotv*(1.-cos_rot_angle);

  double mod_rot=sqrt(rot_x*rot_x+rot_y*rot_y+rot_z*rot_z);
  cout << " axis[0]= " << axis[0] << " axis[1]= " << axis[1] << " axis[2]= " << axis[2] << endl;

  cout << " orig mod= " << mod_prog << " rot mod= " << mod_rot << endl;
  cout << " prog[0]= " << prog[0] << " rot[0]= " << rot_x << endl;
  cout << " prog[1]= " << prog[1] << " rot[1]= " << rot_y << endl;
  cout << " prog[2]= " << prog[2] << " rot[2]= " << rot_z << endl;

}
