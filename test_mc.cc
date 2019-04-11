#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <assert.h>
#include "Parton.h"
#include <fstream>
#include "global.h"

using namespace std;

double Sudakov(double t, double t0);
double generate_t2(double t0, double t1, double R, double numer);
double generate_x2(double t0, double t2, double R);

int alphas_constant;
int simple_splitting;

int main(int argc, char** argv) {

  assert(argc==4);

  alphas_constant=atoi(argv[2]);
  simple_splitting=atoi(argv[3]);

  int Nev=atoi(argv[1]);

  //Parameters  
  double t_max=1000.;
  double x_max=1.;
  double x_min=0.00001;
  double t0=2.;
  double Ein=100.;
  double Pxin=50.;
  double Pyin=0.;
  double Pzin=0.;
  double Pminusin=1./sqrt(2.)*(Ein-Pzin);
  double Pplusin=1./sqrt(2.)*(Ein+Pzin);

  //double biggest_angle=sqrt(t_max/Ein/Ein/zmin/(1.-zmin));

  int lund_nbins=100;
  double lund_hist[100][100]={{0.}};
  double t_norm[100]={0.};
  double t_binsize=log(t_max/2./t0)/double(lund_nbins);
  double z_binsize=log(1./x_min)/double(lund_nbins); 

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);

  //Parton Shower
  for (int iEv=1; iEv<=Nev; iEv++) {
 
    int nsplits=0;
 
    vector<Parton> parton_list;
    
    //Initial Hard Parton
    FourVector phard(Pxin, Pyin, Pzin, Ein);    
    FourVector x;
    Parton hard_parton( Parton(21,2,phard,x) );
    hard_parton.set_mom(-3);

    parton_list.push_back(hard_parton);

    int iter=0;
   
    bool changes=1;  
    while (true) {
      changes=0;
      //Scan the parton list to see whether there are partons that can split (status=1)
      int size_now=parton_list.size();
      for (int iP=0; iP<size_now; iP++)
      {
	//cout << " ip= " << iP << " status= " << parton_list[iP].stat() << endl; 
	//Parton parton=parton_list[iP];
        if (parton_list[iP].stat()<0) continue;
        changes=1;

	//If first iteration
        if (iter==0)
	{
	  t1=t_max;				//Max virtuality is that of hard scattering (Ein^2) 
          EvolveParton(parton_list[iP],1);	//Evolve the parton knowing its energy
	}
        else
        {
          int mom=parton_list[iP].mom();
          int sister=parton_list[mom].d1();
          if (sister==iP) sister=parton_list[mom].d2();
        
	  double ma2=parton_list[mom].virt();  
	  double mbmax=min(sqrt(ma2),parton_list[iP].p.t());
          double mcmax=min(sqrt(ma2),parton_list[sister].p.t());

	  double Eb, Ec, pb, pc, zb, zc;

	  bool evolve_b=1, evolve_c=1;
	  do {
            //Evolve this parton
	    if (evolve_b) EvolveParton(parton_list[iP],mbmax*mbmax,0);
          
	    //Evolve its sister
            if (evolve_c) EvolveParton(parton_list[sister],mcmax*mcmax,0);

	    evolve_b=0, evolve_c=0;

	    //Construct actual energies
            double mb2=parton_list[iP].virt();
            double mc2=parton_list[sister].virt();

	    double rb=(ma2+(mc2-mb2)-sqrt(pow(ma2-mb2-mc2,2.)-4.*mb2*mc2))/2./ma2;
	    double rc=(ma2-(mc2-mb2)-sqrt(pow(ma2-mb2-mc2,2.)-4.*mb2*mc2))/2./ma2;
        
            Eb=Eb0+(rc*Ec0-rb*Eb0);
            Ec=Ec0-(rc*Ec0-rb*Eb0);

	    pb=sqrt(Eb*Eb-mb2);
            pc=sqrt(Ec*Ec-mc2);

	    //Check z constraints
            bool zb_ok=1;
            bool zc_ok=1;

            zb=parton_list[iP].z();
            double zb_min=1./2.*(1.-pb/Eb);
            double zb_max=1./2.*(1.+pb/Eb);
	    if (zb<zb_min || zb>zb_max) zb_ok=0; 

            zc=parton_list[sister].z();
            double zc_min=1./2.*(1.-pc/Ec);
            double zc_max=1./2.*(1.+pc/Ec);
	    if (zc<zc_min || zc>zc_max) zc_ok=0;

	    //Check m constraints
	    bool mb_ok=1;
            bool mc_ok=1;
	    if (sqrt(mb2)+sqrt(mc2)>sqrt(ma2)) {
              double ratb=sqrt(mb2)/mbmax;
	      double ratc=sqrt(mc2)/mcmax;
	      if (ratb>ratc) mb_ok=0;
	      else mc_ok=0;
	    }

	    //Determine whether need to evolve further
            if (zb_ok==0 && zc_ok==1) evolve_b=1;
	    if (zb_ok==1 && zc_ok==0) evolve_c=1;
	    if (zb_ok==0 && zc_ok==0) {
              double ap_b=SOMETHING;
	      double ap_c=SOMETHING_ELSE;
	      if (ap_b>ap_c) evolve_b=1;
	      else evolve_c=1; 
	    }
	    if (zb_ok==1 && zc_ok==1) {
              if (mb_ok==0) evolve_b=1;
	      if (mc_ok==0) evolve_c=1;
	    }

	    //Update maximum virtualities
	    mbmax=sqrt(mb2);
	    mcmax=sqrt(mc2);

            //Make sure we don't evolve frozen partons
	    if (parton_list[iP].status()<0) evolve_b=0;
	    if (parton_list[sister].status()<0) evolve_c=0;

	  } while (evolve_b==1 || evolve_c==1);

	  //Update energy of splitters
	  parton_list[iP].set_en(Eb);
	  parton_list[sister].set_en(Ec); 
	   
          //Introduce daughters of b
	  if (parton_list[iP].status()>0) {
            FourVector p1(0.,0.,0.,zb*Eb);
	    FourVector p2(0.,0.,0.,(1.-zb)*Eb);
	    FourVector x1;
	    FourVector x2;

	    Parton d1 = Parton(21,1,p1,x1);
	    Parton d2 = Parton(21,1,p2,x2);
	  
            //Follow the heir
	    if (parton_list[iP].stat()==2) {
              if (zb>0.5) {
	        d1.set_stat(2);
	      }
	      else { 	
	        d2.set_stat(2);
              }
            }
	
	    d1.set_mom(iP);
	    d2.set_mom(iP);

	    parton_list.push_back(d1);         
	    parton_list[iP].set_d1(parton_list.size()-1); 
          
	    parton_list.push_back(d2);         
	    parton_list[iP].set_d2(parton_list.size()-1);
	  }

	  //Introduce daughters of c 
	  if (parton_list[sister].status()>0) {
            FourVector p1(0.,0.,0.,zc*Ec);
	    FourVector p2(0.,0.,0.,(1.-zc)*Ec);
	    FourVector x1;
	    FourVector x2;

	    Parton d1 = Parton(21,1,p1,x1);
	    Parton d2 = Parton(21,1,p2,x2);
	  
            //Follow the heir
	    if (parton_list[sister].stat()==2) {
              if (zc>0.5) {
	        d1.set_stat(2);
	      }
	      else { 	
	        d2.set_stat(2);
              }
            }
	
	    d1.set_mom(sister);
	    d2.set_mom(sister);

	    parton_list.push_back(d1);         
	    parton_list[sister].set_d1(parton_list.size()-1); 
          
	    parton_list.push_back(d2);         
	    parton_list[sister].set_d2(parton_list.size()-1);
	  }

	  //Deactivate splitted partons
	  parton_list[iP].set_status(-1);
          parton_list[sister].set_status(-1); 
	  
        } //End if iter!=0 
      
      } //End parton_list loop
        
      iter++;

      if (changes==0) break;

    } //Infinite changes loop
    
    //cout << " Parton List Size= " << parton_list.size() << endl;
    //cout << " Nsplits= " << nsplits << endl;
    
    parton_list.clear();
    if (iEv % 1 == 0) {
      cout << "#EVENT= " << iEv << endl;
      
      //Print Hist
      std::ostringstream flund;
      flund << "primary_lund_alphasconstant_" << alphas_constant << "simplesplitting_" << simple_splitting << ".dat";
      std::ofstream lundfile(flund.str().c_str(),std::ios_base::binary);
  
      for (unsigned int a=0; a<lund_nbins; a++) {
        double t_bincen = double(a)*t_binsize+t_binsize/2.;
        double z_spread = log(1.-exp(-(t_bincen+log(2.)))) + (t_bincen+log(2.));
        //cout << " z_spread= " << z_spread << " t_bincen= " << t_bincen << endl;
        for (unsigned int b=0; b<lund_nbins; b++) {
          double normal_bin_area = t_binsize * z_binsize;
          lundfile << double(a)*t_binsize+t_binsize/2. << " " << double(b)*z_binsize+z_binsize/2. << " " << lund_hist[a][b]/double(iEv) << " " << lund_hist[a][b]/double(iEv)/normal_bin_area << endl;
        }
        lundfile << endl;
      }
      lundfile.close();
      
    }
 
  //End Event Loop
  }

  //Print Hist
  std::ostringstream flund;
  flund << "primary_lund_alphasconstant_" << alphas_constant << "simplesplitting_" << simple_splitting << ".dat";
  std::ofstream lundfile(flund.str().c_str(),std::ios_base::binary);
  
  double maxspread=6.14344;
  for (unsigned int a=0; a<lund_nbins; a++) {
    double t_bincen = double(a)*t_binsize+t_binsize/2.;
    double z_spread = log(1.-exp(-(t_bincen+log(2.)))) + (t_bincen+log(2.));
    //cout << " z_spread= " << z_spread << " t_bincen= " << t_bincen << endl;
    for (unsigned int b=0; b<lund_nbins; b++) {
      double bin_area=t_binsize * z_binsize * maxspread/z_spread;
      double normal_bin_area = t_binsize * z_binsize;
      lundfile << double(a)*t_binsize+t_binsize/2. << " " << double(b)*z_binsize+z_binsize/2. << " " << lund_hist[a][b]/double(Nev) << " " << lund_hist[a][b]/double(Nev)/normal_bin_area << " " << lund_hist[a][b]/double(Nev)/bin_area << endl;
    }
      lundfile << endl;
  }
  lundfile.close();
  
  return 0;
}
