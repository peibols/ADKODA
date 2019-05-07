#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <assert.h>
#include "Parton.h"
#include <fstream>
#include "global.h"

using namespace std;

int alphas_constant;
int simple_splitting;
double t0;

double int_Pgg_kernel(double z);
void EvolveParton(Parton& parton, double tmax);
void SetKinematics(Parton pmom, Parton& pb, Parton& pc);
void SetKinematicsSudakovBasis(Parton pmom, Parton& pb, Parton& pc, double prog[4]);
void rotate_prog(double prog[4]);
void AddDaughters(vector<Parton> &parton_list, int iP);

int main(int argc, char** argv) {

  assert(argc==4);

  alphas_constant=atoi(argv[2]);
  simple_splitting=atoi(argv[3]);

  int Nev=atoi(argv[1]);

  //Parameters  
  t0=1.;
  double Ein=1000.;
  double Pxin=0.;
  double Pyin=0.;
  double Pzin=0.;
  double Pminusin=1./sqrt(2.)*(Ein-Pzin);
  double Pplusin=1./sqrt(2.)*(Ein+Pzin);

  //Set tmax
  double t_max=2.*Pplusin*Pminusin;

  double x_min=0.00001;

  int lund_nbins=100;
  double lund_hist[100][100]={{0.}};
  double t_binsize=log(t_max/2./t0)/double(lund_nbins);
  double z_binsize=log(1./x_min)/double(lund_nbins); 

  //Parton Shower
  for (int iEv=1; iEv<=Nev; iEv++) {
 
    int nsplits=0;
 
    vector<Parton> parton_list;
    
    //Initial Hard Parton
    //double prog[4];
    //double test_prog[4]={2.,-13.2,4.1,0.};
    //rotate_prog(test_prog);
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
          //cout << " First iteration ";
	  double mmax2=t_max;				//Max virtuality is that of hard scattering (Ein^2) 
	  //cout << " t_max= " << t_max << endl;   
 
          double mb2, zb, pp, pm;
          bool zb_ok=0;
          do {
            
            EvolveParton(parton_list[iP],mmax2);
          
            zb=parton_list[iP].z();
	    pp=parton_list[iP].pplus();
	    pm=pp;
	    mb2=parton_list[iP].virt();
            if (mb2>2.*pp*pm) { cout << " mb2 > 2 pp pm  at FIRST ITERATION " << endl; exit(0); }         
 
	    zb_ok=1;
            //Check z constraints
            if (parton_list[iP].stat()>0) {
              double zb_min=1./2.*(1.-sqrt(1.-4.*t0/mb2));
              double zb_max=1./2.*(1.+sqrt(1.-4.*t0/mb2)); 
              if (zb>zb_max || zb<zb_min) { cout << " SHOULD NOT TRIGGER Z CONSTRAIN IN FIRST ITERATION "; zb_ok=0; } 	
            }

            mmax2=mb2;

          } while (zb_ok==0);

	  //Update its pm
	  pm=mb2/2./pp;
	  FourVector p(0.,0.,1./sqrt(2.)*(pp-pm),1./sqrt(2.)*(pp+pm));
	  //prog[0]=0., prog[1]=0., prog[2]=1./sqrt(2.)*(pp-pm), prog[3]=1./sqrt(2.)*(pp+pm);
          parton_list[iP].set_p(p);
          parton_list[iP].set_alpha(1.);

	  //cout << " zb= " << zb << " virt= " << parton_list[iP].virt() << endl << endl;
	  //Add daughters
	  if (parton_list[iP].stat()>0) AddDaughters(parton_list, iP);

	}
        else
        {
	  //cout << " Doing sisters iter= " << iter << "\n \n";
          //cout << " I am iP= " << iP << endl;
	  int mom=parton_list[iP].mom();
          int sister=parton_list[mom].d1();
          if (sister==iP) sister=parton_list[mom].d2();
        
	  double mom_z=parton_list[mom].z();

          double ma2=parton_list[mom].virt();  
	  double tbmax=mom_z*ma2-t0/(1.-mom_z);
	  double tcmax=(1.-mom_z)*ma2-t0/mom_z;
	  //cout << " mom= " << mom << endl;
 	  //cout << " mom_z = " << mom_z << endl;

          int siscounter=0;
	  bool evolve_b=1, evolve_c=1;
	  do {
	    //Evolve the two sisters
	    if (evolve_b) EvolveParton(parton_list[iP],tbmax);
	    if (evolve_c) EvolveParton(parton_list[sister],tcmax);

	    evolve_b=0, evolve_c=0;

            double mb2=parton_list[iP].virt();
            double mc2=parton_list[sister].virt();

	    //Check m constraints
	    bool mb_ok=1;
            bool mc_ok=1;
	    if (pow(sqrt(mb2+t0)+sqrt(mc2+t0),2.)>ma2) {
	      double ratb=sqrt(mb2)/sqrt(tbmax);
	      double ratc=sqrt(mc2)/sqrt(tcmax);
	      if (ratb>ratc) mb_ok=0;
	      else mc_ok=0;
	      //cout << " masses not ok mb_ok= " << mb_ok << " mc_ok= " << mc_ok << endl;
	    }

	    //Check z constraints
            if (mb_ok==1 && mc_ok==1) {
	      double zconsmin=1./2./ma2*(ma2+mb2-mc2)*(1.-sqrt(1.-4.*ma2*(t0+mb2)/pow(ma2+mb2-mc2,2.)));
	      double zconsmax=1./2./ma2*(ma2+mb2-mc2)*(1.+sqrt(1.-4.*ma2*(t0+mb2)/pow(ma2+mb2-mc2,2.)));
              //cout << " zconsmax= " << zconsmax << endl;
              if (mom_z<zconsmin || mom_z>zconsmax) {
		double ratb=sqrt(mb2)/sqrt(tbmax);
	        double ratc=sqrt(mc2)/sqrt(tcmax);
	        if (ratb>ratc) mb_ok=0;
	        else mc_ok=0;
	        //cout << " mom z not ok mb_ok= " << mb_ok << " mc_ok= " << mc_ok << endl;
	      }  
	    }
	   
	    if (mb_ok==0) evolve_b=1;
	    if (mc_ok==0) evolve_c=1;

	    //Update maximum virtualities
	    tbmax=mb2;
	    tcmax=mc2;

            //Make sure we don't evolve frozen partons
	    if (parton_list[iP].stat()<0) evolve_b=0;
	    if (parton_list[sister].stat()<0) evolve_c=0;

	    siscounter++;

	  } while (evolve_b==1 || evolve_c==1);

	  //Update momentum of splitters
	  SetKinematics(parton_list[mom],parton_list[iP],parton_list[sister]);
	  //SetKinematicsSudakovBasis(parton_list[mom],parton_list[iP],parton_list[sister],prog);
   
          //Introduce daughters of b
          if (parton_list[iP].stat()>0) AddDaughters(parton_list, iP);	
          //Introduce daughters of c
          if (parton_list[sister].stat()>0) AddDaughters(parton_list, sister);
	  
        } //End if iter!=0 
      
      } //End parton_list loop
        
      iter++;

      if (changes==0) break;

    } //Infinite changes loop
    
    //cout << " Parton List Size= " << parton_list.size() << endl;
    //cout << " Nsplits= " << nsplits << endl;

    for (unsigned int i=0; i<parton_list.size(); i++) {
      //cout << i << " " << parton_list[i].p().x() << " " << parton_list[i].p().y() << " " << parton_list[i].p().z() << " " << parton_list[i].p().t() << " " << parton_list[i].d1() << " " << parton_list[i].d2() << " " << parton_list[i].mom() << " " << parton_list[i].z() << " " << parton_list[i].virt() << endl;
    }
    
    parton_list.clear();
  
    if (iEv % 1000 == 0) {

      cout << "#EVENT= " << iEv << endl;
/*      
      //Print Hist
      std::ostringstream flund;
      flund << "primary_lund_alphasconstant_" << alphas_constant << "simplesplitting_" << simple_splitting << ".dat";
      std::ofstream lundfile(flund.str().c_str(),std::ios_base::binary);
  
      for (int a=0; a<lund_nbins; a++) {
        //cout << " z_spread= " << z_spread << " t_bincen= " << t_bincen << endl;
        for (int b=0; b<lund_nbins; b++) {
          double normal_bin_area = t_binsize * z_binsize;
          lundfile << double(a)*t_binsize+t_binsize/2. << " " << double(b)*z_binsize+z_binsize/2. << " " << lund_hist[a][b]/double(iEv) << " " << lund_hist[a][b]/double(iEv)/normal_bin_area << endl;
        }
        lundfile << endl;
      }
      lundfile.close();
*/      
    }

 
  //End Event Loop
  }

  //Print Hist
  std::ostringstream flund;
  flund << "primary_lund_alphasconstant_" << alphas_constant << "simplesplitting_" << simple_splitting << ".dat";
  std::ofstream lundfile(flund.str().c_str(),std::ios_base::binary);
  
  double maxspread=6.14344;
  for (int a=0; a<lund_nbins; a++) {
    double t_bincen = double(a)*t_binsize+t_binsize/2.;
    double z_spread = log(1.-exp(-(t_bincen+log(2.)))) + (t_bincen+log(2.));
    //cout << " z_spread= " << z_spread << " t_bincen= " << t_bincen << endl;
    for (int b=0; b<lund_nbins; b++) {
      double bin_area=t_binsize * z_binsize * maxspread/z_spread;
      double normal_bin_area = t_binsize * z_binsize;
      lundfile << double(a)*t_binsize+t_binsize/2. << " " << double(b)*z_binsize+z_binsize/2. << " " << lund_hist[a][b]/double(Nev) << " " << lund_hist[a][b]/double(Nev)/normal_bin_area << " " << lund_hist[a][b]/double(Nev)/bin_area << endl;
    }
      lundfile << endl;
  }
  lundfile.close();
  
  return 0;
}
