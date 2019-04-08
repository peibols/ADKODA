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
  //readjust tmax
  t_max=2.*Pplusin*Pminusin-Pxin*Pxin-Pyin*Pyin;
  double zmin=t0/t_max;
  cout << " tmax= " << t_max << endl;

  //double biggest_angle=sqrt(t_max/Ein/Ein/zmin/(1.-zmin));

  int lund_nbins=100;
  double lund_hist[100][100]={{0.}};
  double t_norm[100]={0.};
  double t_binsize=log(t_max/2./t0)/double(lund_nbins);
  double z_binsize=log(1./x_min)/double(lund_nbins); 

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);

  double lowest_sudakov=Sudakov(2.*t0,t0);
  double highest_sudakov=Sudakov(t_max,t0);
  cout << " lowest_sudakov= " << lowest_sudakov << endl;
  //cout << " biggest angle= " << biggest_angle << endl;

  //Parton Shower
  for (int iEv=1; iEv<=Nev; iEv++) {
 
    int nsplits=0;
 
    vector<Parton> parton_list;
    
    //Initial Hard Parton
    FourVector phard(Pxin, Pyin, Pminusin, Pplusin);    
    FourVector x;
    Parton hard_parton( Parton(21,2,phard,x) );
    hard_parton.set_z(1.);
    hard_parton.set_mom(-3);

    parton_list.push_back(hard_parton);

    double t1=t_max;
    double x1=x_max;

    double t2=-1000.;
    double x2=-1000.;
   
    double prev_sudakov=highest_sudakov;

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

        if (iter==0) { t1=t_max; }
        else {
          int mom=parton_list[iP].mom();
	  t1=parton_list[mom].virt();
        }
        double unconst_t1=t1;


        //Correct maximum allowed virtualities (sister dependent!)
        double kin_maxt=pow(1./sqrt(2.)*(parton_list[iP].p().t()+parton_list[iP].p().z()),2.);
        if (kin_maxt<t1) {
	  //cout << " kin_maxt = " << kin_maxt << "t1= " << t1 << endl;
          t1=kin_maxt;
        }
        if (iter!=0) {
          int mom=parton_list[iP].mom();
          int sister=parton_list[mom].d1();
          if (sister==iP) sister=parton_list[mom].d2();
          double sis_virt=2.*t0;
          if (parton_list[sister].virt()!=-1000 && parton_list[sister].d1()!=-1) sis_virt=parton_list[sister].virt();
          double sum_maxt=parton_list[iP].z()*parton_list[mom].virt()
			  -parton_list[iP].z()/(1.-parton_list[iP].z())*sis_virt
			  -1./(1.-parton_list[iP].z())*t0;
          //if (sum_maxt<t1) t1=sum_maxt;
          //cout << endl;
	  if (sum_maxt<0) { 
	    sum_maxt=2.*t0;
	    //cout << " z= " << parton_list[iP].z() << endl;
	    //cout << " NEGATIVE SUM MAX T= " << sum_maxt << endl;
	    //cout << " 1st piece= " << parton_list[iP].z()*parton_list[mom].virt() << endl; 
	    //cout << " 2nd piece= " << -parton_list[iP].z()/(1.-parton_list[iP].z())*sis_virt << endl;
	    //cout << " 3rd piece= " << -1./(1.-parton_list[iP].z())*t0;
	    //exit(0);
	  }
	  if (sum_maxt<t1) t1=sum_maxt;
	  double eb=1./sqrt(2.)*(parton_list[iP].p().t()+parton_list[iP].p().z());
          double ea=1./sqrt(2.)*(parton_list[mom].p().t()+parton_list[mom].p().z());
          if (sum_maxt>ea*ea*parton_list[iP].z() && sum_maxt>2.*t0) {
	    cout << " z= " << parton_list[iP].z() << endl; 
	    cout << " kin_maxt= " << kin_maxt << " sum_maxt= " << sum_maxt << endl;
            cout << " eb*ea= " << eb*ea << endl; 
            cout << " z ea^2 = " << ea*ea*parton_list[iP].z() << endl;
	    //cout << " first piece= " << parton_list[iP].z()*parton_list[mom].virt() << " second piece= " << parton_list[iP].z()/(1.-parton_list[iP].z())*sis_virt << endl;
	    cout << " mom virt= " << parton_list[mom].virt() << " sis virt= " << sis_virt << endl;
	  }
          cout << " unconst_t1= " << unconst_t1 << " kinmaxt= " << kin_maxt << " sum_maxt= " << sum_maxt << endl;
	}

        double R=dis(gen);

	//cout << " t1= " << t1 << endl;
        if (iter!=0 && t1>2.*t0) prev_sudakov=Sudakov(t1,t0);
        if (R > prev_sudakov/lowest_sudakov && t1>2.*t0)
        {
          //Generate t
	  // Use hit and miss method, to compare against CDF method
          double t2_temp=generate_t2(t0,t1,R,log(prev_sudakov));
	 
/*
	  double mcsuda=0.;
          double randian=1.;
          double t2_rand=-1.;
          do {
            randian=dis(gen);
            t2_rand=2.*t0+dis(gen)*(t1-2.*t0);
            mcsuda=prev_sudakov/Sudakov(t2_rand,t0);
            cout << " mcsuda= " << mcsuda << endl;
	  } while (mcsuda<randian);
	  cout << "t1= " << t1 << " t2rand= " << t2_rand << endl;
	  double t2_temp=t2_rand;
*/
          if (t2_temp>2.*t0) t2=t2_temp;
          else {
	    cout << " TRIGGERED \n \n";
            parton_list[iP].set_stat(-1);
            parton_list[iP].set_d1(-1);
            parton_list[iP].set_d2(-1);
            continue;
          }

          //Generate z
          double Rp=dis(gen);
          double xi=generate_x2(t0,t2,Rp);

	  //cout << " Splitting happened " << endl;
	  nsplits+=1;

	  //Generate daughters (provisional tri-momentum, until their virtualities and azimuthal angle are generated
          FourVector p1(0.,0.,0.,xi*parton_list[iP].p().t());
	  FourVector p2(0.,0.,0.,(1.-xi)*parton_list[iP].p().t());
	  FourVector x1;
	  FourVector x2;

	  Parton d1 = Parton(21,1,p1,x1);
	  Parton d2 = Parton(21,1,p2,x2);
	  
          //Follow the heir
	  if (parton_list[iP].stat()==2) {
            if (xi>0.5) {
	      d1.set_stat(2);
	    }
	    else { 	
	      d2.set_stat(2);
            }
          }
	
	  d1.set_mom(iP);
	  d2.set_mom(iP);
          d1.set_z(xi);
	  d2.set_z(1.-xi);

	  parton_list.push_back(d1);         
	  parton_list[iP].set_d1(parton_list.size()-1); 
          
	  parton_list.push_back(d2);         
	  parton_list[iP].set_d2(parton_list.size()-1); 
	  
	  //Assign parton virtuality
          //cout << " t= " << t2 << " z= " << xi << endl;
          //Fill Hist
	  if (parton_list[iP].stat()==2) {
	    int tbin=int(log(t2/2./t0)/t_binsize);
	    int zbin=int(log(1./xi)/z_binsize);
	    if (tbin < lund_nbins && zbin < lund_nbins) {
	      lund_hist[tbin][zbin]+=1.;
	      t_norm[zbin]+=1.;
	    }
	    //Angle 
            double angle_dla=sqrt(t2/(xi*(1.-xi))/pow(parton_list[iP].p().t(),2.));
            //cout << " t= " << t2 << " z= " << xi << " angle_dla= " << angle_dla << endl;
          }

	  parton_list[iP].set_virt(t2);
	  parton_list[iP].set_stat(-1);
	 
        }
        else {
          //cout << " freezing with max virt= " << t1 << endl << endl;
          parton_list[iP].set_virt(2.*t0);
          parton_list[iP].set_stat(-1);
          parton_list[iP].set_d1(-1);
          parton_list[iP].set_d2(-1);
          continue;
	}
        //cout << endl;

      }

      iter++;

      if (changes==0) break;

    }
    
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
