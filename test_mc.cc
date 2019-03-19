#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <assert.h>
#include "Parton.h"

using namespace std;

double Sudakov(double t, double t0);
double generate_t2(double t0, double t1, double R, double numer);
double generate_x2(double t0, double t2, double R);

int main(int argc, char** argv) {

  assert(argc==2);

  int Nev=atoi(argv[1]);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.,1.);

  double t_max=1000.;
  double x_max=1.;
  double t0=2.;

  double lowest_sudakov=Sudakov(2.*t0,t0);
  double highest_sudakov=Sudakov(t_max,t0);

  //Parton Shower
  for (int iEv=1; iEv<=Nev; iEv++) {
 
    int nsplits=0;
 
    vector<Parton> parton_list;
    
    //Initial Hard Parton
    double Ein=1000.;
    double Pxin=0.;
    double Pyin=0.;
    double Pzin=0.;
    FourVector phard(Pxin, Pyin, Pzin, Ein);    
    FourVector x;
    Parton hard_parton( Parton(21,1,phard,x) );
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

        if (iter==0) t1=t_max;
        else {
          int mom=parton_list[iP].mom();
	  t1=parton_list[mom].virt();
        }

        double R=dis(gen);

        if (iter!=0) prev_sudakov=Sudakov(t1,t0); 
        if (R>prev_sudakov/lowest_sudakov)
        {
          //Generate t
          double t2_temp=generate_t2(t0,t1,R,log(prev_sudakov));

          if (t2_temp>2.*t0) t2=t2_temp;
          else {
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

	  Parton d1( Parton(21,1,p1,x1) );
	  Parton d2( Parton(21,1,p2,x2) );
		
	  d1.set_mom(iP);
	  d2.set_mom(iP);
          d1.set_z(xi);
	  d2.set_z(1.-xi);

	  parton_list.push_back(d1);         
	  parton_list[iP].set_d1(parton_list.size()-1); 
          
	  parton_list.push_back(d2);         
	  parton_list[iP].set_d2(parton_list.size()-1); 
	  
	  //Assign parton virtuality
          cout << " t= " << t2 << " z= " << xi << endl;
          parton_list[iP].set_virt(t2);
	  parton_list[iP].set_stat(-1);
	 
        }
        else {
          parton_list[iP].set_stat(-1);
          parton_list[iP].set_d1(-1);
          parton_list[iP].set_d2(-1);
          continue;
	}

      }

      iter++;

      if (changes==0) break;

    }
    
    cout << " Parton List Size= " << parton_list.size() << endl;
    cout << " Nsplits= " << nsplits << endl;
    
    parton_list.clear();
    cout << "#EVENT= " << iEv << endl;
  
  //End Event Loop
  }

  
  return 0;
}
