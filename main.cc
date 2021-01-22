#include <string>
#include <fstream>
#include "BerGEN.h"
#include "Tests.h"
#include "HepMC3/WriterAscii.h"

using namespace Adkoda;
using namespace Util;
using namespace std;

int main(int argc, char **argv) {

  //Read arguments and pass as input file
  assert(argc==2);
  string input_file = *(argv+1);
  
  remove("test/test_weights.out");
  remove("test/test_FinalPartons.out");
  remove("test/test_LundPlane.out");
  remove("test/test_kinematics.out");
  remove("test/test_veto.out");
  remove("test/test_HepMC3.out");

  HepMC3::WriterAscii outfile_test_HepMC3("test/test_HepMC3.hepmc");
  ofstream outfile_test_weights;
  if (!outfile_test_weights.is_open()) outfile_test_weights.open("test/test_weights.out", ios_base::app);
  ofstream outfile_test_FinalPartons;
  if (!outfile_test_FinalPartons.is_open()) outfile_test_FinalPartons.open("test/test_FinalPartons.out", ios_base::app);
  ofstream outfile_test_LundPlane;
  if (!outfile_test_LundPlane.is_open()) outfile_test_LundPlane.open("test/test_LundPlane.out", ios_base::app);

  cout << "#Start Program" << endl;

  int ntau_cuts=6;
  double tau_cuts[6]={0.01,0.02,0.1,0.2,1.,2.};
  
  // D(x) histo
  int max_nx=100;
  double xbin = 1./double(max_nx);
  double dx_hist[100][6]={{0.}};

  // Kx histo
  int max_kx=100;
  double kxbin = 20./max_kx;
  double kx_hist[101][6]={{0.}}; 

  // Mass histo
  double njets=0.;
  double mass_hist[2][20]={{0.}};

  BerGEN bergen(input_file);
  
  // Print to file PU14
  int evol_var=bergen.evol_var();
  std::ostringstream ff;
  ff << "final_partons" << evol_var << ".dat";
  std::ofstream pu14file(ff.str().c_str(),std::ios_base::binary);
 
  // Print to file Polish
  std::ostringstream fp;
  fp << "shower_partons_evol_var_" << evol_var << ".dat";
  std::ofstream polishfile(fp.str().c_str(),std::ios_base::binary);
 
  std::ostringstream fm;
  fm << "test/test_FirstMass" << evol_var << ".out";
  remove(fm.str().c_str());
  ofstream outfile_test_FirstMass;
  if (!outfile_test_FirstMass.is_open()) outfile_test_FirstMass.open(fm.str().c_str(), ios_base::app);

  int nEv = bergen.number_events();
  bergen.init();
  for (int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 100 == 0) cout << "\n #Event: " <<  iEv << endl;
    bergen.next();
    //bergen.print();

    //Test modules
    double event_weight = bergen.get_event_weight();
    //double event_xsec   = bergen.get_event_xsec();
    std::vector<Parton> parton_list = bergen.get_parton_list();
    //Test_Weights(parton_list, event_xsec, event_weight, outfile_test_weights);
    //Test_PrintFinalPartons(parton_list, event_xsec, event_weight, outfile_test_FinalPartons);
    //Test_EnergyMomentumConservation(parton_list);
    //Test_PrintLundPlane_history(parton_list, event_xsec, event_weight, outfile_test_LundPlane);

    //Test_FirstMass(parton_list, outfile_test_FirstMass, mass_hist, njets);

    //HepMC3 ASCII writer
    //write_HepMC3_event(parton_list, event_xsec, event_weight, outfile_test_HepMC3, iEv);
/*
    double pplus;
    double alphas_med = 0.3;
    double qhat = 1.5 * 0.1973;
    double angle, k[3];
    for (unsigned int i=0; i<parton_list.size(); i++) {
      Parton p = parton_list[i];
      if (p.stat()==-23) { 
        FourVector pIni = p.p();
        AlignWithZ(pIni, angle, k);
        pplus=1./2.*(pIni.t()+pIni.z());
        continue;
      }
      if (p.x().t()==0.) continue;
    
      double ti = p.x().t();
      double taui = alphas_med * std::sqrt(qhat/pplus) * ti / 0.1973;

      double tf = p.xf().t();
      if (tf==-1000.) tf = 10000000.; 
      double tauf = alphas_med * std::sqrt(qhat/pplus) * tf / 0.1973;

      //cout << " pplus= " << pplus << " taui= " << taui << " tauf= " << tauf << endl;

      for (int c=0; c<ntau_cuts; c++) {
        if (taui<=tau_cuts[c] && tauf>tau_cuts[c]) {
          FourVector pVec = p.p();
          Rotation(pVec, angle, k);
          // x hist
          double t_pplus = 1./2.*(pVec.t()+pVec.z());
          double x = t_pplus/pplus;
          int nx = int(x/xbin);
          //cout << " x= " << x << " nx= " << nx << endl;
          if (nx>=0 && nx<max_nx) {
            dx_hist[nx][c]+=1.;
          }
          // kx hist
          double kx = pVec.x();
          double sign=0.;
          if (kx>0.) sign=1.;
          else sign=-1.;
          int nkx = max_kx/2+int((kx+sign*kxbin/2.)/kxbin);
          if (nkx>=0 && nkx<=max_kx) {
            kx_hist[nkx][c]+=1.;
          }
        }
      }      
    }
*/    
    // Print to PU14 file
    pu14file << "# event " << iEv << endl;
    pu14file << "weight " << event_weight << endl;
    double total[3]={0.};
    for (unsigned int i=0; i<parton_list.size(); i++) {
      Parton p = parton_list[i];
      if (p.stat()<0) continue;
      int label=0;
      total[0]+=p.px();
      total[1]+=p.py();
      total[2]+=p.pz();
      pu14file << p.px() << " " << p.py() << " " << p.pz() << " " << p.mass() << " " << p.id() << " " << label << endl;
    }
    pu14file << "end" << endl;
    //for (unsigned a=0; a<3; a++) cout << " Total " << a << " = " << total[a] << endl;

    // Print to Polish
    polishfile << "# event " << iEv << endl;
    // pu14file << "weight " << event_weight << endl;
    for (unsigned int i=0; i<parton_list.size(); i++) {
      Parton p = parton_list[i];
      if (p.stat()<0) continue;
      polishfile << p.px() << " " << p.py() << " " << p.pz() << " " << p.e() << " "
		 << p.x().x() << " " << p.x().y() << " " << p.x().z() << " " << p.x().t() << " " << p.id() << endl;
    }
    polishfile << "end" << endl;

  }
  cout << "#End Program" << endl;

  pu14file.close();

  // Print mass hist
  for (int a=0; a<20; a++) {
    outfile_test_FirstMass << 1.+double(a)*2. << " " << mass_hist[0][a]/2./njets << " " << mass_hist[1][a]/2./njets << endl;
  }

  // Print Dx hist
  std::ofstream dxfile("Dx_hist.dat");
  for (int a=0; a<max_nx; a++) {
    double xcen = xbin*double(a)+xbin/2.;
    dxfile << xcen << " ";
    for (int b=0; b<ntau_cuts; b++) {
      dxfile << dx_hist[a][b]/xbin*xcen*sqrt(xcen)/nEv << " " << sqrt(dx_hist[a][b])/xbin*xcen*sqrt(xcen)/nEv << " ";
    }
    dxfile << endl;
  }
  dxfile.close();

  // Print Kx hist
  std::ofstream kxfile("Kx_hist.dat");
  for (int a=0; a<=max_kx; a++) {
    double xcen = -10.1+kxbin*double(a)+kxbin/2.;
    kxfile << xcen << " ";
    for (int b=0; b<ntau_cuts; b++) {
      kxfile << kx_hist[a][b]/kxbin/nEv << " " << sqrt(kx_hist[a][b])/kxbin/nEv << " ";
    }
    kxfile << endl;
  }
  kxfile.close();
  
  //Closing outfiles
  outfile_test_weights.close();
  outfile_test_FinalPartons.close();
  outfile_test_LundPlane.close();

  return 0;
} // End program
