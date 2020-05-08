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

  // D(x) histo
  int ntau_cuts=4;
  double tau_cuts[4]={0.01,0.02,0.1,0.2};
  int max_nx=100;
  double xbin = 1./double(max_nx);
  double dx_hist[100][4]={{0.}};

  BerGEN bergen(input_file);
  int nEv = bergen.number_events();
  bergen.init();
  for (int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 1 == 0) cout << "\n #Event: " <<  iEv << endl;
    bergen.next();
    bergen.print();

    //Test modules
    //double event_weight = bergen.get_event_weight();
    //double event_xsec   = bergen.get_event_xsec();
    std::vector<Parton> parton_list = bergen.get_parton_list();
    //Test_Weights(parton_list, event_xsec, event_weight, outfile_test_weights);
    //Test_PrintFinalPartons(parton_list, event_xsec, event_weight, outfile_test_FinalPartons);
    //Test_EnergyMomentumConservation(parton_list);
    //Test_PrintLundPlane_history(parton_list, event_xsec, event_weight, outfile_test_LundPlane);

    //HepMC3 ASCII writer
    //write_HepMC3_event(parton_list, event_xsec, event_weight, outfile_test_HepMC3, iEv);

    // Freezing time hist
    double pplus;
    double alphas_med = 0.3;
    double qhat = 1.5 * 0.1973;
    for (unsigned int i=0; i<parton_list.size(); i++) {
      Parton p = parton_list[i];
      if (p.stat()==-23) { pplus=p.e(); continue; }
      if (p.x().t()==0.) continue;
    
      double ti = p.x().t();
      double taui = alphas_med * std::sqrt(qhat/pplus) * ti / 0.1973;

      double tf = p.xf().t();
      if (tf==-1000.) tf = 10000000.; 
      double tauf = alphas_med * std::sqrt(qhat/pplus) * tf / 0.1973;

      //cout << " pplus= " << pplus << " taui= " << taui << " tauf= " << tauf << endl;

      for (int c=0; c<ntau_cuts; c++) {
        if (taui<=tau_cuts[c] && tauf>tau_cuts[c]) {
          double x = p.e()/pplus;
          int nx = int(x/xbin);
          //cout << " x= " << x << " nx= " << nx << endl;
          if (nx>=0 && nx<max_nx) {
            dx_hist[nx][c]+=1.;
          }
        }
      }      
    }

  }
  cout << "#End Program" << endl;

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

  //Closing outfiles
  outfile_test_weights.close();
  outfile_test_FinalPartons.close();
  outfile_test_LundPlane.close();

  return 0;
} // End program
