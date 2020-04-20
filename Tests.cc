#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>

#include "Tests.h"
#include "BerGEN.h"
#include "Util.h"

namespace Adkoda {

//TEST MODULE: Save the initiator e-e+ kinematics and event_weights
void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  outfile << parton_list[4].e()  << " " << parton_list[4].px() << " " <<
             parton_list[4].py() << " " << parton_list[4].pz() << " " <<
             event_xsec << event_weight << std::endl;
}

//TEST MODULE: Save final particles (px, py, pz, E)
void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int multi = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    //if (parton_list[ip].stat()<0 && abs(parton_list[ip].id())!=11) continue; //MCnet python durham algorithm needs the incoming  particles.
    if (parton_list[ip].stat() < 0) continue;
     outfile << parton_list[ip].px() << " " << parton_list[ip].py() << " "
             << parton_list[ip].pz() << " " << parton_list[ip].e() << std::endl;
     multi++;
  }
  outfile << "#event xsec: " << event_xsec << ", weight: " << event_weight << ", multiplicity: " << multi << std::endl;
}


void Test_EnergyMomentumConservation(std::vector<Parton> parton_list) {
  double total_p[4] = {0., 0., 0., 0.};
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) {
    if (parton_list[ip].stat()<0) continue; //MCnet python durham algorithm needs the incoming  particles.
    total_p[0] += parton_list[ip].e();
    total_p[1] += parton_list[ip].px();
    total_p[2] += parton_list[ip].py();
    total_p[3] += parton_list[ip].pz();
  }
  for (unsigned a = 0; a < 4; a++) std::cout << "Total momentum: p_" << a << " = " << total_p[a] << std::endl;
}

//TEST MODULE: Save the primary Lund plane using the history and FastJet definition of the variables
void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int m1 = -1, m2 = -1;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (abs(parton_list[ip].stat()) == 23 && m1 == -1) m1 = ip;
    if (abs(parton_list[ip].stat()) == 23 && ip != m1 && m2 == -1) { m2 = ip; break; }
  }
  if (parton_list[m1].stat() < 0) { //It decays
    int mom = m1;
    int d1  = parton_list[mom].d1();
    int d2  = parton_list[mom].d2();
    if (d1==d2) mom = m2; //Find the initial splitter
    //Soft Drop parameters
    double z_cut = 0.;
    double beta = 0.;
    double R = 3.;
    while (parton_list[mom].stat() < 0) {
      d1 = parton_list[mom].d1();
      d2 = parton_list[mom].d2();
      if (abs(parton_list[d1].stat()) == 51) { //A splitter
        if (parton_list[d1].pt() < parton_list[d2].pt()) { //Follow the primary
          int temp = d1;
          d1 = d2;
          d2 = temp;
        }
        double Delta    = Util::Delta(parton_list[d1].p(), parton_list[d2].p());
        double z        = parton_list[d2].pt() / (parton_list[d1].pt() + parton_list[d2].pt());
        if (z > z_cut * pow(Delta/R, beta) && Delta < R) { // Soft Drop condition.
          double kt     = parton_list[d2].pt() * Delta;
          double m2     = Util::m2(parton_list[d1].p(), parton_list[d2].p()); //Mass of the splitter
          double qt2    = m2 / z / (1.-z);
          double energy = parton_list[d1].e() + parton_list[d2].e(); //Energy of the splitter
          double tf     = 2. * energy / m2;
          outfile << qt2 << " " << z << " " << energy << " " << kt*kt << " "
                  << Delta << " " << m2 << " " << tf << " " << "1." << std::endl;
        }
        mom = d1;
      }
      else if (parton_list[d1].stat()>0) break;
      else mom = d1;
    }
  }
  outfile << "#xsec: " << event_xsec << "weight: " << event_weight << std::endl;
}

/*
//TEST MODULE: Quench particle history by counting the splitting inside the medium
void Simple_Quenching(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile) {
  int inSplit = 0;
  for (unsigned int ip = 0; ip < parton_list.size(); ip++) { //Find the hardest scattering
    if (parton_list[ip].x().t()*0.1973 < Rmed_in_fm) inSplit++;
  }
  Test_PrintFinalPartons(parton_list, event_xsec, std::pow(Quench, inSplit), outfile);
}
*/

} // end namespace Adkoda
