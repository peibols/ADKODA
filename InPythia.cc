#include "Pythia8/Event.h"

#include "InPartons.h"

using namespace Pythia8;
using namespace Adkoda;

void InPartons::PythiaInit() {

  pythia.readFile("setup_pythia.cmnd");

  std::cout << "Initialize PYTHIA." << std::endl;
  pythia.init();

}
std::vector<Parton> InPartons::PythiaPartonList() {

  pythia.next();

  //Event Info
  event_weight = pythia.info.weight();
  event_xsec = pythia.info.sigmaGen();
  hard_pt_max = pythia.info.pTHat();
  std::cout << " Pythia cross= " << event_xsec << std::endl;

  std::vector<Parton> hard_list;
  FourVector x;
  std::vector<int> incoming_pair;
  for (int i = 0; i < pythia.event.size(); i++) {
    //if (pythia.event[i].status()!=63 && pythia.event[i].status()!=-21 && pythia.event[i].status()!=-23) continue;
    if (pythia.event[i].status()!=63 && pythia.event[i].status()!=-21 && pythia.event[i].status()!=62) continue;

    FourVector p ( pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz(), pythia.event[i].e() );

    int stat=-1000;
    if (pythia.event[i].status() == 63) stat = 63;		//Beam remnant
    //if (pythia.event[i].status() == -23) stat = 23;		//Hard outgoing parton
    if (pythia.event[i].status() == 62) stat = 23;		//Hard outgoing parton
    if (pythia.event[i].status() == -21) stat = -21; 		//Hard incoming parton
    if (stat==-1000) {
      std::cout << "Unexpected stat." << std::endl;
      exit(1);
    }

    int cols[2] = {pythia.event[i].col(), pythia.event[i].acol()};

    Parton hard_parton( Parton(pythia.event[i].id(),stat,p,x) );
    hard_parton.set_cols(cols);

    hard_list.push_back(hard_parton);
    
    if (stat==-21) incoming_pair.push_back(hard_list.size()-1);
    
    int m1, m2;
    if (stat==-21 || stat==63) m1=-1, m2=-1;
    else {
      if (incoming_pair.size()<2 ) { 
        std::cout << "Did not find hard incoming pair yet." << std::endl;
	exit(1); 
      }
      m1=incoming_pair[0];
      m2=incoming_pair[1];
    }
    
    hard_list[hard_list.size()-1].set_mom1(m1); 
    hard_list[hard_list.size()-1].set_mom2(m2); 

  }

  return hard_list;
}
