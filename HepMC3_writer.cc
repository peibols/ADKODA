// Following status codes, as much as possible, from Appendix A in https://arxiv.org/pdf/1912.08005.pdf

#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
#include "HepMC3/Relatives.h"

#include "Parton.h"

using namespace HepMC3;

namespace Adkoda {

void write_HepMC3_event(std::vector<Parton> parton_list, double event_weight, HepMC3::WriterAscii &outfile) { 

  vector<GenParticlePtr> GenVector;
  for (int p=0; p<parton_list.size(); p++) {
    int stat=2;											//Decayed physical particle
    if (parton_list[p].d1() == parton_list[p].d2() && parton_list[p].d1() != 0) stat=52; 	//Decays into carbon copy
    if (parton_list[p].d1() == 0 && parton_list[p].d2() == 0) stat=1;				//Undecayed physical particle
    GenParticlePtr p1 = make_shared<GenParticle>(HepMC3::FourVector(parton_list[p].px(),parton_list[p].py(),parton_list[p].pz(),parton_list[p].e()),parton_list[p].id(),stat);
    GenVector.push_back(p1);
  }
  
  GenEvent evt(Units::GEV,Units::MM);

  for (int p=0; p<parton_list.size(); p++) {
    int d1=parton_list[p].d1();
    int d2=parton_list[p].d2();
    if (d1==0 && d2==0) continue; //Final Particle
    GenVertexPtr v1 = make_shared<GenVertex>();
    v1->add_particle_in (GenVector[p]);
    evt.add_vertex(v1);
    v1->add_particle_out (GenVector[d1]);
    int stat;
    if (d1 != d2 && d2 != 0) v1->add_particle_out (GenVector[d2]), stat=1; //Add second daughter if not carbon copy
    else stat=2;							   //Carbon copy
    v1->set_status(stat);
  }

  evt.weights().push_back(event_weight);

  outfile.write_event(evt);

  return;

}

} //end namespace Adkoda
