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

void write_HepMC3_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC3::WriterAscii &outfile, int iEv) {
  //Fill particles
  std::vector<GenParticlePtr> GenVector;
  GenVector.reserve( parton_list.size() );
  for (unsigned int i = 0; i < parton_list.size(); ++i) { //FIXME I shouldn't skip the system particle
    int stat = 0;
    if (parton_list[i].stat() == -12)    stat = 4;                      //Incoming beam
    else if (parton_list[i].stat() > 0)  stat = 1;                      //Final-state particle (no further decay)
    else if (parton_list[i].stat() <= -11) stat = -parton_list[i].stat(); //Intermediate virtual particle
    GenParticlePtr p1 = std::make_shared<GenParticle>(HepMC3::FourVector(parton_list[i].px(),
                                                                         parton_list[i].py(),
                                                                         parton_list[i].pz(),
                                                                         parton_list[i].e()),
                                                                         parton_list[i].id(), stat);
    GenVector.push_back(p1);
    GenVector[i]->set_generated_mass(parton_list[i].mass()); //TODO include mass
  }

  GenEvent evt(Units::GEV,Units::MM);
  evt.set_event_number(iEv);

  //Fill vertices
  std::vector<GenVertexPtr> vertex_cache;
  std::vector<GenParticlePtr> beam_particles;
  for (unsigned int i = 1; i < parton_list.size(); ++i) { //FIXME i shouldnt skip the system particle, but otherwise crashes
    std::vector<int> mothers = parton_list[i].motherList();
    if (mothers.size()) {
      GenVertexPtr prod_vtx = GenVector[mothers[0]]->end_vertex();
      if (!prod_vtx) {
        prod_vtx = std::make_shared<GenVertex>();
        vertex_cache.push_back(prod_vtx);
        for(unsigned int j=0; j<mothers.size(); ++j) {
          prod_vtx->add_particle_in( GenVector[mothers[j]] );
        }
      }
      double GeV2MM = 0.1973 * 1.e-12;
      HepMC3::FourVector prod_pos( parton_list[i].x().x()*GeV2MM, parton_list[i].x().y()*GeV2MM,
                                   parton_list[i].x().z()*GeV2MM, parton_list[i].x().t()*GeV2MM );
      // Update vertex position if necessary
      if(!prod_pos.is_zero() && prod_vtx->position().is_zero()) prod_vtx->set_position( prod_pos );
      prod_vtx->add_particle_out( GenVector[i] );
      //TODO set signal vertex
    }
    else beam_particles.push_back(GenVector[i]);
  }

  // Reserve memory for the event
  evt.reserve( GenVector.size(), vertex_cache.size() );

  // Add particles and vertices in topological order
  if (beam_particles.size()!=2) {
    std::cerr << "There are  " << beam_particles.size() <<"!=2 particles without mothers"<< std::endl;
    exit(1);
  }
  evt.add_tree( beam_particles );

  //Attributes should be set after adding the particles to event
  for(unsigned int i = 1; i < parton_list.size(); ++i) { //FIXME I shouldn't skip the system particle
      // Colour flow uses index 1 and 2.
      int colType = 0;
      if (parton_list[i].col() == 0 && parton_list[i].acol() == 0)      colType = 0; //color-less
      else if (parton_list[i].col() == 0 && parton_list[i].acol() != 0) colType = -1; //anti-quark
      else if (parton_list[i].col() != 0 && parton_list[i].acol() == 0) colType = 1; //quark
      else if (parton_list[i].col() != 0 && parton_list[i].acol() != 0) colType = 2; //gluon
      if (colType ==  -1 ||colType ==  1 || colType == 2) {
          int flow1=0, flow2=0;
          if (colType ==  1 || colType == 2) flow1=parton_list[i].col();
          if (colType == -1 || colType == 2) flow2=parton_list[i].acol();
          GenVector[i]->add_attribute("flow1",std::make_shared<IntAttribute>(flow1));
          GenVector[i]->add_attribute("flow2",std::make_shared<IntAttribute>(flow2));
      }
  }

  // Check for particles had been forgot to add to the event
  for(unsigned int i = 1; i < parton_list.size(); ++i) {
    if ( !GenVector[i] ) {
      std::cerr << "hanging particle " << i << std::endl;
      GenVertexPtr prod_vtx;
      prod_vtx->add_particle_out( GenVector[i] );
      evt.add_vertex(prod_vtx);
    }
  }


  GenCrossSectionPtr xsec = std::make_shared<GenCrossSection>();
  //xsec is in pb
  xsec->set_cross_section( event_xsec, 0. ); //TODO include error
  evt.set_cross_section(xsec);
  evt.weights().push_back(event_weight); //FIXME many weights???

 // Print::listing(evt);
  outfile.write_event(evt);

  //std::cout << "Done hepmc" << std::endl;
  return;

}

} //end namespace Adkoda
