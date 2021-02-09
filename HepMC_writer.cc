#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Units.h"

#include "Parton.h"

using namespace HepMC;

namespace Adkoda {

void write_HepMC_event(std::vector<Parton> parton_list, double event_xsec, double event_weight, HepMC::IO_GenEvent &outfile, int iEv) {

  HepMC::GenEvent* evt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  evt->set_event_number(iEv);

  // 2. Create a particle instance for each entry and fill a map, and
  // a vector which maps from the particle index to the GenParticle address.
  std::vector<GenParticle*> hepevt_particles( parton_list.size() );
  for (unsigned int i = 0; i < parton_list.size(); ++i) { //FIXME I shouldn't skip the system particle
    int stat = 0;
    if (parton_list[i].stat() == -12)    stat = 4;                      //Incoming beam
    else if (parton_list[i].stat() > 0)  stat = 1;                      //Final-state particle (no further decay)
    else if (parton_list[i].stat() <= -11) stat = -parton_list[i].stat(); //Intermediate virtual particle
    // Fill the particle.
    hepevt_particles[i] = new GenParticle(
            HepMC::FourVector( parton_list[i].px(), parton_list[i].py(),
                               parton_list[i].pz(), parton_list[i].e()  ),
                               parton_list[i].id(), stat );
    hepevt_particles[i]->suggest_barcode(i);
    hepevt_particles[i]->set_generated_mass(parton_list[i].mass()); //TODO include mass

    // Colour flow uses index 1 and 2.
    int colType = 0;
    if (parton_list[i].col() == 0 && parton_list[i].acol() == 0)      colType = 0; //color-less
    else if (parton_list[i].col() == 0 && parton_list[i].acol() != 0) colType = -1; //anti-quark
    else if (parton_list[i].col() != 0 && parton_list[i].acol() == 0) colType = 1; //quark
    else if (parton_list[i].col() != 0 && parton_list[i].acol() != 0) colType = 2; //gluon

    // Colour flow uses index 1 and 2.
    if (colType ==  1 || colType == 2)
      hepevt_particles[i]->set_flow(1, parton_list[i].col());
    if (colType == -1 || colType == 2)
      hepevt_particles[i]->set_flow(2, parton_list[i].acol());
  }

  // Here we assume that the first two particles in the list
  // are the incoming beam particles.
  evt->set_beam_particles( hepevt_particles[1], hepevt_particles[2] ); //FIXME Assigning by hand

  // 3. Loop over particles AGAIN, this time creating vertices.
  // We build the production vertex for each entry in hepevt.
  // The HEPEVT pointers are bi-directional, so gives decay vertices as well.
  GenParticle* rootParticle = 0;
  bool append = false;
  for (unsigned int i = 0; i < parton_list.size(); ++i) { //FIXME I shouldn't skip the system particle
    GenParticle* p = hepevt_particles[i];
    // 3a. Search to see if a production vertex already exists.
    std::vector<int> mothers = parton_list[i].motherList();
    unsigned int imother = 0;
    int mother = -1;
    if ( !mothers.empty() ) mother = mothers[imother];
    GenVertex* prod_vtx = p->production_vertex();
    while ( !prod_vtx && mother > 0 ) {
      prod_vtx = (append && mother == 1) ? rootParticle->end_vertex()
               : hepevt_particles[mother]->end_vertex();
      if (prod_vtx) prod_vtx->add_particle_out( p );
      mother = ( ++imother < mothers.size() ) ? mothers[imother] : -1;
    }
    // 3b. If no suitable production vertex exists - and the particle has
    // at least one mother or position information to store - make one.
    double GeV2MM = 0.1973 * 1.e-12;
    HepMC::FourVector prod_pos( GeV2MM * parton_list[i].x().x(), GeV2MM * parton_list[i].x().y(),
                                GeV2MM * parton_list[i].x().z(), GeV2MM * parton_list[i].x().t() );
    if ( !prod_vtx && ( mothers.size() > 0 || prod_pos != HepMC::FourVector() ) ) {
      prod_vtx = new GenVertex();
      prod_vtx->add_particle_out( p );
      evt->add_vertex( prod_vtx );
    }
    // 3c. If prod_vtx doesn't already have position specified, fill it.
    if ( prod_vtx && prod_vtx->position() == HepMC::FourVector() )
          prod_vtx->set_position( prod_pos );
    // 3d. loop over mothers to make sure their end_vertices are consistent.
    imother = 0;
    mother = -1;
    if ( !mothers.empty() ) mother = mothers[imother];
    while ( prod_vtx && mother > 0 ) {
      // If end vertex of the mother isn't specified, do it now.
      GenParticle* ppp = (append && mother == 1) ? rootParticle
                       : hepevt_particles[mother];
      if ( !ppp->end_vertex() ) {
        prod_vtx->add_particle_in( ppp );
        // Problem scenario: the mother already has a decay vertex which
        // differs from the daughter's production vertex. This means there is
        // internal inconsistency in the HEPEVT event record. Print an error.
        // Note: we could provide a fix by joining the two vertices with a
        // dummy particle if the problem arises often.
      } else if (ppp->end_vertex() != prod_vtx )
          std::cout << " BerGenToHepMC inconsistent mother/daugher " << std::endl;

      // End of vertex-setting loops.
      mother = ( ++imother < mothers.size() ) ? mothers[imother] : -1;
    }
  }

  // 4. Check for particles which come from nowhere, i.e. are without
  // mothers or daughters. These need to be attached to a vertex, or else
  // they will never become part of the event.
  for (unsigned int i = 1; i < parton_list.size(); ++i) {
    if ( !hepevt_particles[i]->end_vertex() && !hepevt_particles[i]->production_vertex() ) {
      std::cout << " BerGenToHepMC error: " << "hanging particle " << i << std::endl;
      GenVertex* prod_vtx = new GenVertex();
      prod_vtx->add_particle_out( hepevt_particles[i] );
      evt->add_vertex( prod_vtx );
    }
  }

  //FIXME include hadronization!
  //FIXME indlude pdf weights
  //FIXME include signal process, event scale, alphaQED, alphaQCD
  HepMC::GenCrossSection xsec;
  xsec.set_cross_section( event_xsec * 1e9, 0.); //FIXME include err
  evt->set_cross_section(xsec);
  evt->weights().push_back( event_weight );

  //evt->print();
  outfile << evt;

  //std::cout << "Done hepmc" << std::endl;
  return;

}

} //end namespace Adkoda
