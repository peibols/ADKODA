#include <iostream>
#include <fstream>

#include "Shower.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

bool Shower::evolve() {

  double wmar2, wzp;
  int wSplit, wSpect, wKernel;

  double t = t_min;
  for (unsigned int iSplit=0; iSplit<parton_list.size(); iSplit++) { // Praticle list scan
    Parton p_split = parton_list[iSplit];
    if (p_split.stat()<0) continue; // Skip if inactive parton
    for (unsigned int iSpect=0; iSpect<parton_list.size(); iSpect++) { // Find dipole
      if (iSplit == iSpect) continue;
      Parton p_spect = parton_list[iSpect];
      if (p_spect.stat()<0) continue; // Skip if inactive
      if (!p_split.ColourConnected(p_spect)) continue; // Skip if not colour connected
      for (int iKernel=0; iKernel<kernels.size(); iKernel++) { // Loop over possible kernels
        if (kernels[iKernel]->flav(0) != p_split.id()) continue; // Skip if kernel not applies
        double mar2 = m2(p_split.p(), p_spect.p());
        if (mar2 < 4.*pt_min*pt_min) continue; // Skip if dipole has no enough phase-space
	      // Overshoot z
        double epsilon = pt_min*pt_min / mar2; //This overshoots all ordering (except AO).
        if (DATA.evol_scale == 3) epsilon = std::sqrt(pt_min*pt_min / mar2);
        double zp = 0.5 * (1. + std::sqrt(1. - 4.*epsilon));
	      // Generate overshooted scale
        double g = max_alpha_s / (2.*M_PI) * kernels[iKernel]->Integral(1.-zp, zp);
        double tt = t_max * std::pow(dis(gen), 1./g);
        if ( tt > t ) { // Store if highest scale from all particles all recoilers and all possible splittings.
          t       = tt;
      	  wSplit  = iSplit;
      	  wSpect  = iSpect;
      	  wKernel = iKernel;
      	  wmar2   = mar2;
      	  wzp     = zp;
        }
      } // End Kernel loop
    } // End recoiler loop
  } // End parton_list loop
  t_max = t;
  // Terminate shower if no winner found
  if (t <= t_min) return 0;
  // Generate final z value
  double z = kernels[wKernel]->GenerateZ(1.-wzp, wzp, dis(gen));
  // Compute virtuality
  double Q2 = 0.;
  if      (DATA.evol_scale == 0) Q2 = t / z / ( 1.0 - z );
  else if (DATA.evol_scale == 1) Q2 = t;
  else if (DATA.evol_scale == 2) Q2 = t * 2.0 * std::sqrt(wmar2)/2.0;
  else if (DATA.evol_scale == 3) Q2 = t * z * (1.0 - z);
  else std::cout << "ERROR: evol_scale = 0-3" << std::endl;

  double pt2ev    = Q2 * z * (1.-z);
  double pt2_daug = Q2 * ( z * (1.-z) * std::pow(wmar2 + Q2, 2.) - wmar2 * Q2 ) / std::pow( wmar2 - Q2, 2.);
  double y = 0.;
  if (DATA.shower_kernel==1) y = t/wmar2/z/(1.-z);
  if (y < 1.) {
    // Accept / Reject veto procedure
    double f = (1.-y) * alpha_s(pt2ev) * kernels[wKernel]->Value(z, y);
    double g = max_alpha_s * kernels[wKernel]->Estimate(z);
    if (DATA.shower_kernel==0) { // If does not succeed, have a next try
      if (f/g < dis(gen) || Q2 > wmar2 || pt2ev < pt_min*pt_min || pt2_daug < 0.0 ) {return 1;}  //veto, scale in between the min/max scale, physical pt.
    } else if (DATA.shower_kernel==1) {
      if (f/g < dis(gen) || y > 1.) {return 1;}
    } else std::cout << "ERROR: shower_kernel = 0,1" << std::endl;

    // If success, set new kinematics of system
    Update( wSplit, wSpect, wKernel, wmar2, z, y, Q2 );//FIXME Why is there no additional t --> z^2*t decrision?

/*
    //TEST MODULE: print all (z,t)
    std::ofstream outfile_test;
    if (!outfile_test.is_open()) outfile_test.open("test/test_veto.out", std::ios_base::app);
    outfile_test  << z << " " << Q2 << " " << pt2ev << " " << Q2/z/(1.-z) << " " << Q2/2./(std::sqrt(wmar2)/2.) << " " << wmar2 << " " << y << std::endl;
*/
/*
    //TEST MODULE: stop at the first branch and print veto variables and kinematics
    std::ofstream outfile_test_veto;
    if (!outfile_test_veto.is_open()) outfile_test_veto.open("test/test_veto.out", std::ios_base::app);
    std::ofstream outfile_kinem;
    if (!outfile_kinem.is_open()) outfile_kinem.open("test/test_kinematics.out", std::ios_base::app);
    if (wSplit==5 || wSplit==4) {
    //if (parton_list[wSplit].mom1()==5 || parton_list[wSplit].mom2()==4) {
      outfile_test_veto << z << " " << Q2 << " " << pt2ev << " " << Q2/z/(1.-z) << " " << Q2/2./(std::sqrt(wmar2)/2.) << " " << wmar2 << " " << y << std::endl;

      outfile_kinem << parton_list[parton_list[wSplit].d2()].e()  << " " << parton_list[parton_list[wSplit].d2()].px() << " " <<
                       parton_list[parton_list[wSplit].d2()].py() << " " << parton_list[parton_list[wSplit].d2()].pz() << " " << std::endl;
      return 0;
    }
*/
    return 1;
  }
  return 1;
}

void Shower::Update( int Split, int Spect, int Kernel, double mar2, double z, double y, double Q2 ) {

  FourVector pSplit = parton_list[ Split ].p();
  FourVector pSpect = parton_list[ Spect ].p();
  FourVector UpSpect, pDau1, pDau2;
  double pt_angle = 2. * M_PI * dis(gen); //FIXME include polarization: 10.1146/annurev.ns.36.120186.001345
  if (DATA.shower_kernel==0) {
    // COM frame
    double beta[3] = {0.};
    VelCOM(pSplit, pSpect, beta); // Find boost vector to COM.

    // Boost to COM frame
    FourVector BpSplit = Boost( beta, pSplit );

    // Perform rotation of emitter along z, store angle and axis
    double angle = -1000.;
    double k[3]  = {0.};
    AlignWithZ( BpSplit, angle, k ); // Find proper andle and k to rotate around.

    // Update kinematics of spectator
    UpSpect.Set( 0.0 , 0.0 ,
  			         -( mar2 - Q2 ) / 2.0 / std::sqrt( mar2 ),
  			         ( mar2 - Q2 ) / 2.0 / std::sqrt( mar2 ) );

    // Kinematics of daughters
    double Ea = ( mar2 + Q2 ) / 2.0 / std::sqrt( mar2 );
    double E1 = z * Ea;
    double E2 = ( 1.0 - z ) * Ea;

    double pt2 = Q2 * ( z * ( 1.0 - z ) * std::pow( mar2 + Q2 , 2.0 ) - mar2 * Q2 ) / std::pow( mar2 - Q2 , 2.0 );
    if (pt2<0.) { std::cout << "Error: 0 > pt = " << pt2 << std::endl; exit(0); }

    double Pza = ( mar2 - Q2 ) / 2.0 / std::sqrt( mar2 );
    double pz1 = std::sqrt( Ea * Ea * z * z - pt2 );
    double pz2 = Pza - pz1;

    pDau1.Set( std::sqrt( pt2 ) * std::cos( pt_angle ),
  			       std::sqrt( pt2 ) * std::sin( pt_angle ),
  			       pz1,
  			       E1 );
    pDau2.Set( -pDau1.x(), -pDau1.y(), pz2, E2 );

    // Rotate & Boost back to original frame
    Rotation( UpSpect, -angle, k );
    Rotation( pDau1,  -angle, k );
    Rotation( pDau2,  -angle, k );

    UpSpect = BoostBack( beta, UpSpect);
    pDau1   = BoostBack( beta, pDau1);
    pDau2   = BoostBack( beta, pDau2);
  }
  else if(DATA.shower_kernel==1) {
    double rkt = std::sqrt(mar2*y*z*(1.-z));
    FourVector kt1 = Cross(pSplit, pSpect);
    if (kt1.p3abs() < 1.0e-8) kt1 = Cross(pSplit, FourVector(1.,0.,0.,0.));
    kt1 *= rkt*std::cos(pt_angle)/kt1.p3abs();
    FourVector kt2cms = Cross(BoostForCS(pSplit+pSpect, pSplit), kt1);
    kt2cms *= rkt*std::sin(pt_angle)/kt2cms.p3abs();
    FourVector kt2 = BoostBackForCS(pSplit+pSpect, kt2cms);
    pDau1   = pSplit*z      + pSpect*(1.-z)*y + kt1 + kt2; //pi
    pDau2   = pSplit*(1.-z) + pSpect*z*y      - kt1 - kt2; //pj
    UpSpect = pSpect*(1.-y); //pk
  }

  // Make colours of daughters
  int col1[2] = {0};
  int col2[2] = {0};
  int dau1_id = kernels[Kernel]->flav(1); //"original" particle
  int dau2_id = kernels[Kernel]->flav(2); //"new" emitted paricle
  MakeColours( Split, Spect, dau1_id, col1, col2 );

  // Add daughters to parton list
  FourVector x;
  Parton daughter1( Parton( dau1_id, 51, pDau1, x ) ); //Status: 51, active shower particles
  Parton daughter2( Parton( dau2_id, 51, pDau2, x ) );
  daughter1.set_mom1( Split ), daughter1.set_mom2( 0 );
  daughter2.set_mom1( Split ), daughter2.set_mom2( 0 );
  daughter1.set_cols( col1 );
  daughter2.set_cols( col2 );

  parton_list.push_back ( daughter1 );
  int Dau1 = int(parton_list.size()) - 1;
  parton_list.push_back ( daughter2 );
  int Dau2 = int(parton_list.size()) - 1;

  // Add recoiler to parton list
  Parton recoiler = parton_list[Spect];
  recoiler.reset_momentum( UpSpect );
  recoiler.set_mom1( Spect ), recoiler.set_mom2( Spect );
  recoiler.set_stat( 52 );
  parton_list.push_back( recoiler );
  int Rec = int(parton_list.size()) - 1;

  // Update splitter + recoiler
  parton_list[Split].set_stat( -parton_list[Split].stat() ); //It decayed: state -> decayed
  parton_list[Spect].set_stat( -parton_list[Spect].stat() ); //It recoiled: state -> recoiler

  parton_list[Split].set_d1( Dau1 );
  parton_list[Split].set_d2( Dau2 );

  parton_list[Spect].set_d1( Rec );
  parton_list[Spect].set_d2( Rec );
}

void Shower::MakeColours( int Split, int Spect, int dau_id, int col1[2], int col2[2]){

  max_colour++;
  int split_id = parton_list[Split].id(); //flavor before split
  int cols[2] = { parton_list[Split].col(), parton_list[Split].acol() };
  int colr[2] = { parton_list[Spect].col(), parton_list[Spect].acol() };

  // If quark
  if ( split_id != 21 ) {
    if ( split_id > 0 ) { //q --> q + g
      col1[0] = max_colour, col1[1] = 0;
      col2[0] = cols[0], col2[1] = max_colour;
    } else { //qbar --> qbar +
      col1[0] = 0, col1[1] = max_colour;
      col2[0] = max_colour, col2[1] = cols[1];
    }
  }
  // If gluon
  else if (split_id == 21) {
    if ( dau_id == 21 ) { //g --> g + g
      if ( cols[0] == colr[1] ) {
      	// If recoiler is singlet partner, random assignation
      	if ( cols[1] == colr[0] && dis(gen) > 0.5 ) {
          col1[0] = cols[0], col1[1] = max_colour;
      	  col2[0] = max_colour, col2[1] = cols[1];
          return;
      	}
        col1[0] = max_colour, col1[1] = cols[1];
      	col2[0] = cols[0], col2[1] = max_colour;
      } else {
        col1[0] = cols[0], col1[1] = max_colour;
	      col2[0] = max_colour, col2[1] = cols[1];
      }
    } else {
      if ( dau_id > 0 ) {
        col1[0] = cols[0], col1[1] = 0;
        col2[0] = 0, col2[1] = cols[1];
      } else {
        col1[0] = 0, col1[1] = cols[1];
        col2[0] = cols[0], col2[1] = 0;
      }
    }
  }

  return;
}

} // end namespace Adkoda
