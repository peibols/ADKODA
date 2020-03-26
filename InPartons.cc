#include <vector>

#include "InPartons.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

std::vector<Parton> InPartons::PartonList() {
  int id1, id2;
  int cols1[2], cols2[2];

  FourVector p1, p2, x;
  std::vector<Parton> hard_list;

  if (DATA.parton_gun == 1) { // Back-to-back di-parton gun
    double px = DATA.pt_max/2.;
    double py = 0.;
    double pz = 0.;
    double en = px; // FIXME Assumed massless partons for now.

    p1.Set( en, py, pz, en );
    p2.Set( -en, -py, -pz, en );

    // Assign id and color for the hard scattering.
    if (DATA.hard_partons == 0) { // gg
      cols1[0]=101, cols1[1]=102;
      cols2[0]=102, cols2[1]=101;
      id1=21, id2=21;
    }
    else if (DATA.hard_partons == 1) { // qqbar
      cols1[0]=101, cols1[1]=0;
      cols2[0]=0, cols2[1]=101;
      id1=1, id2=-1;
    }
    else if (DATA.hard_partons == 2) {
      cols1[0]=101, cols1[1]=0;
      cols2[0]=102, cols2[1]=101;
      id1=1, id2=21;
    }
  }
  else if (DATA.parton_gun == 0) { // LO: e-e+ --> jj
    double ecms = DATA.pt_max;
    double ct = 2.*dis(gen)-1.;
    double st = std::sqrt(1.-ct*ct);
    double phi = 2.*M_PI*dis(gen);
    p1.Set( st*std::cos(phi)*ecms/2., st*std::sin(phi)*ecms/2., ct*ecms/2., ecms/2. );
    p2.Set( -p1.x(), -p1.y(), -p1.z(), p1.t() );
    FourVector pa ( 0., 0., ecms/2., ecms/2. );
    FourVector pb ( 0., 0., -ecms/2., ecms/2. );
    cols1[0]=101, cols1[1]=0;
    cols2[0]=0,   cols2[1]=101;
    id1=1, id2=-1; //FIXME random flavor
    double s = (pa+pb)*(pa+pb);
    double t = (pa-p1)*(pa-p1);
    double lome = ME2(id1, s, t);
    double dxs = 5.*lome*3.89379656e8/(8.*M_PI)/(2.*std::pow(ecms, 2.)); //FIXME where to put the event weight?

    Parton hard_parton_a( Parton(-11,-12,pa,x) ); //FIXME assign beam particles and daughters in the right way.
    Parton hard_parton_b( Parton(11,-12,pb,x) );

    hard_list.push_back(hard_parton_a);
    hard_list.push_back(hard_parton_b);
  }

  int stat1 = 23; // Primary, final state
  int stat2 = 23;

  Parton hard_parton1( Parton(id1,stat1,p1,x) );
  Parton hard_parton2( Parton(id2,stat2,p2,x) );
  hard_parton1.set_mom1(-3), hard_parton1.set_mom2(-4); // Mothers -3,-4 = primary
  hard_parton2.set_mom1(-3), hard_parton2.set_mom2(-4);
  hard_parton1.set_cols(cols1);
  hard_parton2.set_cols(cols2);

  hard_list.push_back(hard_parton1);
  hard_list.push_back(hard_parton2);

  return hard_list;
}

double InPartons::ME2(int fl, double s, double t) {
  double qe = -1.;
  double ae = -0.5;
  double ve = ae - 2.*qe*sin2tw;
  double qf, af;
  if (fl == 2 || fl == 4) { qf = 2./3.; af = 0.5; }
  else { qf = -1./3.; af = -0.5; }
  double vf    = af - 2.*qf*sin2tw;
  double kappa = 1./(4.*sin2tw*(1.-sin2tw));
  double chi1  = kappa * s * (s-MZ2)/(std::pow(s-MZ2, 2.) + GZ2*MZ2);
  double chi2  = std::pow(kappa * s, 2.)/(std::pow(s-MZ2, 2.) + GZ2*MZ2);
  double term1 = (1.+std::pow(1.+2.*t/s, 2.))*(std::pow(qf*qe, 2.) + 2.*(qf*qe*vf*ve)*chi1 + (ae*ae+ve*ve)*(af*af+vf*vf)*chi2);
  double term2 = (1.+2.*t/s)*(4.*qe*qf*ae*af*chi1 + 8.*ae*ve*af*vf*chi2);
  return std::pow(4.*M_PI*Alpha_QED, 2.) * 3. * (term1 + term2);
}

} // end namespace Adkoda
