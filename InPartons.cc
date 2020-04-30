#include <vector>

#include "InPartons.h"
#include "Util.h"

using namespace Util;

namespace Adkoda {

std::vector<Parton> InPartons::PartonList() {
  int id1 = -1000;
  int id2 = -1000;
  int cols1[2], cols2[2];

  // Random number generator //FIXME where to initialize and put the random generator?
  std::random_device rd; // obtain a random number from hardware
  std::mt19937 gen(rd()); // seed the generator
  dis = std::uniform_real_distribution<double> { 0.0, 1.0 }; // define the range
  dis_int = std::uniform_int_distribution<> { 1, 5 };

  FourVector p0, p1, p2, x;
  std::vector<Parton> hard_list;

  double ecms = 0.;
  if (DATA.parton_gun==1) ecms = DATA.pt_max;
  else if (DATA.parton_gun==0) ecms = 91.1876;
  else std::cout << "ERROR: parton_gun = 0,1" << std::endl;
  p0.Set( 0., 0., 0., ecms );
  Parton hard_parton_system( Parton(90,-11,p0,x) );
  hard_parton_system.set_mass(ecms);
  hard_list.push_back(hard_parton_system);

  if (DATA.parton_gun == 1) { // Back-to-back di-parton gun random angle in the transverse plane

    event_weight = 1.;
    event_xsec   = 1.;

    double phi = 2.*M_PI*dis(gen);
    phi = 0.;
    p1.Set(std::cos(phi)*ecms/2., std::sin(phi)*ecms/2., 0., ecms/2.);
    p2.Set(-p1.x(), -p1.y(), -p1.z(), p1.t());

    // Assign id and color for the hard scattering.
    if (DATA.hard_partons == 0) { // gg
      cols1[0]=101, cols1[1]=102;
      cols2[0]=102, cols2[1]=101;
      id1=21, id2=21;
    }
    else if (DATA.hard_partons == 1) { // qqbar
      cols1[0]=101, cols1[1]=0;
      cols2[0]=0, cols2[1]=101;
      id1=dis_int(gen), id2=id1;
    }
    else if (DATA.hard_partons == 2) { //nonphysical
      cols1[0]=101, cols1[1]=0;
      cols2[0]=102, cols2[1]=101;
      id1=dis_int(gen), id2=21;
    }
    else std::cout << "ERROR: hard_partons = 0-2" << std::endl;

    Parton hard_parton_a( Parton(id1,-12,p1,x) ); //stat:12 incoming beam
    Parton hard_parton_b( Parton(id2,-12,p2,x) );
    hard_parton_a.set_d1(3); //FIXME I set it by hand.
    hard_parton_b.set_d1(3);
    hard_list.push_back(hard_parton_a);
    hard_list.push_back(hard_parton_b);
    Parton hard_parton_c( Parton(90,-22,p1+p1,x) ); //id:0 (System), stat:22 intermediate
    hard_parton_c.set_mass(ecms);
    hard_parton_c.set_mom1(1); //FIXME set it automatically
    hard_parton_c.set_mom2(2);
    hard_parton_c.set_d1(4);
    hard_parton_c.set_d2(5);
    hard_list.push_back(hard_parton_c);
  }
  else if (DATA.parton_gun == 0) { // LO: e-e+ --> jj
    double ct = 2.*dis(gen)-1.;
    double st = std::sqrt(1.-ct*ct);
    double phi = 2.*M_PI*dis(gen);
    p1.Set( st*std::cos(phi)*ecms/2., st*std::sin(phi)*ecms/2., ct*ecms/2., ecms/2. );
    p2.Set( -p1.x(), -p1.y(), -p1.z(), p1.t() );
    FourVector pa ( 0., 0., ecms/2., ecms/2. );
    FourVector pb ( 0., 0., -ecms/2., ecms/2. );
    cols1[0]=101, cols1[1]=0;
    cols2[0]=0,   cols2[1]=101;
    id1=dis_int(gen), id2=-id1;
    double s = (pa+pb)*(pa+pb);
    double t = (pa-p1)*(pa-p1);
    double lome = ME2(id1, s, t);
    double dxs  = 5.*lome*3.89379656e8/(8.*M_PI)/(2.*std::pow(ecms, 2.)); //xsec is in pb
    event_xsec   = dxs;
    event_weight = 1.;

    Parton hard_parton_a( Parton(-11,-12,pa,x) ); //id:11 (e), stat:12 incoming beam
    Parton hard_parton_b( Parton(11,-12,pb,x) );
    hard_parton_a.set_d1(3); //FIXME I set it by hand.
    hard_parton_b.set_d1(3);
    hard_list.push_back(hard_parton_a);
    hard_list.push_back(hard_parton_b);
    Parton hard_parton_c( Parton(23,-22,pa+pb,x) ); //id:22 (Z0), stat:22 intermediate
    hard_parton_c.set_mass(91.1876);
    hard_parton_c.set_mom1(1); //FIXME set it automatically
    hard_parton_c.set_mom2(2);
    hard_parton_c.set_d1(4);
    hard_parton_c.set_d2(5);
    hard_list.push_back(hard_parton_c);
  }
  else std::cout << "ERROR: parton_gun = 0,1" << std::endl;

  int stat1 = 23; // Primary, final state
  int stat2 = 23;

  Parton hard_parton1( Parton(id1,stat1,p1,x) );
  Parton hard_parton2( Parton(id2,stat2,p2,x) );
  hard_parton1.set_mom1(hard_list.size()-1);
  hard_parton2.set_mom1(hard_list.size()-1);
  hard_parton1.set_cols(cols1);
  hard_parton2.set_cols(cols2);
  hard_parton1.set_scale(ecms);
  hard_parton2.set_scale(ecms);

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
