#ifndef Parton_H
#define Parton_H

#include "utils/FourVector.h"
#include <stdio.h>
#include <math.h>
#include "utils/fjcore.hh"

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace Adkoda {

using std::ostream;

class Parton : protected fjcore::PseudoJet
{
  friend class fjcore::PseudoJet;

  public:

    inline void reset_momentum( const double px, const double py, const double pz, const double e ){
      fjcore::PseudoJet::reset_momentum ( px, py, pz, e );
    }

    inline void reset_momentum( const FourVector& p ){
      fjcore::PseudoJet::reset_momentum ( p.x(), p.y(), p.z(), p.t() );
    }

    fjcore::PseudoJet GetPseudoJet() const{
      return PseudoJet ( *this );
    }

    using fjcore::PseudoJet::px;
    using fjcore::PseudoJet::py;
    using fjcore::PseudoJet::pz;
    using fjcore::PseudoJet::e;
    using fjcore::PseudoJet::E;

    using fjcore::PseudoJet::phi;
    using fjcore::PseudoJet::phi_std;
    using fjcore::PseudoJet::phi_02pi;
    using fjcore::PseudoJet::rap;
    using fjcore::PseudoJet::rapidity;
    using fjcore::PseudoJet::pseudorapidity;
    using fjcore::PseudoJet::eta;
    using fjcore::PseudoJet::pt2;
    using fjcore::PseudoJet::pt;
    using fjcore::PseudoJet::perp2;
    using fjcore::PseudoJet::perp;
    using fjcore::PseudoJet::kt2;

    using fjcore::PseudoJet::modp2;
    using fjcore::PseudoJet::modp;
    using fjcore::PseudoJet::Et;
    using fjcore::PseudoJet::Et2;

    using fjcore::PseudoJet::kt_distance;
    using fjcore::PseudoJet::plain_distance;
    using fjcore::PseudoJet::squared_distance;
    using fjcore::PseudoJet::delta_R;
    using fjcore::PseudoJet::delta_phi_to;
    using fjcore::PseudoJet::beam_distance;

    using fjcore::PseudoJet::operator*=;
    using fjcore::PseudoJet::operator/=;
    using fjcore::PseudoJet::operator+=;
    using fjcore::PseudoJet::operator-=;

    using fjcore::PseudoJet::user_index;
    using fjcore::PseudoJet::set_user_index;
    using fjcore::PseudoJet::UserInfoBase;
    using fjcore::PseudoJet::InexistentUserInfo;

    using fjcore::PseudoJet::user_info;
    using fjcore::PseudoJet::set_user_info;
    using fjcore::PseudoJet::has_user_info;
    using fjcore::PseudoJet::user_info_ptr;
    using fjcore::PseudoJet::user_info_shared_ptr;

    using fjcore::PseudoJet::description;

  public:

    Parton() : PseudoJet() {};
    Parton (int id, int stat, const FourVector& p, const FourVector& x);
    //Parton (const Parton& srp);

    ~Parton();

    void set_cols( int cols[2] ) { _col = cols[0], _acol = cols[1]; }
    int col() { return _col; }
    int acol() { return _acol; }

    void set_id( int id );
    int id();

    void set_stat(int stat);
    int stat();

    void set_mass(double mass);
    double mass();

    void set_scale(double scale);
    double scale();

    void set_virt(double virt);
    double virt();

    void set_mom1( int mom ) { _mom1 = mom; }
    int mom1() { return _mom1; }
    void set_mom2( int mom ) { _mom2 = mom; }
    int mom2() { return _mom2; }

    void set_d1(int d1);
    int d1();

    void set_d2(int d2);
    int d2();

    FourVector p();

    void set_x(const FourVector& x); //Creation point in the lab frame
    const FourVector x();

    bool ColourConnected(Parton& p);

    void display();

    std::vector<int> motherList() const;

  private:

    FourVector _x;
    int _col, _acol;
    int _id;
    int _stat;
    int _mom1, _mom2, _d1, _d2;
    double _mass, _scale;
};

} // end namespace Adkoda

#endif // Parton_H
