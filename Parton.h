#include "utils/FourVector.h"
#include <stdio.h>
#include <math.h>
#include "utils/fjcore.hh"

#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>

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
    Parton (const Parton& srp);
	  
    virtual ~Parton();

  protected:

    void set_id(int id);
    int id();
    
    void set_stat(int stat);
    int stat();

    void set_p(const FourVector& p);
    const FourVector p();

    void set_x(const FourVector& x);
    const FourVector x();
  
  private:

    FourVector _p;
    FourVector _x;
    int _id;
    int _stat;
};
