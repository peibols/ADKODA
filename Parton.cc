#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "Parton.h"

namespace Adkoda {

Parton::~Parton(){}

Parton::Parton (int id, int stat, const FourVector& p, const FourVector& x) {
  _id=id;
  _stat=stat;
  reset_momentum (p);
  _x=x; //Creation point in the lab frame
  _col=0;
  _acol=0;

  _d1=0;
  _d2=0;
  _mom1=0;
  _mom2=0;
  _mass=0.; //FIXME Assume all partons massless for the moment
  _scale=0.; //Scale pt [GeV] in which was produced

  _xfrac=1.; //Energy fraction wrt cascade initiator

  _is_frozen=0;
}

void Parton::display() {
  cout << _id << "\t "
	<< _stat << "\t "
	<< _mom1 << "\t " << _mom2 << "\t "
	<< _d1 << "\t " << _d2 << "\t "
	<< _col << "\t " << _acol << "\t "
	<< p().x() << "\t " << p().y() << "\t " << p().z() << "\t " << p().t() << "\t "
	<< _mass << "\t "
  << _scale << "\t "
  << x().x() << "\t " << x().y() << "\t " << x().z() << "\t " << x().t()
  << endl;
}

bool Parton::ColourConnected(Parton& p) {
  if ( (_col < 0 || _acol < 0) || (_col == 0 && _acol == 0) ) { std::cout << "This parton has not assigned colours! \n"; exit(0); }
  if ( (p.col() < 0 || p.acol() < 0) || (p.col() == 0 && p.acol() == 0) ) { std::cout << "Recoil parton has not assigned colours! \n"; exit(0); }
  if ( (p.col() == _acol && _acol != 0) || (p.acol() == _col && _col != 0) ) return 1;
  else return 0;
}

void Parton::set_id(int id) { _id=id; }
int Parton::id(){ return _id; }

void Parton::set_stat(int stat) { _stat=stat; }
int Parton::stat() { return _stat; }

FourVector Parton::p() {
  fjcore::PseudoJet psj = GetPseudoJet();
  FourVector pmu (psj.px(), psj.py(), psj.pz(), psj.e());
  return pmu;
}

void Parton::set_x(const FourVector& x) { _x=x; } //Creation point in the lab frame
FourVector Parton::x() { return _x; }

void Parton::set_d1(int d1) { _d1=d1; }
int Parton::d1(){ return _d1; }

void Parton::set_d2(int d2) { _d2=d2; }
int Parton::d2() { return _d2; }

void Parton::set_mass(double mass) { _mass=mass; }
double Parton::mass() { return _mass; }

void Parton::set_scale(double scale) { _scale=scale; } //The pt scale in which was produced.
double Parton::scale() { return _scale; }


std::vector<int> Parton::motherList() const { // Find complete list of mothers.
  // Vector of all the mothers; created empty. Done if no event pointer.
  std::vector<int> motherVec;
  int statusSaveAbs = abs(_stat);
  if  (statusSaveAbs == 11 || statusSaveAbs == 12) ;   // Special cases in the beginning, where the meaning of zero is unclear.
  else if (_mom1 == 0 && _mom2 == 0) motherVec.push_back(0);
  else if (_mom2 == 0 || _mom2 == _mom1) motherVec.push_back(_mom1); // One mother or a carbon copy.
  else if ( (statusSaveAbs >  80 && statusSaveAbs <  90) || (statusSaveAbs > 100 && statusSaveAbs < 107) ) // A range of mothers from string fragmentation.
    for (int iRange = _mom1; iRange <= _mom2; ++iRange) motherVec.push_back(iRange);
  else {    // Two separate mothers.
    motherVec.push_back( std::min(_mom1, _mom1) );
    motherVec.push_back( std::max(_mom1, _mom2) );
  }
  return motherVec;
}

} //end namespace Adkoda
