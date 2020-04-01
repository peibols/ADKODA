#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "Parton.h"

namespace Adkoda {

Parton::~Parton(){}

/*
Parton::Parton (const Parton& srp) : PseudoJet (srp)
{
  _id=srp._id;
  _stat=srp._stat;
  _p=srp._p;
  _x=srp._x;

  _d1=srp._d1;
  _d1=srp._d2;
  _mom=srp._mom;

  _mass=srp._mass;
  _virt=srp._virt;
  _z=srp._z;
}
*/

Parton::Parton (int id, int stat, const FourVector& p, const FourVector& x)
{
  _id=id;
  _stat=stat;
  reset_momentum (p);
  _x=x;

  _d1=0;
  _d2=0;
  _mom1=0;
  _mom2=0;

  _mass=0.; //FIXME Assume all partons massless for the moment

  _col=0;
  _acol=0;
}

void Parton::display()
{
  cout << _id << "\t "
	<< _stat << "\t "
	<< _mom1 << "\t " << _mom2 << "\t "
	<< _d1 << "\t " << _d2 << "\t "
	<< _col << "\t " << _acol << "\t "
	<< p().x() << "\t "
	<< p().y() << "\t "
	<< p().z() << "\t "
	<< p().t() << "\t "
	<< _mass << endl;
}

bool Parton::ColourConnected(Parton& p)
{
  if ( (_col < 0 || _acol < 0) || (_col == 0 && _acol == 0) ) { std::cout << "This parton has not assigned colours! \n"; exit(0); }
  if ( (p.col() < 0 || p.acol() < 0) || (p.col() == 0 && p.acol() == 0) ) { std::cout << "Recoil parton has not assigned colours! \n"; exit(0); }
  if ( (p.col() == _acol && _acol != 0) || (p.acol() == _col && _col != 0) ) return 1;
  else return 0;
}

void Parton::set_id(int id)
{
  _id=id;
}

int Parton::id()
{
  return _id;
}

void Parton::set_stat(int stat)
{
  _stat=stat;
}

int Parton::stat()
{
  return _stat;
}

FourVector Parton::p()
{
  fjcore::PseudoJet psj = GetPseudoJet();
  FourVector pmu (psj.px(), psj.py(), psj.pz(), psj.e());
  return pmu;
}

void Parton::set_x(const FourVector& x)
{
  _x=x;
}

const FourVector Parton::x()
{
  return _x;
}

void Parton::set_d1(int d1)
{
  _d1=d1;
}

int Parton::d1()
{
  return _d1;
}

void Parton::set_d2(int d2)
{
  _d2=d2;
}

int Parton::d2()
{
  return _d2;
}

void Parton::set_mass(double mass)
{
  _mass=mass;
}

double Parton::mass()
{
  return _mass;
}

} //end namespace Adkoda
