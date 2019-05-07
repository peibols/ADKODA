#ifndef PARTON_H
#define PARTON_H

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "Parton.h"

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
  _p=p;
  _x=x;

  _d1=-1000;
  _d2=-1000;
  _mom=-1000;

  _mass=0.; //Assume gluon only for now
  _virt=-1000.;
  _z=-1000.;
  _alpha=-1000.;
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

void Parton::set_p(const FourVector& p)
{
  _p=p;
}

const FourVector Parton::p()
{
  return _p;
}

double Parton::en()
{
  return _p.t();
}

double Parton::pplus()
{
  return 1./sqrt(2.)*(_p.t()+_p.z());
}

double Parton::pminus()
{
  return 1./sqrt(2.)*(_p.t()-_p.z());
}

void Parton::set_x(const FourVector& x)
{
  _x=x;
}

const FourVector Parton::x()
{
  return _x;
}

void Parton::set_mom(int mom)
{
  _mom=mom;
}

int Parton::mom()
{
  return _mom;
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

void Parton::set_virt(double virt)
{
  _virt=virt;
}

double Parton::virt()
{
  return _virt;
}

void Parton::set_z(double z)
{
  _z=z;
}

double Parton::z()
{
  return _z;
}

void Parton::set_alpha(double alpha)
{
  _alpha=alpha;
}

double Parton::alpha()
{
  return _alpha;
}

#endif
