#ifndef PARTON_H
#define PARTON_H

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <assert.h>
#include "Parton.h"

Parton::~Parton(){}

Parton::Parton (const Parton& srp) : PseudoJet (srp)
{
  _id=srp._id;
  _stat=srp._stat;
  _p=srp._p;
  _x=srp._x;
}

Parton::Parton (int id, int stat, const FourVector& p, const FourVector& x)
{
  _id=id;
  _stat=stat;
  _p=p;
  _x=x;
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

void Parton::set_x(const FourVector& x)
{
  _x=x;
}

const FourVector Parton::x()
{
  return _x;
}

#endif
