#ifndef Kernel_H
#define Kernel_H

#include "Util.h"

namespace Adkoda {

class Kernel {

  public:

    Kernel (int f1, int f2, int f3) { flavs[0] = f1, flavs[1] = f2, flavs[2] = f3 ;} // 1 --> 2 + 3
    virtual ~Kernel () {};

    virtual double Value ( double z, double y ) = 0;
    virtual double Estimate ( double z ) = 0;
    virtual double Integral ( double zm, double zp ) = 0;
    virtual double GenerateZ ( double zm, double zp , double R) = 0;

    int flav ( int i_flav ) { return flavs[i_flav]; }

  private:

    int flavs[3];

};

//Altarelli-Parisi Splitting Functions
class Pgg : public Kernel { //g --> g + g
  public:
    using Kernel::Kernel;
    double Value ( double z, double y ) {return CA * ( z / (1.-z) + (1.-z) / z + z * (1.-z) );}
    double Estimate ( double z ) {return CA * ( 1. / z + 1. / (1.-z) );}
    double Integral ( double zm, double zp ) {return CA * std::log( (zp*(1.-zm)) / (zm*(1.-zp)) );}
    double GenerateZ ( double zm, double zp , double R) { return 1. / ( 1. + (1.-zm)/zm * std::pow((zp*(1.-zm)) / (zm*(1.-zp)), -R) );}
};

class Pqq : public Kernel { //q --> q + g
  public:
    using Kernel::Kernel;
    double Value ( double z, double y ) {return CF * (1.+z*z) / (1.-z);} //FIXME this kernel should be fine.
    double Estimate ( double z ) {return CF * 2. / (1.-z);}
    double Integral ( double zm, double zp ) {return CF * 2. * std::log((1.-zm) / (1.-zp));}
    double GenerateZ ( double zm, double zp , double R) {return 1.0 + (zm - 1.0) * std::pow( (1.0 - zp) / (1.0 - zm), R );}
};

class Pgq : public Kernel { //g --> q + qbar
  public:
    using Kernel::Kernel;
    double Value ( double z, double y ) {return TR * (z*z + (1.-z)*(1.-z)); } //FIXME there is a 1/2 difference in AP and CS.
    double Estimate ( double z ) {return TR;}
    double Integral ( double zm, double zp ) {return TR * (zp - zm);}
    double GenerateZ ( double zm, double zp , double R) {return zm + R * (zp - zm); }
};

//Catani-Seymour Splitting Functions
class Pgg_CS : public Kernel {//g --> g + g
  public:
    using Kernel::Kernel;
    double Value (double z, double y) {return CA / 2. * ( 2./(1.-z*(1.-y)) - 2. + z*(1.-z) );} //FIXME extra 1/2, not symmetrized in z?
    double Estimate (double z) {return CA / (1.-z);}
    double Integral (double zm, double zp) {return CA * std::log((1.-zm)/(1.-zp));}
    double GenerateZ (double zm, double zp, double R) {return 1. + (zp-1.) * std::pow((1.-zm)/(1.-zp), R);}
};

class Pqq_CS : public Kernel {//q --> q + g
  using Kernel::Kernel;
  public:
    double Value (double z, double y) {return CF * ( 2./(1.-z*(1.-y)) - (1.+z) );}
    double Estimate (double z) {return CF * 2./(1.-z);}
    double Integral (double zm, double zp) {return CF * 2. * std::log((1.-zm)/(1.-zp));}
    double GenerateZ (double zm, double zp, double R) {return 1. + (zp-1.) * std::pow((1.-zm)/(1.-zp), R);}
};

class Pgq_CS : public Kernel {//g --> q + qbar
  public:
    using Kernel::Kernel;
    double Value (double z, double y) {return TR/2. * (1. - 2.*z*(1.-z));} //FIXME extra 1/2?
    double Estimate (double z) {return TR/2.;}
    double Integral (double zm, double zp) {return TR/2. * (zp-zm);}
    double GenerateZ (double zm, double zp, double R) {return zm + (zp-zm) * R;}
};

//Medium cascade Splitting Functions
class Kgg : public Kernel {//g --> g + g
  public:
    using Kernel::Kernel;
    double Value (double z, double y) {return CA * std::pow(1.+z*(z-1.),2.) * std::sqrt(CA*(z/(1.-z)+1./z)) /z /(1.-z);}
    double Estimate (double z) {return std::pow(CA /z /(1.-z),3./2.);}
    double Integral (double zm, double zp) {return 2. * CA * (1.-2.*zm) * std::sqrt(CA/(zm-zm*zm));} // Includes the 1/2
    double GenerateZ (double zm, double zp, double R) {return (R*(2.-4.*zm) + std::sqrt(1.+4.*R*(R-1.)*std::pow(1.-2.*zm,2.)) + 2.*zm -1.)/(2.*std::sqrt(1.+4.*R*(R-1.)*std::pow(1.-2.*zm,2.)));} 
};

} // end namespace Adkoda

#endif // Kernel_H
