#ifndef Kernel_H
#define Kernel_H

#include "Util.h"

namespace Adkoda {

class Kernel {

  public:

    Kernel (int f1, int f2, int f3) { flavs[0] = f1, flavs[1] = f2, flavs[2] = f3 ;}
    virtual ~Kernel () {};

    virtual double Value ( double z ) = 0;
    virtual double Estimate ( double z ) = 0;
    virtual double Integral ( double zm, double zp ) = 0;
    virtual double GenerateZ ( double zm, double zp , double R) = 0;

    int flav ( int i_flav ) { return flavs[i_flav]; }

  private:

    int flavs[3];

};

class Pgg : public Kernel {

  using Kernel::Kernel;

  public:

    double Value ( double z ) {return CA * ( z / (1.0 - z) + (1.0 - z) / z + z * (1.0 - z) );}
    double Estimate ( double z ) {return CA * ( 1.0 / z + 1.0 / (1.0 - z) );}
    double Integral ( double zm, double zp ) {return CA * std::log( (zp * (1.0 - zm)) / (zm * (1.0 - zp)) );}
    double GenerateZ ( double zm, double zp , double R) { return 1.0 / ( 1.0 + (1.0 - zm)/zm * std::pow((zp * (1.0 - zm)) / (zm * (1.0 - zp)), -R) );}

};

class Pqg : public Kernel {

  public:

    using Kernel::Kernel;

    double Value ( double z ) {return CF * (1.0 + z * z) / (1.0 - z);}
    double Estimate ( double z ) {return CF * 2.0 / (1.0 - z);}
    double Integral ( double zm, double zp ) {return CF * 2.0 * std::log( (1.0 - zm) / (1.0 - zp) );}
    double GenerateZ ( double zm, double zp , double R) {return 1.0 + (zm - 1.0) * std::pow( (1.0 - zp) / (1.0 - zm), R );}

};

class Pqq : public Kernel {

  public:

    using Kernel::Kernel;

    // FIXME should be solved analytically.
    double Value ( double z ) {return 0.5 * ( z * z + ( 1.0 - z ) * ( 1.0 - z ) ); } //This is not proportional with Nf, because it is summed up later for all flavors.
    double Estimate ( double z ) {return 0.5;}
    double Integral ( double zm, double zp ) {return 0.5 * ( zp - zm );}
    double GenerateZ ( double zm, double zp , double R) {return zm + R * ( zp - zm ); }

};

} // end namespace Adkoda

#endif // Kernel_H
