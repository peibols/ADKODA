#ifndef Util_H
#define Util_H

#ifndef NC
#define NC (3.0)
#endif

#ifndef TR
#define TR (0.5)
#endif

#ifndef CA
#define CA (3.0)
#endif

#ifndef CF
#define CF ((NC*NC-1.)/(2.*NC))
#endif

#ifndef Lambda
#define Lambda (0.2)
#endif

#ifndef MZ2
#define MZ2 (91.1876*91.1876)
#endif

#ifndef GZ2
#define GZ2 (2.4952*2.4952)
#endif

#ifndef Alpha_QED
#define Alpha_QED (1./128.802)
#endif

#ifndef sin2tw
#define sin2tw (0.22293)
#endif

#include "utils/FourVector.h"

namespace Util {

  double Delta(FourVector v1, FourVector v2);
  double m(const FourVector& v1, const FourVector& v2);
  double m2(const FourVector& v1, const FourVector& v2);
  FourVector Cross(const FourVector& v1, const FourVector& v2);
  FourVector BoostForCS(FourVector p, FourVector q);
  FourVector BoostBackForCS(FourVector p, FourVector q);
  FourVector Boost(double v[3], FourVector p);
  FourVector BoostBack(double v[3], FourVector p);

  void VelCOM( FourVector p1, FourVector p2, double v[3]);

  void rotate_parton( double prog[4], double u_z[3]);
  void Rotation( FourVector &vector, double angle, double axis[3]);
  void AlignWithZ( FourVector &vec, double &angle, double k[3]);

  std::string StringFind4(std::string file_name, std::string str_in);

}

#endif // Util_H
