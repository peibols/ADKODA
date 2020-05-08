#include "Cascade.h"
#include "Util.h"

namespace Adkoda {

double Cascade::BroadInt()
{

  double alphabar_med = alphas_med;
  double tstar = 1./alphabar_med/std::sqrt(qhat/pplus);
  double n = 1.;
  return 8.*M_PI*alphas_med*alphas_med*n*NC*tstar/qmin/qmin; 

}

double Cascade::GenerateQ(double R)
{

 return qmin/std::sqrt(1.-R);

}

} // end namespace Adkoda
