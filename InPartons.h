#ifndef InPartons_H
#define InPartons_H

#include <vector>
#include <random>

#include "Parton.h"
#include "data.h"

namespace Adkoda {

class InPartons {

  public:

    InPartons(const InitData &DATA_in) : DATA(DATA_in) {}
    ~InPartons() {}

   double event_weight;
   std::vector<Parton> PartonList();
   double ME2(int id, double s, double t);

 protected:

   std::mt19937 gen;
   std::uniform_real_distribution<double> dis;
   std::uniform_int_distribution<int> dis_int;

  private:

    const InitData &DATA;
};

} // end namespace Adkoda

#endif // InPartons_H
