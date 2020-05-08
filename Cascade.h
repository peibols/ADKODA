#ifndef Cascade_H
#define Cascade_H

#include <vector>
#include <random>

#include "Parton.h"
#include "Kernel.h"
#include "data.h"

namespace Adkoda {

class Cascade {

  public:

    Cascade(const InitData &DATA_in);
    ~Cascade() {}

    void init(std::vector<Parton> parton_list);

    void run();

    bool evolve(double &start_time, std::vector<Parton> &cascade_list);

    double BroadInt();
    double GenerateQ(double R);

    double tau(double t /*in fm*/);
    double tkin(double tau /*dim less*/); 
    
    std::vector<Parton> get_parton_list() { return parton_list; }

  protected:

    std::mt19937 gen; // seed the generator
    std::uniform_real_distribution<double> dis;
    std::uniform_int_distribution<> dis_int; // define the range: use dis_int

    std::vector<Kernel*> kernels;
    std::vector<Parton> parton_list;

    double L_med;
    double alphas_med;
    double qhat;
    double eps_med;    
    double xmin_med;

    double qmin;

    // Updated by each shower initiator
    double pplus;
    double angle;
    double k[3];

  private:

    const InitData &DATA;
};

} // end namespace Adkoda

#endif // Cascade_H

