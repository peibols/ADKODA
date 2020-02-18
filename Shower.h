#ifndef Shower_H
#define Shower_H

#include <vector>
#include <random>

#include "Parton.h"
#include "Kernel.h"
#include "data.h"
#include "InPartons.h"

namespace Adkoda {

class Shower {

  public:

    Shower(const InitData &DATA_in);
    ~Shower() {}

    void init( InPartons inpartons);

    void run();

    bool evolve();

    void print();

    double alpha_s( double t );

    void Update( int Split, int Spect, int Kernel, double mar2, double z, double Q2 );

    void MakeColours( int Split, int Spect, int dau_id, int col1[2], int col2[2]);

    std::vector<Parton> get_parton_list() { return parton_list; }

  protected:

    std::uniform_real_distribution<double> dis;
    std::mt19937 gen;

    std::vector<Kernel*> kernels;
    std::vector<Parton> parton_list;
    int max_colour;
    double pt_min;
    double t_min;
    double t_max;
    double max_alpha_s;

  private:

    const InitData &DATA;
};

} // end namespace Adkoda

#endif // Shower_H
