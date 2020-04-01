#include "Shower.h"
#include "BerGEN.h"
#include "InPartons.h"
#include "data.h"

namespace Adkoda {

void Test_Weights(std::vector<Parton> parton_list, double event_weight, std::ofstream &outfile);

void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_weight, std::ofstream &outfile);

void Test_EnergyMomentumConservation(std::vector<Parton> parton_list);

void Test_PrintLundPlane(std::vector<Parton> parton_list, double event_weight, std::ofstream &outfile);

}
