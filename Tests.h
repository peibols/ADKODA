#ifndef Tests_H
#define Tests_H

#include "Shower.h"
#include "BerGEN.h"
#include "InPartons.h"
#include "data.h"

namespace Adkoda {

class Tests {

  public:

    Tests();
    ~Tests() {}

    void run(std::vector<Parton> parton_list, double event_weight, double event_xsec, int iEv);
    void close(int nEv);

  protected:

    void Test_Weights(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

    void Test_PrintFinalPartons(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

    void Test_EnergyMomentumConservation(std::vector<Parton> parton_list);

    void Test_FirstSplit(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

    void Test_PrintLundPlane_history(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight);

    void Test_PrintLundPlane_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile, std::ofstream &outfile_weight);

    void Test_JetMass_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

    void Test_IteratedSoftDrop_FJ(std::vector<Parton> parton_list, double event_xsec, double event_weight, std::ofstream &outfile);

    std::vector<int> ListFinalDescendants(std::vector<Parton> parton_list, int ip);

    std::vector<int> ListSplits(std::vector<Parton> parton_list, int ip);

    std::vector<int> CountRealSplits(std::vector<std::vector<int>> split_arr);

    int CountSplitsInMedium(std::vector<Parton> parton_list, std::vector<int> splits);

    void CascadeDist(std::vector<Parton> parton_list);

  private:

    std::ofstream outfile_test_weights;
    std::ofstream outfile_test_FinalPartons;
    std::ofstream outfile_test_FirstSplit;
    std::ofstream outfile_test_LundPlane_history;
    std::ofstream outfile_test_LundPlane_history_weight;
    std::ofstream outfile_test_LundPlane_FJ;
    std::ofstream outfile_test_LundPlane_FJ_weight;
    std::ofstream outfile_test_JetMass_FJ;
    std::ofstream outfile_test_IteratedSoftDrop_FJ;

    int max_nx;
    double x_bin;
    double dx_hist[100][6];

    int max_kx;
    double kx_bin;
    double kx_hist[101][6];

    int ntau_cuts;
    double tau_cuts[6];

};

} // end namespace Adkoda

#endif // Tests_H
