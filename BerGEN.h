#ifndef BerGEN_H
#define BerGEN_H

#include "Shower.h"
#include "Cascade.h"
#include "InPartons.h"
#include "data.h"

namespace Adkoda {

class BerGEN {

  private:

    Shower    *shower    =	nullptr;
    Cascade   *cascade   =	nullptr;
    InPartons *inpartons = 	nullptr;
    InitData  DATA;

  protected:

    double event_weight;
    double event_xsec;
    std::vector<Parton> parton_list;

  public:

    BerGEN(std::string input_file);
    ~BerGEN() {}

    // Read parameters
    InitData read_in_parameters(std::string input_file);

    // Initialization
    void init();

    // Next Event
    void next();

    // Get #Events
    int number_events() { return DATA.number_events; }

    void print();
    double get_event_weight() { return event_weight; }
    double get_event_xsec() { return event_xsec; }
    std::vector<Parton> get_parton_list() { return parton_list; }

};

} // end namespace Adkoda

#endif // BerGEN_H
