#ifndef BerGEN_H
#define BerGEN_H

#include "Shower.h"
#include "InPartons.h"
#include "data.h"

namespace Adkoda {

class BerGEN {

  private:

    Shower    *shower    =	nullptr;
    InPartons *inpartons = 	nullptr;
    InitData  DATA;

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

};

} // end namespace Adkoda

#endif // BerGEN_H
