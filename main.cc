#include <string>
#include <fstream>
#include "BerGEN.h"

using namespace Adkoda;
using namespace Util;

int main(int argc, char **argv) {

  assert(argc==2);

  std::string input_file;
  input_file = *(argv+1);

  std::cout << "Start Program" << std::endl;

  // Read input file
  BerGEN bergen(input_file);

  int nEv=bergen.number_events();

  bergen.init();

  for ( int iEv = 0; iEv < nEv; iEv++) {

    if (iEv % 1000 == 0) std::cout << "Event: " <<  iEv << std::endl;
    bergen.next();

  }
  std::cout << "End Program" << std::endl;

  return 0;

// End program
}
