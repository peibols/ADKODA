# to compile: ./compiler.sh main
#!/bin/bash
g++ -g -std=c++11 -Wall BerGEN.cc utils/fjcore.cc Util.cc Shower.cc Parton.cc InPartons.cc Evolution.cc $1.cc -o $1
