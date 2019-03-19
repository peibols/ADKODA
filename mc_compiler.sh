#!/bin/bash

g++ -std=c++11 -Wall test_mc.cc Parton.cc utils/fjcore.cc sudakovs.cc zfraction.cc -o test_mc -I/usr/local/include/ -L /usr/local/lib/ -lgsl -lgslcblas
