#!/bin/bash

g++ -g -std=c++11 -Wall $1.cc Parton.cc utils/fjcore.cc sudakovs.cc zfraction.cc -o $1 -I/usr/local/include/ -L /usr/local/lib/ -lgsl -lgslcblas
