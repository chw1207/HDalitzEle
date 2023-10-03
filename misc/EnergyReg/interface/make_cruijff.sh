#!/usr/bin/env bash 

g++ RooCruijff.cpp -o libCruijff.so -shared -fPIC -I$CONDA_PREFIX/include -std=c++17 -Wall -O3 $(root-config --glibs --cflags --libs) -lRooFitCore -lRooFit 

