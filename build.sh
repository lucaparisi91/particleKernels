#! /usr/bin/sh
mpic++ -O3  -std=c++17 -Wfatal-errors randomKernel.cpp -lmkl_core  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -fopenmp   main.cpp  -o main
