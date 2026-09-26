#!/bin/sh
set -eu

c++ -O2 -std=gnu++17 -fopenmp 3dFEexp.cpp -o a.out
c++ -O2 -std=c++17 -Wall -Wextra 2d_axisymmetric_fe.cpp \
  -o a_2d_axisymmetric_fe
