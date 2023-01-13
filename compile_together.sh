#! /bin/bash

g++ -O2 -ffast-math -march=native -Wall -Wextra -Wno-unused-parameter \
-o $1 $1.cpp \
overlap.cpp  \
effective_current.cpp \
cylinder_modes.cpp \
integration.cpp \
-lgsl -lgslcblas -lm -lfmt -std=c++17
