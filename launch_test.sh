#!/bin/zsh

cmake -B build -DCMAKE_BUILD_TYPE=Release
make -C build
./build/main

gnuplot plot_energy.gp
