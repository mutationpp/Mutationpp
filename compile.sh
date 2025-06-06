#!/bin/bash
rm -r build; mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=/Users/raghava/Desktop/CHANL_HPC/Mutationpp/install ..
make -j 20 install
