#!/bin/bash

rm -rf build
mkdir build
cd build
cmake ..
make

cp ../test.input src
cp ../Grids/overset/part1.grd src
cp ../Grids/overset/part2.grd src 
cp ../Grids/overset/background.grd src
cd src
mpirun -np 4 ./myAPP