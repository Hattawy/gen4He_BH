#!/bin/bash

rm -rf gen4He_*.lund
rm -rf gen4He_cohdvcs_R.root
./compile.sh gen4He.C
./gen4He -n 50001 -r1 -t0
cd test-alu/
root -b run.cc
