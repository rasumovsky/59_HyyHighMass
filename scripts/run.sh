#!/bin/bash

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

./GlobalP0Toys_Faster settings_Scalar.cfg none $1 $2 0

cp ./Scalar_Mar9/GlobalP0Toys/single_files/toy_mu0_$1.root ./
