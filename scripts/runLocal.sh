#!/bin/bash

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

./LocalP0Toys settings_Scalar.cfg none $1 $2 0
./LocalP0Toys settings_Scalar.cfg none $1 $2 1

tar zcf toy_$1.root.tar Scalar_Mar9/GlobalP0Toys/single_files/*
