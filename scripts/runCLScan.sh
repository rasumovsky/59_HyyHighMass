#!/bin/bash

# Set up GCC and ROOT:
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

# For now, use a constant width:
width=10
toyOption="ToyForScan"

# Loop over mass values:
#for mass in 500 1000 1500 2000 2500
for mass in 500
do
    # Loop over cross-sections (for CL/p0 scan):
    #for xs in 100 200 400 800 1600 3200 6400 12800 25600 51200
    for xs in 6400
    do
	echo "mass $mass and width $width and cross-section $xs"
	./GenericToys settings_Graviton.cfg $toyOption $1 $2 0 $mass $width $xs
	./GenericToys settings_Graviton.cfg $toyOption $1 $2 1 $mass $width $xs
    done
done

# Then copy all output files to a single .tar file:
tar zcf toy_$1.root.tar Graviton_Mar30/GenericToys/single_files/*
