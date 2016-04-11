#!/bin/bash

# Set up GCC and ROOT:
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

# For now, use a constant width:
width=10
toyOption="ToyForScan_ToyImportSamp"

# Loop over mass values:
for mass in 730 740 750 760 770
do
    # Only use single cross-section:
    xs=20000
    echo "mass $mass and width $width and cross-section $xs"
    ./GenericToys settings_Graviton.cfg $toyOption $1 $2 0 $mass $width $xs
done

# Then copy all output files to a single .tar file:
tar zcf toy_$1.root.tar Graviton_Mar30/GenericToys/single_files/*
