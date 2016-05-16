#!/bin/bash

# Input directory containing all MxAOD folders:
directoryIn=$1
directoryOut=$2
# Loop over mass values:
for d in $directoryIn/*
do
    echo "d=" $d
    for f in $d/* 
    do
	echo "f=" $f
	cp $f $directoryOut/
    done
done

hadd -f combinedMxAOD.root $directoryOut/*.root