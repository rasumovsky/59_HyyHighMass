#!/bin/bash

# For now, use a constant width:
directory=$1

# Loop over mass values:
for f in $directory/*
do
    tar zxvf $f
    rm $f
done
