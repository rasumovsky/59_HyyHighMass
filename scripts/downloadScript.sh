#!/bin/bash

taskName=$1
outputDir=$2

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
localSetupRucioClients
voms-proxy-init --voms atlas

python /afs/cern.ch/user/w/wguan/public/toy/toyoutput.py  --outputDir=${outputDir} --taskname=${taskName}
