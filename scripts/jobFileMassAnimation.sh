#!/bin/bash

if [[ $# -lt 7 ]]; then
    echo "USAGE: toy_jobfile.sh <jobname> <configfile> <input_file> <exe_name> <option> <minFrame> <maxFrame>"
    
else
    jobname=$1
    configfile=$2
    input_file=$3
    exe_name=$4
    option=$5
    minFrame=$6
    maxFrame=$7

    date
    echo
    
    echo $jobname $configfile $input_file $exe_name $option $minFrame $maxFrame
    
    out="min${minFrame}"
    output_dir="/afs/cern.ch/work/a/ahard/files_HighMass/FullAnalysis/${jobname}/MassAnimation"
    
    # HSG7 patched ROOT on cvmfs:
    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh
    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/project/eos/installation/pro/lib64/
    
    # setup GCC:
    #source  /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
    
    # Make output directories:    
    mkdir -vp ${output_dir}/single_files
    mkdir -vp ${output_dir}/log
    mkdir -vp ${output_dir}/err
    
####################
# copying necessary inputs into working dir
    cp $input_file .
    tar zxvf Cocoon.tar
    mv ForJob/* .

####################    
# Go!
    echo "Printing directory contents before running."
    ls
    
    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh

    source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh

    #./bin/${exe_name} ${configfile} ${option} ${minFrame} ${maxFrame} 1> ${out}.log 2>${out}.err;
    ./${exe_name} ${configfile} ${option} ${minFrame} ${maxFrame} 1> ${out}.log 2>${out}.err;
    
    mv *.log ${output_dir}/log/
    mv *.err ${output_dir}/err/
    rm * -rf
    
####################
# End
    echo "----------All done------------"
    date
    
fi
