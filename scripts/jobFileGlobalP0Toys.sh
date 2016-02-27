#!/bin/bash

if [[ $# -lt 7 ]]; then
    echo "USAGE: toy_jobfile.sh <jobname> <configfile> <input_file> <exe_name> <option> <seed> <toysperjob>"
    
else
    jobname=$1
    configfile=$2
    input_file=$3
    exe_name=$4
    option=$5
    seed=$6
    toysperjob=$7

    date
    echo
    
    echo $jobname $configfile $input_file $exe_name $option $seed $toysperjob
    
    out="${jobname}_${seed}"
    output_dir="/afs/cern.ch/work/a/ahard/files_HighMass/FullAnalysis/${jobname}/GlobalP0Toys"
    
    # setup ROOT:
    #source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh 
    #source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin/thisroot.sh
    
    # setup the HSG7 patched ROOT:
    source /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6/setup.sh 
    source /afs/cern.ch/atlas/project/HSG7/root/root_v5-34-32/x86_64-slc6-gcc48/bin/thisroot.sh
    
    export PATH=$ROOTSYS/bin:$PATH
    export LD_LIBRARY_PATH=$ROOTSYS/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/project/eos/installation/pro/lib64/
    
    # setup GCC:
    source  /afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
    
    # Make output directories:    
    mkdir -vp ${output_dir}/single_files
    mkdir -vp ${output_dir}/log
    mkdir -vp ${output_dir}/err
    
####################
# copying necessary inputs into working dir
    cp $input_file .
    tar zxvf Cocoon.tar

####################    
# Go!
    echo "Printing directory contents before running."
    ls
    
    ./bin/${exe_name} ${configfile} ${option} ${seed} ${toysperjob} 0 1> ${out}_mu0.log 2>${out}_mu0.err;
#    ./bin/${exe_name} ${configfile} ${option} ${seed} ${toysperjob} 1 1> ${out}_mu1.log 2>${out}_mu1.err;
    
    mv *.log ${output_dir}/log/
    mv *.err ${output_dir}/err/
    rm * -rf
    
####################
# End
    echo "----------All done------------"
    date
    
fi
