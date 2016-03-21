# A search for resonances in the diphoton final state with ATLAS

### Introduction
This package implements statistical tools for a search for diphoton decays of 
scalar and spin-2 particles into the diphoton final state using the ATLAS 
detector at the LHC. 

The code has been structured so that all general analysis settings are specified
in the configuration files residing in the `data/` directory. 

If you have any questions about the code, please contact the author, who will be
happy to assist you in getting the code running.

##### Input files

The analysis requires RooWorkspace objects stored in `.root` files as inputs to 
the analysis software. 

##### Config files

The configuration files for the analysis are included in the `data/` directory.
The baseline config files are currently:
 - settings_Graviton_TightIso.cfg (graviton analysis with tight isolation)
 - settings_Scalar.cfg (scalar analysis)

The config files allow the user to specify everything that might change, such as
analysis luminosity. Only developers should need to modify the code.

##### Compiling and running the package

Several executables should be compiled prior to running. The package currently
uses GCC and ROOT versions hosted on cvmfs, in order to facilitate grid and 
cluster jobs:
```
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/x86_64-slc6-gcc49-opt/setup.sh
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.14-x86_64-slc6-gcc49-opt/bin/thisroot.sh
```

Alternatively, one can run the setup script:
```
source scripts/package_setup.sh
```

To compile the executables:
```
make bin/HMMaster  
make bin/GlobalP0Toys_Faster  
make bin/GlobalP0Analysis  
make bin/LocalP0Toys  
make bin/LocalP0Analysis  
make bin/TossToys_NoFit
```

Every program in this package can be run via the "HMMaster" main method. The 
user only needs to specify the step of the analysis to run and the config file
with all of the analysis settings of interest. For instance, to perform the 
analysis of the global p0 for the scalar analysis, simply run the command below:
```
./bin/DHMaster GlobalP0Toys data/settings_Scalar.cfg
```

The analysis steps are listed and explained below.
 - GlobalP0Toys: create a pseudo-experiment ensemble to get the global Z value.
 - GlobalP0Analysis: analyze the global toy ensemble to get the global z value.
 - LocalP0Analysis: analyze pseudo-experiments to cross-check asymptotic qmu,q0.

