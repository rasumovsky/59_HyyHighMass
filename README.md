# A search for resonances in the diphoton final state with ATLAS

### Introduction
This package implements statistical tools for a search for diphoton decays of 
scalar and spin-2 particles into the diphoton final state using the ATLAS 
detector at the LHC. 

The code has been structured so that all general analysis settings are specified
in the configuration files residing in the `data/` directory. 

If you have any questions about the code, please contact the author, who will be
happy to assist you in getting the code running.


### Input files

The analysis requires RooWorkspace objects stored in `.root` files as inputs to 
the analysis software. 

### Config files

The configuration files for the analysis are included in the `data/` directory.
The baseline config files are currently:
 - settings_Graviton.cfg (graviton analysis with tight isolation)
 - settings_Scalar.cfg (scalar analysis)

The config files allow the user to specify everything that might change, such as
analysis luminosity. Only developers should need to modify the code.

### Executing the package locally

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
source scripts/makeAll.sh
```

This command compiles the following macros:
 - HMMaster  
 - GlobalP0Toys_Faster  
 - GlobalP0Analysis  
 - GenericToys  
 - MassAnimation
 - LocalP0Analysis  
 - TossToys_NoFit
 - AddDataToWorkspace
 - PlotWS

Every program in this package can be run via the `HMMaster` main method. The 
user only needs to specify the step of the analysis to run and the config file
with all of the analysis settings of interest. For instance, to perform the 
analysis of the global p0 for the scalar analysis, simply run the command below:
```
./bin/DHMaster GlobalP0Toys data/settings_Scalar.cfg
```

The analysis steps are listed and explained below.
 - Workspace: create a workspace using the model specified in the config file.
 - AddDataToWS: add a dataset in MxAODs to the workspace. 
 - PlotWS: plot the signal and background fits to data from a workspace.
 - GlobalP0Toys: create a pseudo-experiment ensemble to get the global Z value.
 - GlobalP0Analysis: analyze the global toy ensemble to get the global z value.
 - LocalP0Analysis: analyze pseudo-experiments to cross-check asymptotic qmu,q0.
 - StatScan: plot the limits or p0 vs mass (toys or asymptotics)
 - ExtrapolateSig: extrapolate the 2015 signal into 2016 data with Asimov data
 - MassAnimation: animated GIF of the diphoton mass spectrum.

### Preparation for remote processing

Since large pseudo-experiment ensembles are CPU-intensive but highly 
parallelizable, the code can be easily configured to run on the GRID or other
remote resources. The step by step instructions below will guide you through
setting up a package for remote processing. 

#### Copy inputs and executables

Start in the package directory, and copy the necessary executables into a new 
directory:

```
mkdir forTar
cp *.pcm forTar/
cp bin/GlobalP0Toys_Faster forTar/
cp bin/GenericToys forTar/
cp data/settings_Graviton.cfg forTar/
cp $WORKSPACE forTar/
cp scripts/run.sh forTar/
```

You should substitute `data/settings_Graviton.cfg` for the config file of 
interest to you. Similarly, you should use the "WorkspaceFile" entry in the .cfg
file in place of "$WORKSPACE". Finally, if the executable changes, you should 
substitute that for `bin/GenericToys`. 

#### Modify the config file

The next few steps cannot be automated as easily. Open the settings file that
you just copied to the `forTar/` directory. Find the `WorkspaceFile` setting,
and copy that file to the `forTar/` directory. Then change the `WorkspaceFile`
to point just to the filename, as in the example below.

Then update your config file:
```
MasterOutput:		.
PackageLocation:	.
```
Also take a gander at the `JobName` in the config file before closing it, as it 
will be important for the next step. 

#### Modify the job script

Open the file `forTar/run.sh`. Make sure that the executable name 
`GlobalP0Toys_Faster` and config file name `settings_Graviton.cfg` are set as 
desired. Then check that the output location `Graviton_Mar9` corresponds to the
`JobName` setting from the config file.

#### Create and test the tar file

Create a `.tar` file with everything for the remote job. Then, test it using 
the `run.sh` script with a seed of 95487 and 1 toy. WARNING! Be absolutely sure
that the number of retries is set to the desired value in your config file 
(usually >=50 for remote jobs). NumRetries=1 is good enough for a test. 
```
cd forTar/
tar zcf Cocoon.tar *
source run.sh 95487 1
```

Assuming everything worked, you are ready to waste thousands of CPU hours!
