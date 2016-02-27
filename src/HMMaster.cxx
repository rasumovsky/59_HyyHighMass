////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: HMMaster.cxx                                                        //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 27/02/2016                                                          //
//                                                                            //
//  This program is useful as an interface to the high-mass diphoton analysis //
//  tools. It centralizes the commands for creating inputs, plots, workspaces,//
//  and statistical results. Some of the commands will rely on accessing      //
//  classes (mass points, signal parameterization), while others will use     //
//  system commands to submit jobs to various clusters.                       //
//                                                                            //
//  MasterOption:                                                             //
//    - TossGlobalP0Toys                                                      //
//    - PlotGlobalP0Toys                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HMMaster.h"

/**
   -----------------------------------------------------------------------------
   Submits the mu limit jobs to the lxbatch server. 
   @param exeConfigFile - the config file.
   @param exeOption - the job options for the executable.
   @param exeSeed - the seed for the randomized dataset generation.
   @param exeToysPerJob - the number of toy datasets to create per job.
*/
void submitPEViaBsub(TString exeConfigFile, TString exeOption, int exeSeed, 
		     int exeToysPerJob) {
  //@param exeSignal - the signal to process in the executable.

  // Make directories for job info:
  TString dir = Form("%s/%s_PseudoExp",
		     (m_config->getStr("ClusterFileLocation")).Data(),
		     (m_config->getStr("JobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
  
  TString exeAna = m_config->getStr("AnalysisType");
  
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exePseudoExp")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("jobScriptPseudoExp")).Data()));
    system(Form("chmod +x %s/%s/DHWorkspace/rootfiles/workspaceDH_%s.root", 
		(m_config->getStr("MasterOutput")).Data(), 
		(m_config->getStr("JobName")).Data(), exeAna.Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptPseudoExp")).Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
  }
  
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/out/%s_%d.out", dir.Data(),
			     (m_config->getStr("JobName")).Data(), exeSeed);
  TString nameErrFile = Form("%s/err/%s_%d.err", dir.Data(),
			     (m_config->getStr("JobName")).Data(), exeSeed);
  
  // Here you define the arguments for the job script:
  TString nameJScript = Form("%s/jobFilePseudoExp.sh %s %s %s %s %s %d %d", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exePseudoExp")).Data(),
			     exeOption.Data(), exeSeed, exeToysPerJob);
  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   This is the main HMMaster method. Runs every analysis component.
   @param masterOption - The task to run (see document header or README).
   @param configFileName - The name of the config file. 
*/
int main (int argc, char **argv) {
  // Check arguments:
  if (argc < 3) {
    printf("\nUsage: %s <option> <configFileName>\n\n", argv[0]);
    exit(0);
  }
  
  // The job name and options (which analysis steps to perform):
  TString masterOption = argv[1];
  TString configFileName = argv[2];
  
  // For grid/bsub submission:
  m_isFirstJob = true;
  
  // Load the config class and file:
  std::cout << "HMMaster: Loading the global config file." << std::endl;
  m_config = new Config(configFileName);
  m_config->printDB();
  TString fullConfigPath = Form("%s/%s",
				(m_config->getStr("PackageLocation")).Data(),
				configFileName.Data());
  
  //--------------------------------------//
  // Step 5.0: Create pseudoexperiment ensemble:
  if (masterOption.Contains("TossGlobalP0Toys")) {
    std::cout << "HMMaster: Step 5.1 - Create pseudoexperiments." << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(fullConfigPath, m_config->getStr("PseudoExpOptions"),
		      i_s, nToysPerJob);
      m_isFirstJob = false;
    }
    std::cout << "HMMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 5.2: Plot pseudo-experiment ensemble results:
  if (masterOption.Contains("PlotGlobalP0Toys")) {
    std::cout << "HMMaster: Step 5.2 - Plot pseudoexperiment results."
	      << std::endl;
    //ToyAnalysis *ta = new ToyAnalysis(configFileName, "NONE");
  }
  
  return 0;
}
