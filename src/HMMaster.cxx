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
//    - Workspace                                                             //
//    - GlobalP0Toys                                                          //
//    - GlobalP0Analysis                                                      //
//    - LocalP0Analysis                                                       //
//    - StatScan                                                              //
//    - ExtrapolateSig                                                        //
//    - MassAnimation + BSUB, GIF                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "HMMaster.h"

/**
   -----------------------------------------------------------------------------
   Submits the animated mass GIF jobs to the lxbatch server. 
   @param exeConfigFile - The config file.
   @param exeOption - The job options for the executable.
   @param exeMinFrame - The first frame to print.
   @param exeMaxFrame - The upper limit on the frame number for this job.
*/
void submitGIFViaBsub(TString exeConfigFile, TString exeOption, int exeMinFrame,
		      int exeMaxFrame) {
  // Make directories for job info:
  TString dir = Form("%s/%s_MassAnimation",
		     (m_config->getStr("ClusterFileLocation")).Data(),
		     (m_config->getStr("JobName")).Data());
  TString out = Form("%s/out", dir.Data());
  TString err = Form("%s/err", dir.Data());
  TString exe = Form("%s/exe", dir.Data());
  system(Form("mkdir -vp %s", out.Data()));
  system(Form("mkdir -vp %s", err.Data()));
  system(Form("mkdir -vp %s", exe.Data()));
    
  // create .tar file with everything:
  if (m_isFirstJob) {
    system("mkdir ForJob");
    system(Form("cp %s/*.pcm ForJob/", 
		(m_config->getStr("PackageLocation")).Data()));
    system(Form("cp %s/bin/%s ForJob/",
		(m_config->getStr("PackageLocation")).Data(),
		(m_config->getStr("exeMassAnimation")).Data()));
    //system(Form("cp %s ForJob/",
    //		(m_config->getStr("jobScriptMassAnimation")).Data()));
    system(Form("cp %s ForJob/",(m_config->getStr("WorkspaceFile")).Data()));
    system(Form("cp -f %s/%s %s/jobFileMassAnimation.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptMassAnimation")).Data(),
		exe.Data()));
    system("tar zcf Cocoon.tar ForJob/*");
    system(Form("mv Cocoon.tar %s", exe.Data()));
    system("rm -rf ForJob");
    m_isFirstJob = false;
  }
  /*
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exeMassAnimation")).Data()));
    system(Form("chmod +x %s",
		(m_config->getStr("jobScriptMassAnimation")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("WorkspaceFile")).Data()));
    system(Form("cp -f %s/%s %s/jobFileMassAnimation.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptMassAnimation")).Data(),
		exe.Data()));
    system(Form("chmod +x %s/jobFileMassAnimation.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
    m_isFirstJob = false;
  }
  */
  TString inputFile = Form("%s/Cocoon.tar", exe.Data());
  TString nameOutFile = Form("%s/%s_%d.out", out.Data(),
			     (m_config->getStr("JobName")).Data(), exeMinFrame);
  TString nameErrFile = Form("%s/%s_%d.err", err.Data(),
			     (m_config->getStr("JobName")).Data(), exeMinFrame);
  
  // Here you define the arguments for the job script:
  
  TString nameJScript = Form("%s/jobFileMassAnimation.sh %s %s %s %s %s %d %d", 
			     exe.Data(), (m_config->getStr("JobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeMassAnimation")).Data(),
			     exeOption.Data(), exeMinFrame, exeMaxFrame);
  /*
  TString nameJScript = Form("jobFileMassAnimation.sh %s %s %s %s %s %d %d", 
			     (m_config->getStr("JobName")).Data(),
			     exeConfigFile.Data(), inputFile.Data(),
			     (m_config->getStr("exeMassAnimation")).Data(),
			     exeOption.Data(), exeMinFrame, exeMaxFrame);
  */  
  // submit the job:
  system(Form("bsub -q wisc -o %s -e %s %s", nameOutFile.Data(), 
	      nameErrFile.Data(), nameJScript.Data()));
}

/**
   -----------------------------------------------------------------------------
   Submits the pseudo-experiment jobs to the lxbatch server. 
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
    
  // create .tar file with everything:
  if (m_isFirstJob) {
    system(Form("tar zcf Cocoon.tar bin/%s", 
		(m_config->getStr("exePseudoExp")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("jobScriptPseudoExp")).Data()));
    system(Form("chmod +x %s",(m_config->getStr("WorkspaceFile")).Data()));
    system(Form("cp -f %s/%s %s/jobFilePseudoExp.sh", 
		(m_config->getStr("PackageLocation")).Data(), 
		(m_config->getStr("jobScriptPseudoExp")).Data(), exe.Data()));
    system(Form("chmod +x %s/jobFilePseudoExp.sh", exe.Data()));
    system(Form("mv Cocoon.tar %s", exe.Data()));
    m_isFirstJob = false;
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
  // Step 0.0: Create a workspace for fitting:
  if (masterOption.Contains("Workspace")) {
    std::cout << "HMMaster: Step 0.0 - Making a workspace." << std::endl;
    Workspace *ws = new Workspace(configFileName, 
				  m_config->getStr("WorkspaceOptions"));
    if (!ws->fitsAllConverged()) {
      std::cout << "HMMaster: Problem with workspace fit!" << std::endl;
      exit(0);
    }
    std::cout << "HMMaster: Successfully built workspace!" << std::endl;
  }
  
  //--------------------------------------//
  // Step 1.0: Create pseudoexperiment ensemble for global p0 study:
  if (masterOption.Contains("GlobalP0Toys")) {
    std::cout << "HMMaster: Step 1.0 - Create pseudoexperiments." << std::endl;
    
    int toySeed = m_config->getInt("toySeed");
    int nToysTotal = m_config->getInt("nToysTotal");
    int nToysPerJob = m_config->getInt("nToysPerJob");
    int increment = nToysPerJob;
    int highestSeed = toySeed + nToysTotal;
    
    for (int i_s = toySeed; i_s < highestSeed; i_s += increment) {
      submitPEViaBsub(fullConfigPath, m_config->getStr("PseudoExpOptions"),
		      i_s, nToysPerJob);
    }
    std::cout << "HMMaster: Submitted " << (int)(nToysTotal/nToysPerJob) 
	      << " total pseudo-experiments." << std::endl;
  }
  
  //--------------------------------------//
  // Step 1.1: Plot pseudo-experiment ensemble results for global significance:
  if (masterOption.Contains("GlobalP0Analysis")) {
    std::cout << "HMMaster: Step 1.1 - Plot global significance toy results."
	      << std::endl;
    system(Form("./bin/GlobalP0Analysis %s %s", fullConfigPath.Data(),
		m_config->getStr("GlobalP0AnalysisOptions").Data()));
  }
  
  //--------------------------------------//
  // Step 2.1: Plot pseudo-experiment ensemble results for local significance:
  if (masterOption.Contains("LocalP0Analysis")) {
    std::cout << "HMMaster: Step 2.1 - Plot local significance toy results."
	      << std::endl;
    system(Form("./bin/LocalP0Analysis %s %s", fullConfigPath.Data(),
		m_config->getStr("LocalP0AnalysisOptions").Data()));
  }
  
  //--------------------------------------//
  // Step 3.1: Plot limits or p0 from asymptotics or pseudo-experiments:
  if (masterOption.Contains("StatScan")) {
    std::cout << "HMMaster: Step 3.0 - Plot asymptotic or toy limits or p0."
	      << std::endl;
    system(Form("./bin/PlotStatScan %s %s", fullConfigPath.Data(),
		m_config->getStr("StatScanOptions").Data()));
  }

  //--------------------------------------//
  // Step 4.0: Extrapolate the 2015 excess using S+B or B-only hypothesis:
  if (masterOption.Contains("ExtrapolateSig")) {
    std::cout << "HMMaster: Step 4.0 - Plot toy limits." << std::endl;
    system(Form("./bin/ExtrapolateSig %s %s 0", fullConfigPath.Data(),
		m_config->getStr("ExtrapSigOptions").Data()));
    system(Form("./bin/ExtrapolateSig %s %s 1", fullConfigPath.Data(),
    		m_config->getStr("ExtrapSigOptions").Data()));
  }
  
  //--------------------------------------//
  // Step 5.0: Create the mass animation:
  if (masterOption.Contains("MassAnimation")) {
    TString animationOptions = m_config->getStr("AnimationOptions");
    if (masterOption.Contains("MassAnimationBSUB")) {
      int framesPerJob = m_config->getInt("AnimationFramesPerJob");
      for (int i_f = 0; i_f < m_config->getInt("AnimationFrames");
	   i_f += framesPerJob) {
	submitGIFViaBsub(fullConfigPath, animationOptions, i_f, 
			 (i_f+framesPerJob));
      }
    }
    else if (masterOption.Contains("MassAnimationGIF")) {
      MassAnimation *animation = new MassAnimation(fullConfigPath, 
						   animationOptions);
      animation->setNFrames(m_config->getInt("AnimationFrames"));
      animation->makeGIF();
    }
    else {
      std::cout << "HMMaster: Step 5.0 - Plot mass GIF." << std::endl;
      MassAnimation *animation = new MassAnimation(fullConfigPath, 
						   animationOptions);
      animation->setNFrames(m_config->getInt("AnimationFrames"));
      animation->getDataForFrames();
      animation->makeAllFrames();
      animation->makeGIF();
    }
  }
  
  return 0;
}
