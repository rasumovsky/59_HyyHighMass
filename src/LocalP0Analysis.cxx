////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: LocalP0Analysis.cxx                                                //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 09/03/2016                                                          //
//                                                                            //
//  Calculates the global significance from a background-only toy MC ensemble.//
//                                                                            //
//  Macro options:                                                            //
//  - PlotGauss                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "ToyAnalysis.h"
#include "TPolyLine.h"

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options: "New","FromFile","toy","asymptotic","NEvents"
   @param resMass - The resonance mass.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <options>" 
	      << std::endl;
    exit(0);
  }
  
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString jobName = config->getStr("JobName");
  TString anaType = config->getStr("AnalysisType");
  TString outputDir = Form("%s/%s/LocalP0Analysis", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Load the toy analysis class, which does the analysis of toy MC jobs:
  ToyAnalysis *toyAna = new ToyAnalysis(configFile, "None");
  toyAna->setOutputDir(outputDir);
  
  std::vector<TString> fitTypes;
  fitTypes.clear();
  fitTypes.push_back("0");
  fitTypes.push_back("1");
  fitTypes.push_back("Free");
  toyAna->setFitTypes(fitTypes);
  
  // ToyImportSamp
  // ToyForScan
  // ToyMLPoint
  
  if (options.Contains("ToyForScan")) {
    toyAna->loadToy(0, Form("%s/%s/GenericToys/single_files/toy_mu0*ForScan*", 
			    (config->getStr("MasterOutput")).Data(),
			    (config->getStr("JobName")).Data()));
    toyAna->loadToy(1, Form("%s/%s/GenericToys/single_files/toy_mu1*ForScan*", 
			    (config->getStr("MasterOutput")).Data(),
			    (config->getStr("JobName")).Data()));
  }
  else {
    toyAna->loadToy(0, Form("%s/%s/GenericToys/single_files/toy_mu0*", 
			    (config->getStr("MasterOutput")).Data(),
			    (config->getStr("JobName")).Data()));
    toyAna->loadToy(1, Form("%s/%s/GenericToys/single_files/toy_mu1*", 
			    (config->getStr("MasterOutput")).Data(),
			    (config->getStr("JobName")).Data()));
  }
  
  if (!(toyAna->areInputFilesOK())) {
    std::cout << "LocalP0Analysis: ERROR loading toys." << std::endl;
    exit(0);
  }
  
  // Plot the toy MC nuisance parameter, global observables, and PoI:
  std::vector<TString> namesGlobs = toyAna->getNamesGlobalObservables();
  std::vector<TString> namesNuis = toyAna->getNamesNuisanceParameters();
  std::vector<TString> namesPars = toyAna->getNamesPoI();
  for (int i_g = 0; i_g < (int)namesGlobs.size(); i_g++) {
    if (!(namesGlobs[i_g]).Contains("gamma_stat_channel_bin")) {
      toyAna->plotHist(namesGlobs[i_g], 0);// Mu=0 toy data
      toyAna->plotHist(namesGlobs[i_g], 1);// Mu=1 data
    }
  }
  for (int i_n = 0; i_n < (int)namesNuis.size(); i_n++) {
    if (!(namesNuis[i_n]).Contains("gamma_stat_channel_bin")) {
      toyAna->plotHist(namesNuis[i_n], 0);
      toyAna->plotHist(namesNuis[i_n], 1);
    }
  }
  for (int i_p = 0; i_p < (int)namesPars.size(); i_p++) {
    toyAna->plotHist(namesPars[i_p], 0);
    toyAna->plotHist(namesPars[i_p], 1);
  }
  
  // Then plot statistics:
  toyAna->plotProfiledMu();
  toyAna->plotTestStat("QMu");
  toyAna->plotTestStat("Q0");
  toyAna->plotTestStatComparison("QMu");
  toyAna->plotTestStatComparison("Q0");
  
  delete config;
  return 0;
}
