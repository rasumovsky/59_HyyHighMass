////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: LocalP0Toys.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 09/03/2016                                                          //
//                                                                            //
//  This program tosses signal-plus-background or background-only toy MC.     //
//                                                                            //
//  Options:                                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "TestStat.h"
#include <time.h>
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

/**
   -----------------------------------------------------------------------------
   Convert key and value in a map to two separate vectors that are passed by 
   reference. Stupid, yeah...
   @param map - The input map with keys and values to be split.
   @names - A vector of names.
   @values - A vector of values.
*/
void mapToVectors(std::map<std::string,double> map, 
		  std::vector<std::string>& names, std::vector<double>& values){
  names.clear(); 
  values.clear();
  for (std::map<std::string,double>::iterator mapIter = map.begin(); 
       mapIter != map.end(); mapIter++) {
    names.push_back(mapIter->first);
    values.push_back(mapIter->second);
  }
}

/**
   -----------------------------------------------------------------------------
   The main method tosses toys and saves data in a TTree.
   @param configFile - The name of the analysis config file.
   @param options - The options (see header note).
   @param seed - The random seed for pseudoexperiment creation.
   @param toysPerJob - The number of pseudoexperiments to create per job.
   @param poiVal - The value of the parameter of interest to use.
*/
int main(int argc, char **argv) {
  if (argc < 6) {
    std::cout << "Usage: " << argv[0] 
	      << " <configFile> <options> <seed> <toysPerJob> <poiVal> " 
	      << std::endl;
    exit(0);
  }
  
  // Clock the toys:
  clock_t time;
  time = clock();
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString options = argv[2];
  int seed = atoi(argv[3]);
  int nToysPerJob = atoi(argv[4]);
  int inputPoIVal = atoi(argv[5]);
  
  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  TString anaType = config->getStr("AnalysisType");
  
  // Copy the input workspace file locally:
  TString originFile = config->getStr("WorkspaceFile");
  TString copiedFile = Form("workspace_%s.root", anaType.Data());
  system(Form("cp %s %s", originFile.Data(), copiedFile.Data()));
  
  // Construct the output directories:
  TString outputDir = Form("%s/%s/LocalP0Toys", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/err", outputDir.Data()));
  system(Form("mkdir -vp %s/log", outputDir.Data()));
  system(Form("mkdir -vp %s/single_files", outputDir.Data()));
  
  // Output TTree file:
  TString tempOutputFileName = Form("%s/single_files/toy_mu%i_%i.root",
				    outputDir.Data(), inputPoIVal, seed);
  TFile fOutputFile(tempOutputFileName, "recreate");
  TTree fOutputTree("toy", "toy");
  
  // Variables to store in the TTree:
  int bestFitUpdate;
  double numEvents;
  bool convergedMu1, convergedMu0, convergedMuFree;
  double profiledPOIVal;
  double nllMu0, nllMu1, nllMuFree, llrL1L0, llrL0Lfree, llrL1Lfree;
  std::vector<std::string> namesNP; namesNP.clear();
  std::vector<double> valuesNPMu0; valuesNPMu0.clear();
  std::vector<double> valuesNPMu1; valuesNPMu1.clear();
  std::vector<double> valuesNPMuFree; valuesNPMuFree.clear();
  std::vector<string> namesGlobs; namesGlobs.clear();
  std::vector<double> valuesGlobsMu0; valuesGlobsMu0.clear();
  std::vector<double> valuesGlobsMu1; valuesGlobsMu1.clear();
  std::vector<double> valuesGlobsMuFree; valuesGlobsMuFree.clear();
  std::vector<string> namesPoIs; namesPoIs.clear();
  std::vector<double> valuesPoIsMu0; valuesPoIsMu0.clear();
  std::vector<double> valuesPoIsMu1; valuesPoIsMu1.clear();
  std::vector<double> valuesPoIsMuFree; valuesPoIsMuFree.clear();
  std::vector<double> numEventsPerCate; numEventsPerCate.clear();
  std::vector<double> nllPerRetry; nllPerRetry.clear();
  
  fOutputTree.Branch("seed", &seed, "seed/I");
  fOutputTree.Branch("bestFitUpdate", &bestFitUpdate, "bestFitUpdate/I");
  fOutputTree.Branch("numEvents", &numEvents, "numEvents/D");
  fOutputTree.Branch("numEventsPerCate", &numEventsPerCate);
  fOutputTree.Branch("nllPerRetry", &nllPerRetry);
  fOutputTree.Branch("profiledPOIVal", &profiledPOIVal, "profiledPOIVal/D");
  fOutputTree.Branch("convergedMu0", &convergedMu0, "convergedMu0/O");
  fOutputTree.Branch("convergedMu1", &convergedMu1, "convergedMu1/O");
  fOutputTree.Branch("convergedMuFree", &convergedMuFree, "convergedMuFree/O");
  fOutputTree.Branch("nllMu0", &nllMu0, "nllMu0/D");
  fOutputTree.Branch("nllMu1", &nllMu1, "nllMu1/D");
  fOutputTree.Branch("nllMuFree", &nllMuFree, "nllMuFree/D");
  fOutputTree.Branch("llrL1L0", &llrL1L0, "llrL1L0/D");
  fOutputTree.Branch("llrL0Lfree", &llrL0Lfree, "llrL0Lfree/D");
  fOutputTree.Branch("llrL1Lfree", &llrL1Lfree, "llrL1Lfree/D");
  fOutputTree.Branch("namesNP", &namesNP);
  fOutputTree.Branch("valuesNPMu0", &valuesNPMu0);
  fOutputTree.Branch("valuesNPMu1", &valuesNPMu1);
  fOutputTree.Branch("valuesNPMuFree", &valuesNPMuFree);
  fOutputTree.Branch("namesGlobs", &namesGlobs);
  fOutputTree.Branch("valuesGlobsMu1", &valuesGlobsMu1);
  fOutputTree.Branch("valuesGlobsMu0", &valuesGlobsMu0);
  fOutputTree.Branch("valuesGlobsMuFree", &valuesGlobsMuFree);
  fOutputTree.Branch("namesPoIs", &namesPoIs);
  fOutputTree.Branch("valuesPoIsMu1", &valuesPoIsMu1);
  fOutputTree.Branch("valuesPoIsMu0", &valuesPoIsMu0);
  fOutputTree.Branch("valuesPoIsMuFree", &valuesPoIsMuFree);
  
  // Get the name of the normalization parameter:
  TString normalizationPoI = config->getStr("PoIForNormalization");
  
  //----------------------------------------//
  // Loop to generate pseudo experiments:
  std::cout << "LocalP0Toys: Generating " << nToysPerJob
	    << " toys with mu_DH = " << inputPoIVal << endl;
  for (int i_t = 0; i_t < nToysPerJob; i_t++) {
    std::cout << "LocalP0Toys: Starting toy " << i_t << " of " << nToysPerJob
	      << std::endl;
    
    // Load model, data, etc. from workspace:
    TFile inputFile(copiedFile, "read");
    RooWorkspace *workspace
      = (RooWorkspace*)inputFile.Get(config->getStr("WorkspaceName"));
    
    // The statistics class, for calculating qMu etc. 
    TestStat *testStat = new TestStat(configFile, "new", workspace);
        
    // Load the background-only snapshot (signal + background snapshot)
    // for background-only toys (signal + background toys).
    TString snapshotName = (inputPoIVal == 0) ? 
      config->getStr("WorkspaceSnapshotMu0") :
      config->getStr("WorkspaceSnapshotMu1");
    testStat->setNominalSnapshot(snapshotName);
    
    // For setting the PoI values in toy creation and fitting:
    std::vector<TString> listPoI = config->getStrV("WorkspacePoIs");
    std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
    std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
    
    // Set the parameters of interest using mu = 1 snapshot. The mu=1 fit 
    // parameters are used even for the mu=0 parameters of interest because the
    // background-only fits and toys should still have the spurious signal 
    // randomized and fitted at the proper mass.
    const RooArgSet *snapshotMu1
      = workspace->getSnapshot(config->getStr("WorkspaceSnapshotMu1"));
    TIterator *snapMu1Iter = snapshotMu1->createIterator();
    RooRealVar *snapParMu1 = NULL;
    while ((snapParMu1 = (RooRealVar*)snapMu1Iter->Next())) {
      TString currName = snapParMu1->GetName();
      for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
	if (currName.EqualTo(listPoI[i_p])) {
	  mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
	  std::cout << "From snapshot, setting " << listPoI[i_p] 
		    << " equal to " << mapPoIMu1[listPoI[i_p]]
		    << " for Mu=1 map " << std::endl;
	  if (currName.EqualTo(normalizationPoI)) {
	    mapPoIMu0[listPoI[i_p]] = 0.0;
	    std::cout << "From snapshot, setting " << listPoI[i_p] 
		      << " equal to " << mapPoIMu0[listPoI[i_p]]
		      << " for Mu=0 map " << std::endl;
	  }
	  else {
	    mapPoIMu0[listPoI[i_p]] = snapParMu1->getVal();
	    std::cout << "From snapshot, setting " << listPoI[i_p] 
		      << " equal to " << mapPoIMu0[listPoI[i_p]]
		      << " for Mu=0 map " << std::endl;
	  }
	}
      }
    }
    
    /*
      I don't think we want to set PoI ranges for the local p0 study

    // Set the PoI ranges for this study:
    for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
      std::vector<double> currRange
	= config->getNumV(Form("PoIRange_%s", (listPoI[i_p]).Data()));
      if (testStat->theWorkspace()->var(listPoI[i_p])) {
	testStat->theWorkspace()->var(listPoI[i_p])
	  ->setRange(currRange[0], currRange[1]);
      }
      else {
	std::cout << "LocalP0Toys: Workspace has no variable " << listPoI[i_p]
		  << std::endl;
	exit(0);
      }
    }
    */
        
    // Also turn off MC stat errors if requested for Graviton jobs:
    if (config->isDefined("TurnOffTemplateStat") && 
	config->getBool("TurnOffTemplateStat")) {
      testStat->theWorkspace()->loadSnapshot(snapshotName);
      const RooArgSet *nuisanceParameters
	= testStat->theModelConfig()->GetNuisanceParameters();
      TIterator *nuisIter = nuisanceParameters->createIterator();
      RooRealVar *nuisCurr = NULL;
      while ((nuisCurr = (RooRealVar*)nuisIter->Next())) {
	TString currName = nuisCurr->GetName();
	if (currName.Contains("gamma_stat_channel_bin")) {
	  testStat->setParam(currName, nuisCurr->getVal(), true);
	}
      }
    }
    
    // Create the pseudo data:
    RooDataSet *newToyData = NULL;
    if (inputPoIVal == 0) {
      //newToyData = 
      testStat->createPseudoData(seed, inputPoIVal, snapshotName, mapPoIMu0);
    }
    else {
      //newToyData = 
      testStat->createPseudoData(seed, inputPoIVal, snapshotName, mapPoIMu1);
    }
    
    numEvents = workspace->data("toyData")->sumEntries();
    numEventsPerCate = testStat->getNEventsToys();
    
    // Globs are only randomized once for each dataset:
    mapToVectors(testStat->getGlobalObservables(), namesGlobs, valuesGlobsMu0);
    mapToVectors(testStat->getGlobalObservables(), namesGlobs, valuesGlobsMu1);
    mapToVectors(testStat->getGlobalObservables(), namesGlobs, 
		 valuesGlobsMuFree);
    
    // Strategy settings (if defined):
    if (config->isDefined("FitOptions")) {
      testStat->setFitOptions(config->getStr("FitOptions"));
    }
    
    // Mu = 0 fits:
    std::cout << "LocalP0Toys: Mu=0 fit starting" << std::endl;
    nllMu0 = testStat->getFitNLL("toyData", 0, true, mapPoIMu0, false);
    convergedMu0 = testStat->fitsAllConverged();
    mapToVectors(testStat->getNuisanceParameters(), namesNP, valuesNPMu0);
    mapToVectors(testStat->getPoIs(), namesPoIs, valuesPoIsMu0);
    
    // Mu = 1 fits:
    std::cout << "LocalP0Toys: Mu=1 fit starting" << std::endl;
    nllMu1 = testStat->getFitNLL("toyData", 1, true, mapPoIMu1, false);
    convergedMu1 = testStat->fitsAllConverged();
    mapToVectors(testStat->getNuisanceParameters(), namesNP, valuesNPMu1);
    mapToVectors(testStat->getPoIs(), namesPoIs, valuesPoIsMu1);
    double nominalNormalization
      = (testStat->getPoIs())[(std::string)(normalizationPoI)];
    
    // Fix the mass and width before the mu-free fit:
    // Loop over the parameters of interest, and set all of those that are
    // not the cross-section constant and to the proper values
    for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
      // Set every parameter constant to mu=1 value EXCEPT cross-section:
      if (!(listPoI[i_p]).EqualTo(config->getStr("PoIForNormalization"))) {
	// get the value of the PoI to set from the previous fit:
	for (int i_n = 0; i_n < (int)namesPoIs.size(); i_n++) {
	  if (TString(namesPoIs[i_n]).EqualTo(listPoI[i_p])) {
	    testStat->setParam(listPoI[i_p], valuesPoIsMu1[i_n], true);
	    break;
	  }
	}
      }
    }
    
    // Mu free fits:
    std::cout << "LocalP0Toys: Mu-free fit starting" << std::endl;
    nllMuFree = testStat->getFitNLL("toyData", 1, false, mapPoIMu1, false);
    convergedMuFree = testStat->fitsAllConverged();
    mapToVectors(testStat->getNuisanceParameters(), namesNP, valuesNPMuFree);
    mapToVectors(testStat->getPoIs(), namesPoIs, valuesPoIsMuFree);
    double profiledNormalization
      = (testStat->getPoIs())[(std::string)(normalizationPoI)];
    
    // Get the profiled PoI value:
    profiledPOIVal = (profiledNormalization / nominalNormalization);
    
    // Calculate profile likelihood ratios:
    llrL1L0 = nllMu1 - nllMu0;
    llrL1Lfree = profiledPOIVal > 1.0 ? 0.0 : (nllMu1 - nllMuFree);
    llrL0Lfree = profiledPOIVal < 0.0 ? 0.0 : (nllMu0 - nllMuFree);
    
    // Fill the tree:
    fOutputTree.Fill();
    
    // Count the toys:
    seed++;
    fOutputTree.AutoSave("SaveSelf");
    
    // Close the input file before the loop repeats:
    inputFile.Close();
  }
  
  // Write the output file, delete local file copies:
  fOutputFile.cd();
  fOutputTree.Write();
  fOutputFile.Close();
  system(Form("rm %s", copiedFile.Data()));
  
  // Clock the toys:
  time = clock() - time;
  
  if (config->getBool("Verbose")) {
    std::cout << "\nLocalP0Toys_Faster: Toy procedure concluded." << std::endl;
    printf("\t%d toys required %d clock cycles (%f seconds).\n\n",
	   nToysPerJob, (int)time, ((float)time/CLOCKS_PER_SEC));
  }
  return 0;
}
