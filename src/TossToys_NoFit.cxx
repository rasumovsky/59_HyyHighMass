////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: TossToys_NoFit.cxx                                                  //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 02/03/2016                                                          //
//                                                                            //
//  This program just tosses and saves a bunch of toys without trying to fit. //
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
  TString outputDir = Form("%s/%s/TossToys_NoFit", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/err", outputDir.Data()));
  system(Form("mkdir -vp %s/log", outputDir.Data()));
  system(Form("mkdir -vp %s/single_files", outputDir.Data()));
  
  // Set the mass and width according to the given hypothesis:
  std::vector<TString> listPoI = config->getStrV("WorkspacePoIs");
  std::vector<double> inValPoIMu0 = config->getNumV("PoIValuesMu0");
  std::vector<double> inValPoIMu1 = config->getNumV("PoIValuesMu1");
  std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
  std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
  for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
    mapPoIMu0[listPoI[i_p]] = inValPoIMu0[i_p];
    mapPoIMu1[listPoI[i_p]] = inValPoIMu1[i_p];
  }
  
  // Load model, data, etc. from workspace:
  TFile inputFile(copiedFile, "read");
  RooWorkspace *workspace
    = (RooWorkspace*)inputFile.Get(config->getStr("WorkspaceName"));
  
  // The statistics class, for calculating qMu etc. 
  TestStat *testStat = new TestStat(configFile, "new", workspace);
  
  //----------------------------------------//
  // Loop to generate pseudo experiments:
  std::cout << "TossToys_NoFit: Generating " << nToysPerJob
	    << " toys with mu_DH = " << inputPoIVal << endl;
  for (int i_t = 0; i_t < nToysPerJob; i_t++) {
    std::cout << "TossToys_NoFit: Starting toy " << i_t << " of " << nToysPerJob
	      << std::endl;
    
    // Then create snapshot for Mu=1 or Mu=0 hypothesis! This must be re-done
    // for every toy job, since the signal hypothesis can change!
    // Set the PoI ranges for this study:
    for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
      std::vector<double> currRange
	= config->getNumV(Form("PoIRange_%s", (listPoI[i_p]).Data()));
      if (testStat->theWorkspace()->var(listPoI[i_p])) {
	testStat->theWorkspace()->var(listPoI[i_p])
	  ->setRange(currRange[0], currRange[1]);
      }
      else {
	std::cout << "TossToys_NoFit: Workspace has no variable "
		  << listPoI[i_p] << std::endl;
	exit(0);
      }
    }
    
    // Create the pseudo data:
    TString snapshotName = config->getStr("WorkspaceSnapshot");
    RooDataSet *newToyData = testStat->createPseudoData(seed, inputPoIVal,
							snapshotName, mapPoIMu0,
							i_t);
    // Count the toys:
    seed++;
  }
  workspace->importClassCode();
  workspace->writeToFile(Form("%s/workspaceWithToys_%s.root",
			      outputDir.Data(), anaType.Data()));
  
  inputFile.Close();
  
  // Write the output file, delete local file copies:
  system(Form("rm %s", copiedFile.Data()));
  
  // Clock the toys:
  time = clock() - time;

  if (config->getBool("Verbose")) {
    std::cout << "\nTossToys_NoFit: Toy procedure concluded." << std::endl;
    printf("\t%d toys required %d clock cycles (%f seconds).\n\n",
	   nToysPerJob, (int)time, ((float)time/CLOCKS_PER_SEC));
  }
  return 0;
}
