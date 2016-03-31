////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: GenericToys.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 31/03/2016                                                          //
//                                                                            //
//  This program extrapolates a given excess to higher luminosities in order  //
//  to determine when it might be discoverable (or excludable).               //
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
   The main method tosses toys and saves data in a TTree.
   @param configFile - The name of the analysis config file.
   @param options - The options (see header note).
*/
int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <configFile> <options>" << std::endl;
    exit(0);
  }
  
  // Clock the toys:
  clock_t time;
  time = clock();
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString options = argv[2];
  
  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  
  // Construct the output directory:
  TString outputDir = Form("%s/%s/GenericToys", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/", outputDir.Data()));
  
  // Get the luminosity points:
  std::vector<double> lumiValues = config->getNumV("SigExtrapLuminosities");
  
  // Create a graph:
  TGraph *gZ0vsLumi = new TGraph();
  gZ0vsLumi->SetNameTitle("gZOvsLumi","gZOvsLumi");
  
  //----------------------------------------//
  // Loop over luminosities to get expectation.
  std::cout << "ExtrapolateSig: Looping over lumi for extraploation" << endl;
  for (int i_l = 0; i_l < (int)lumiValues.size(); i_l++) {
    std::cout << "ExtrapolateSig: Lumi point" << i_l << " of " 
	      << lumiValues.size() << std::endl;
    
    // Use maximum-likelihood snapshot: 
    TString snapshotName = config->getStr("WorkspaceSnapshotMu1");
    
    
    // Load model, data, etc. from workspace:
    TFile inputFile(config->getStr("WorkspaceFile"), "read");
    RooWorkspace *workspace 
      = (RooWorkspace*)inputFile.Get(config->getStr("WorkspaceName"));
    
    // The statistics class, for calculating qMu etc. 
    TestStat *testStat = new TestStat(configFile, "new", workspace);
    
    // Turn off MC stat errors if requested for Graviton jobs:
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
    
    // Get the names of the parameters of interest:
    TString poiForNorm = config->getStr("PoIForNormalization");
    TString poiForMass = config->getStr("PoIForMass");
    TString poiForWidth = config->getStr("PoIForWidth");
    
    
    // For setting the PoI values in toy creation and fitting:
    std::vector<TString> listPoI = config->getStrV("WorkspacePoIs");
    std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
    std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
    
    // For a maximum-likelihood fit, set the parameters of interest using mu = 1
    // snapshot. The mu=1 fit parameters of interest are used even for the mu=0 
    // parameters of interest because the background-only fits and toys should 
    // still have the spurious signal randomized and fitted at the proper mass.
    testStat->setNominalSnapshot(snapshotName);
    const RooArgSet *snapshotMu1 = workspace->getSnapshot(snapshotName);
    TIterator *snapMu1Iter = snapshotMu1->createIterator();
    RooRealVar *snapParMu1 = NULL;
    while ((snapParMu1 = (RooRealVar*)snapMu1Iter->Next())) {
      TString currName = snapParMu1->GetName();
      for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
	if (currName.EqualTo(listPoI[i_p])) {
	  mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
	  if (currName.EqualTo(poiForNorm)) mapPoIMu0[listPoI[i_p]] = 0.0;
	  else mapPoIMu0[listPoI[i_p]] = snapParMu1->getVal();
	}
      }
    }
        
    // Create Asimov data:
    testStat->createAsimovData(1.0, snapshotName, mapPoIMu1, lumiValues[i_l]);
    
    //----------------------------------------//
    // Perform the fits:
    
    // Mu = 0 fits:
    std::cout << "ExtrapolateSig: Mu=0 fit starting" << std::endl;
    double nllMu0 
      = testStat->getFitNLL("asimovDataMu1", 0, true, mapPoIMu0, false);
    // Mu free fits:
    std::cout << "ExtrapolateSig: Mu-free fit starting" << std::endl;
    double nllMuHat
      = testStat->getFitNLL("asimovDataMu1", 1, false, mapPoIMu1, false);
    double profiledNormalization
      = (testStat->getPoIs())[(std::string)(poiForNorm)];
    
    // Significance calculation:
    double nominalNormalization = mapPoIMu1[poiForNorm];
    double muHat = (profiledNormalization / nominalNormalization);
    double currQ0 = testStat->getQ0FromNLL(nllMu0, nllMuHat, muHat);
    double currZ0 = testStat->getZ0FromQ0(currQ0);
    
    gZ0vsLumi->SetPoint(i_l, lumiValues[i_l], currZ0);

    // Close the input file before the loop repeats:
    inputFile.Close();
    delete testStat;
    delete workspace; 
  }
  
  //----------------------------------------//
  // Print TGraph of sensitivity:
  gZ0vsLumi->SetLineWidth(2);
  gZ0vsLumi->SetLineColor(kBlue+1);
  gZ0vsLumi->GetXaxis()->SetTitle("13 TeV Luminosity [fb-1]");
  gZ0vsLumi->GetYaxis()->SetTitle("z_{0}^{Local} [#sigma]");
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  gZ0vsLumi->Draw("ALP");
  can->Print(Form("%s/plot_Z0vsLumi.eps",outputDir.Data()));

  // Clock the toys for final print summary:
  time = clock() - time;
  std::cout << "\nExtrapolateSig: procedure concluded." << std::endl;
  printf("\tProgram required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return 0;
}
