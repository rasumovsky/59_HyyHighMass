////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: ExtrapolateSig.cxx                                                  //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 31/03/2016                                                          //
//                                                                            //
//  This program extrapolates a given excess to higher luminosities in order  //
//  to determine when it might be discoverable (or excludable).               //
//                                                                            //
//  Options:                                                                  //
//  - "FromFile": Load values from text file.                                 //
//  - "Only2016": Look at standalone significance of 2016 dataset.            //
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
  if (argc < 4) {
    std::cout << "Usage: " << argv[0] << " <configFile> <options> <muHypothesis"
	      << std::endl;
    exit(0);
  }
  
  // Clock the toys:
  clock_t time;
  time = clock();
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString options = argv[2];
  int muHypothesis = atoi(argv[3]);

  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  
  // Construct the output directory:
  TString outputDir = Form("%s/%s/ExtrapolateSig", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Get the luminosity points:
  std::vector<double> lumiValues = config->getNumV("SigExtrapLuminosities");
  
  // Store Z0 values to be used in TGraph:
  double xValues[100] = {0.0};
  double xValuesSqrt[100] = {0.0};
  double yValues_Z0[100] = {0.0};
  double yErrorHi_Z0[100] = {0.0};
  double yErrorLo_Z0[100] = {0.0};
  
  double yValues_CL[100] = {0.0};
  double yErrorHi_CL[100] = {0.0};
  double yErrorLo_CL[100] = {0.0};
  
  TString fileTag = options.Contains("Only2016") ? "Only2016" : "Combined";
  TString textFileName
    = Form("%s/extrapolationValues_Mu%d_%s_%s.txt",outputDir.Data(),
	   muHypothesis,(config->getStr("AnalysisType").Data()),fileTag.Data());
  
  //----------------------------------------//
  // Load luminosity scaling from file:
  if (options.Contains("FromFile")) {
    std::ifstream inputLumiText(textFileName);
    if (inputLumiText.is_open()) {
      TString text = ""; int val0 = 0; double val1 = 0.0; double val2 = 0.0;
      while (inputLumiText >> text >> val0 >> val1 >> val2) {
	xValues[val0] = val1;
	xValuesSqrt[val0] = sqrt(xValues[val0]);
	if (text.Contains("z0Nominal")) yValues_Z0[val0] = val2;
	else if (text.Contains("z0Hi")) yErrorHi_Z0[val0] = val2;
	else if (text.Contains("z0Lo")) yErrorLo_Z0[val0] = val2;
	else if (text.Contains("CLNominal")) yValues_CL[val0] = val2;
	else if (text.Contains("CLHi")) yErrorHi_CL[val0] = val2;
        else if (text.Contains("CLLo")) yErrorLo_CL[val0] = val2;
	else {
	  std::cout << "ERROR! No matching type" << std::endl;
	  exit(0);
	}
      }
    }
    inputLumiText.close();
  }
  
  //----------------------------------------//
  // Create luminosity scaling from scratch:
  else {
    // Create an output text file:
    std::ofstream outputLumiText(textFileName);
    
    // Loop over luminosities to get expectation.
    std::cout << "ExtrapolateSig: Looping over lumi for extraploation" << endl;
    for (int i_l = 0; i_l < (int)lumiValues.size(); i_l++) {
      std::cout << "ExtrapolateSig: Lumi point" << i_l << " of " 
		<< lumiValues.size() << std::endl;
      
      double lumiSF = (lumiValues[i_l] / config->getNum("AnalysisLuminosity"));
      double lumiRelSF = ((lumiValues[i_l]-config->getNum("AnalysisLuminosity"))
			  / config->getNum("AnalysisLuminosity"));
      // Loop over errors:
      for (int i_e = -1; i_e < 2; i_e++) {
	
	// Use maximum-likelihood snapshot: 
	TString snapshotName = config->getStr("WorkspaceSnapshotMu1");
	
	// Load model, data, etc. from workspace:
	TFile inputFile(config->getStr("WorkspaceFile"), "read");
	RooWorkspace *workspace 
	  = (RooWorkspace*)inputFile.Get(config->getStr("WorkspaceName"));
	
	// The statistics class, for calculating qMu etc. 
	TestStat *testStat = new TestStat(configFile, "new", workspace);
	
	// Set the PoI ranges for this study:
	std::vector<TString> listPoI = config->getStrV("WorkspacePoIs");
	for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
	  std::vector<double> currRange
	    = config->getNumV(Form("ScanPoIRange_%s", (listPoI[i_p]).Data()));
	  if (testStat->theWorkspace()->var(listPoI[i_p])) {
	    testStat->theWorkspace()->var(listPoI[i_p])
	      ->setRange(currRange[0], currRange[1]);
	  }
	  else {
	    std::cout << "StatScan: Workspace has no variable " 
		      << listPoI[i_p].Data() << std::endl;
	    exit(0);
	  }
	}
	
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
	std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
	std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
	double snapshotXSection = 1.0;
	
	// For a maximum-likelihood fit, set the parameters of interest using 
	// mu=1 snapshot. The mu=1 fit params of interest are used even for the 
	// mu=0 parameters of interest because the background-only fits and toys
	// should still have the spurious signal randomized and fitted at that
	// mass.
	testStat->setNominalSnapshot(snapshotName);
	const RooArgSet *snapshotMu1 = workspace->getSnapshot(snapshotName);
	TIterator *snapMu1Iter = snapshotMu1->createIterator();
	RooRealVar *snapParMu1 = NULL;
	while ((snapParMu1 = (RooRealVar*)snapMu1Iter->Next())) {
	  TString currName = snapParMu1->GetName();
	  for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
	    if (currName.EqualTo(listPoI[i_p])) {
	      if (currName.EqualTo(poiForNorm)) {
		mapPoIMu0[listPoI[i_p]] = 0.0;
		
		// Do the extrapolation for 2016 significance only:
		if (options.Contains("Only2016")) {
		  // There is no data in the new dataset:
		  if (muHypothesis == 0) {
		    mapPoIMu1[listPoI[i_p]] = 0.0;
		    if (i_e == 0) {
		      snapshotXSection = lumiSF * snapParMu1->getVal();
		    }
		    else if (i_e == -1) {
		      snapshotXSection = lumiSF * (snapParMu1->getVal() - 
						   snapParMu1->getError());
		    }
		    else if (i_e == 1) {
		      snapshotXSection = lumiSF * (snapParMu1->getVal() +
						   snapParMu1->getError());
		    }
		  }
		  // The signal also appears in the new data:
		  else {
		    if (i_e == 0) {
		      mapPoIMu1[listPoI[i_p]] = lumiSF * snapParMu1->getVal();
		    }
		    else if (i_e == -1) {
		      mapPoIMu1[listPoI[i_p]]
			= lumiSF*(snapParMu1->getVal()-snapParMu1->getError());
		    }
		    else if (i_e == 1) {
		      mapPoIMu1[listPoI[i_p]]
			= lumiSF*(snapParMu1->getVal()+snapParMu1->getError());
		    }
		  }
		}
		// Do the extrapolation for 2015+2016 significance:
		else {
		  // There is no data in the new dataset:
		  if (muHypothesis == 0) {
		    if (i_e == 0) {
		      mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
		    }
		    else if (i_e == -1) {
		      mapPoIMu1[listPoI[i_p]] = (snapParMu1->getVal() -
						 snapParMu1->getError());
		    }
		    else if (i_e == 1) {
		      mapPoIMu1[listPoI[i_p]] = (snapParMu1->getVal() +
						 snapParMu1->getError());
		    }
		  }
		  
		  // The signal also appears in the new data:
		  else {
		    if (i_e == 0) {
		      mapPoIMu1[listPoI[i_p]] = lumiSF * snapParMu1->getVal();
		    }
		    else if (i_e == -1) {
		      mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal() + 
			(lumiRelSF*((snapParMu1->getVal() - 
				     snapParMu1->getError())));
		    }
		    else if (i_e == 1) {
		      mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal() + 
			(lumiRelSF*(snapParMu1->getVal() + 
				    snapParMu1->getError()));
		    }
		  }
		}
	      }
	      else {
		mapPoIMu0[listPoI[i_p]] = snapParMu1->getVal();
		mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
	      }
	    }
	  }
	}
	
	// If we are assuming no signal in 2016 data, then the amount of signal
	// should be constant (only what we saw in 2015):
	testStat->scaleAsimovData(lumiSF, config->getStrV("ParamsToScale"));
	testStat->createAsimovData(1.0, snapshotName, mapPoIMu1);
	
	//----------------------------------------//
	// Perform the fits:
	
	// Also force the mass and width to be constant always:
	testStat->setParam(poiForMass, mapPoIMu1[poiForMass], true);
	testStat->setParam(poiForWidth, mapPoIMu1[poiForWidth], true);
	
	// Also update the mu=1 cross-section to reflect the S+B expectation:
	if (muHypothesis == 0) mapPoIMu1[poiForNorm] = snapshotXSection;
	
	
	/*
	// Calculate the p0:
	std::vector<double> p0Values
	  = testStat->asymptoticP0(mapPoIMu1, "asimovDataMu1", snapshotName,
				   poiForNorm);
	double currZ0 = testStat->getZFromP(p0Values[0]);
	
	// Calculate the CL:
	//std::vector<double> CLValues 
	//= testStat->asymptoticCL(mapPoIMu1, "asimovDataMu1", snapshotName,
	//			   poiForNorm, config->getBool("UseQMuTilde"));
	double currCL = 0.0;//CLValues[0];
	*/
	
	
	// Mu = 0 fits:
	std::cout << "ExtrapolateSig: Mu=0 fit starting" << std::endl;
	double nllMu0 
	  = testStat->getFitNLL("asimovDataMu1", 0, true, mapPoIMu0, false);
	
	// Mu = 1 fits:
	std::cout << "ExtrapolateSig: Mu=1 fit starting" << std::endl;
	double nllMu
	  = testStat->getFitNLL("asimovDataMu1", 1, true, mapPoIMu1, false);
	
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
	
	// CL calculation:
	double currQMu = testStat->getQMuFromNLL(nllMu, nllMuHat, muHat, 1.0);
	double sigma = testStat->getSigma(currQMu, 1.0, muHat, false);
	double currCL = testStat->getCLFromQMu(currQMu, sigma, 0);
	
	std::cout << "DEBUGGER" << std::endl;
	std::cout << "\tnllMu0 " << nllMu0 << std::endl;
	std::cout << "\tnllMu " << nllMu << std::endl;
	std::cout << "\tnllMuHat " << nllMuHat << std::endl;
	std::cout << "\tprofiledNorm " << profiledNormalization << std::endl;
	std::cout << "\tnominalNorm " << nominalNormalization << std::endl;
	std::cout << "\tmuHat " << muHat << std::endl;
	std::cout << "\tcurrQ0 " << currQ0 << std::endl;
	std::cout << "\tcurrZ0 " << currZ0 << std::endl;
	std::cout << "\tcurrQMu " << currQMu << std::endl;
	std::cout << "\tcurrCL " << currCL << std::endl;
	std::cout << std::endl;
	
	
	// Fill arrays for graphs:
	if (options.Contains("Only2016")) xValues[i_l] = lumiValues[i_l]/1000.0;
	else {
	  // Total uminosity - 2015 luminosity = 2016 luminosity:
	  xValues[i_l] = ((lumiValues[i_l] - 
			   config->getNum("AnalysisLuminosity")) / 1000.0);
	}
	xValuesSqrt[i_l] = sqrt(xValues[i_l]);
	if (i_e == 0) { // Nominal signal
	  yValues_Z0[i_l] = currZ0;
	  yValues_CL[i_l] = currCL;
	  outputLumiText << "z0Nominal " << i_l << " " << xValues[i_l] << " " 
			 << yValues_Z0[i_l] << std::endl;
	  outputLumiText << "CLNominal " << i_l << " " << xValues[i_l] << " " 
			 << yValues_CL[i_l] << std::endl;
	}
	else if (i_e == 1) { // Signal + fit error
	  yErrorHi_Z0[i_l] = currZ0;
	  yErrorHi_CL[i_l] = currCL;
	  outputLumiText << "z0Hi " << i_l << " " << xValues[i_l] << " " 
			 << yErrorHi_Z0[i_l] << std::endl;
	  outputLumiText << "CLHi " << i_l << " " << xValues[i_l] << " " 
			 << yErrorHi_CL[i_l] << std::endl;
	}
	else if (i_e == -1) { // Signal - fit error
	  yErrorLo_Z0[i_l] = currZ0;
	  yErrorLo_CL[i_l] = currCL;
	  outputLumiText << "z0Lo " << i_l << " " << xValues[i_l] << " " 
			 << yErrorLo_Z0[i_l] << std::endl;
	  outputLumiText << "CLLo " << i_l << " " << xValues[i_l] << " " 
			 << yErrorLo_CL[i_l] << std::endl;
	}
	
	// Close the input file before the loop repeats:
	inputFile.Close();
	delete testStat;
	delete workspace; 
      }
    }
    outputLumiText.close();
  }
  
  //----------------------------------------//
  // Print TGraph of sensitivity:
  
  // First re-format the y-errors:
  for (int i_l = 0; i_l < (int)lumiValues.size(); i_l++) {
    // Format Z0:
    yErrorHi_Z0[i_l] = fabs(yErrorHi_Z0[i_l] - yValues_Z0[i_l]);
    yErrorLo_Z0[i_l] = fabs(yValues_Z0[i_l] - yErrorLo_Z0[i_l]);
    if (std::isnan(yErrorLo_Z0[i_l]) && std::isnan(yErrorHi_Z0[i_l])) {
      yErrorLo_Z0[i_l] = 0.0;
      yErrorHi_Z0[i_l] = 0.0;
    }
    else if (std::isnan(yErrorLo_Z0[i_l])) yErrorLo_Z0[i_l] = yErrorHi_Z0[i_l];
    else if (std::isnan(yErrorHi_Z0[i_l])) yErrorHi_Z0[i_l] = yErrorLo_Z0[i_l];
    
    
    // Format CL:
    yErrorHi_CL[i_l] = fabs(yErrorHi_CL[i_l] - yValues_CL[i_l]);
    yErrorLo_CL[i_l] = fabs(yValues_CL[i_l] - yErrorLo_CL[i_l]);
    if (std::isnan(yErrorLo_CL[i_l]) && std::isnan(yErrorHi_CL[i_l])) {
      yErrorLo_CL[i_l] = 0.0;
      yErrorHi_CL[i_l] = 0.0;
    }
    else if (std::isnan(yErrorLo_CL[i_l])) yErrorLo_CL[i_l] = yErrorHi_CL[i_l];
    else if (std::isnan(yErrorHi_CL[i_l])) yErrorHi_CL[i_l] = yErrorLo_CL[i_l];
  }
  
  // Create a graph with error bars:
  TGraphAsymmErrors *gZ0vsLumi_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValues, yValues_Z0, 0, 0,
			    yErrorLo_Z0, yErrorHi_Z0);
  TGraphAsymmErrors *gCLvsLumi_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValues, yValues_CL, 0, 0,
			    yErrorLo_CL, yErrorHi_CL);
  TGraphAsymmErrors *gZ0vsLumiSqrt_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValuesSqrt, yValues_Z0, 0,
			    0, yErrorLo_Z0, yErrorHi_Z0);
  TGraphAsymmErrors *gCLvsLumiSqrt_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValuesSqrt, yValues_CL, 0,
			    0, yErrorLo_CL, yErrorHi_CL);
  gZ0vsLumi_err->SetNameTitle(Form("gZOvsLumi_err_Mu%d", muHypothesis),
			      Form("gZOvsLumi_err_Mu%d", muHypothesis));
  gCLvsLumi_err->SetNameTitle(Form("gCLvsLumi_err_Mu%d", muHypothesis),
			      Form("gCLvsLumi_err_Mu%d", muHypothesis));
  gZ0vsLumiSqrt_err->SetNameTitle(Form("gZOvsLumiSqrt_err_Mu%d", muHypothesis),
				  Form("gZOvsLumiSqrt_err_Mu%d", muHypothesis));
  gCLvsLumiSqrt_err->SetNameTitle(Form("gCLvsLumiSqrt_err_Mu%d", muHypothesis),
				  Form("gCLvsLumiSqrt_err_Mu%d", muHypothesis));
  if (muHypothesis == 1) {
    gZ0vsLumi_err->SetFillColor(kBlue+1);
    gZ0vsLumi_err->SetFillStyle(3245);
    gZ0vsLumiSqrt_err->SetFillColor(kBlue+1);
    gZ0vsLumiSqrt_err->SetFillStyle(3245);
    gCLvsLumi_err->SetFillColor(kBlue+1);
    gCLvsLumi_err->SetFillStyle(3245);
    gCLvsLumiSqrt_err->SetFillColor(kBlue+1);
    gCLvsLumiSqrt_err->SetFillStyle(3245);
  }
  else {
    gZ0vsLumi_err->SetFillColor(kRed+1);
    gZ0vsLumi_err->SetFillStyle(3254);
    gZ0vsLumiSqrt_err->SetFillColor(kRed+1);
    gZ0vsLumiSqrt_err->SetFillStyle(3254);
    gCLvsLumi_err->SetFillColor(kRed+1);
    gCLvsLumi_err->SetFillStyle(3254);
    gCLvsLumiSqrt_err->SetFillColor(kRed+1);
    gCLvsLumiSqrt_err->SetFillStyle(3254);
  }
  TString xAxisName = "Luminosity at 13 TeV in 2016 [fb^{-1}]";
  TString xAxisSqrtName = "#sqrt{L} at 13 TeV in 2016 [#sqrt{fb^{-1}}]";
  gZ0vsLumi_err->GetXaxis()->SetTitle(xAxisName);
  gZ0vsLumiSqrt_err->GetXaxis()->SetTitle(xAxisSqrtName);
  gCLvsLumi_err->GetXaxis()->SetTitle(xAxisName);
  gCLvsLumiSqrt_err->GetXaxis()->SetTitle(xAxisSqrtName);
  if (options.Contains("Only2016")) {
    gZ0vsLumi_err->GetYaxis()->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gZ0vsLumiSqrt_err->GetYaxis()
      ->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gCLvsLumi_err->GetYaxis()->SetTitle("#it{CL} exclusion with 2016 data");
    gCLvsLumiSqrt_err->GetYaxis()->SetTitle("#it{CL} exclusion with 2016 data");
  }
  else {
    gZ0vsLumi_err->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gZ0vsLumiSqrt_err->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gCLvsLumi_err->GetYaxis()->SetTitle("Combined #it{CL} exclusion");
    gCLvsLumiSqrt_err->GetYaxis()->SetTitle("Combined #it{CL} exclusion");
  }
  
  // Create a graph with the median value:
  TGraph *gZ0vsLumi_nom 
    = new TGraph((int)lumiValues.size(), xValues, yValues_Z0);
  TGraph *gCLvsLumi_nom 
    = new TGraph((int)lumiValues.size(), xValues, yValues_CL);
  TGraph *gZ0vsLumiSqrt_nom
    = new TGraph((int)lumiValues.size(), xValuesSqrt, yValues_Z0);
  TGraph *gCLvsLumiSqrt_nom
    = new TGraph((int)lumiValues.size(), xValuesSqrt, yValues_CL);
  gZ0vsLumi_nom->SetNameTitle(Form("gZOvsLumi_nom_Mu%d", muHypothesis),
			      Form("gZOvsLumi_nom_Mu%d", muHypothesis));
  gZ0vsLumiSqrt_nom->SetNameTitle(Form("gZOvsLumiSqrt_nom_Mu%d", muHypothesis),
				  Form("gZOvsLumiSqrt_nom_Mu%d", muHypothesis));
  gCLvsLumi_nom->SetNameTitle(Form("gCLvsLumi_nom_Mu%d", muHypothesis),
			      Form("gCLvsLumi_nom_Mu%d", muHypothesis));
  gCLvsLumiSqrt_nom->SetNameTitle(Form("gCLvsLumiSqrt_nom_Mu%d", muHypothesis),
				  Form("gCLvsLumiSqrt_nom_Mu%d", muHypothesis));
  if (muHypothesis == 1) {
    gZ0vsLumi_nom->SetLineColor(kBlue+1);
    gZ0vsLumiSqrt_nom->SetLineColor(kBlue+1);
    gCLvsLumi_nom->SetLineColor(kBlue+1);
    gCLvsLumiSqrt_nom->SetLineColor(kBlue+1);
  }
  else {
    gZ0vsLumi_nom->SetLineColor(kRed+1);
    gZ0vsLumiSqrt_nom->SetLineColor(kRed+1);
    gCLvsLumi_nom->SetLineColor(kRed+1);
    gCLvsLumiSqrt_nom->SetLineColor(kRed+1);
  }
  gZ0vsLumi_nom->SetLineWidth(2);
  gZ0vsLumiSqrt_nom->SetLineWidth(2);
  gCLvsLumi_nom->SetLineWidth(2);
  gCLvsLumiSqrt_nom->SetLineWidth(2);
  gZ0vsLumi_nom->GetXaxis()->SetTitle(xAxisName);
  gZ0vsLumiSqrt_nom->GetXaxis()->SetTitle(xAxisSqrtName);
  gCLvsLumi_nom->GetXaxis()->SetTitle(xAxisName);
  gCLvsLumiSqrt_nom->GetXaxis()->SetTitle(xAxisSqrtName);
  if (options.Contains("Only2016")) {
    gZ0vsLumi_nom->GetYaxis()->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gZ0vsLumiSqrt_nom->GetYaxis()
      ->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gCLvsLumi_nom->GetYaxis()->SetTitle("#it{CL} exclusion with 2016 data");
    gCLvsLumiSqrt_nom->GetYaxis()->SetTitle("#it{CL} exclusion with 2016 data");
  }
  else {
    gZ0vsLumi_nom->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gZ0vsLumiSqrt_nom->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gCLvsLumi_nom->GetYaxis()->SetTitle("Combined #it{CL} exclusion");
    gCLvsLumiSqrt_nom->GetYaxis()->SetTitle("Combined #it{CL} exclusion");
  }
  
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  gZ0vsLumi_err->Draw("A3");
  gZ0vsLumi_nom->Draw("LSAME");
  
  // 5 sigma line:
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(kRed+1);
  line->DrawLine(gZ0vsLumi_err->GetXaxis()->GetXmin(), 5.0,
		 gZ0vsLumi_err->GetXaxis()->GetXmax(), 5.0);
  
  // Create a TLegend:
  TLegend leg(0.20, 0.75, 0.47, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if ((config->getStr("AnalysisType")).EqualTo("Graviton")) {
    leg.AddEntry(gZ0vsLumi_nom,
		 "Profiled #sigma#timesBR(G*#rightarrow#gamma#gamma)", "l");
  }
  else {
    leg.AddEntry(gZ0vsLumi_nom,
		 "Profiled #sigma#timesBR(X#rightarrow#gamma#gamma)", "l");
  }
  leg.AddEntry(gZ0vsLumi_err, "Profiled #pm error", "F");
  leg.Draw("SAME");
  
  // Print Z0 canvas:
  can->Print(Form("%s/plot_Z0vsLumi_Mu%d_%s_%s.eps", outputDir.Data(),
		  muHypothesis, (config->getStr("AnalysisType").Data()),
		  fileTag.Data()));

  // Print CL canvas:
  can->Clear();
  gCLvsLumi_err->Draw("A3");
  gCLvsLumi_nom->Draw("LSAME");
  line->DrawLine(gCLvsLumi_err->GetXaxis()->GetXmin(), 0.95,
		 gCLvsLumi_err->GetXaxis()->GetXmax(), 0.95);
  leg.Draw("SAME");
  can->Print(Form("%s/plot_CLvsLumi_Mu%d_%s_%s.eps", outputDir.Data(),
		  muHypothesis, (config->getStr("AnalysisType").Data()),
		  fileTag.Data()));
  
  // Also save a ROOT file:
  TFile *outputFile = new TFile(Form("%s/sigExtrap_Mu%d_%s_%s.root",
				     outputDir.Data(), muHypothesis,
				     (config->getStr("AnalysisType").Data()),
				     fileTag.Data()), "RECREATE");
  gZ0vsLumi_err->Write();
  gZ0vsLumi_nom->Write();
  gZ0vsLumiSqrt_err->Write();
  gZ0vsLumiSqrt_nom->Write();
  gCLvsLumi_err->Write();
  gCLvsLumi_nom->Write();
  gCLvsLumiSqrt_err->Write();
  gCLvsLumiSqrt_nom->Write();
  can->Write();
  outputFile->Close();
  
  // Clock the toys for final print summary:
  time = clock() - time;
  std::cout << "\nExtrapolateSig: procedure concluded." << std::endl;
  printf("\tProgram required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return 0;
}
