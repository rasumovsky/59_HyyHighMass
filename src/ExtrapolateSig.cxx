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
  
  // Store values to be used in TGraph:
  double xValues[100] = {0.0};
  double xValuesSqrt[100] = {0.0};
  double yValues[100] = {0.0};
  double yErrorHi[100] = {0.0};
  double yErrorLo[100] = {0.0};
  
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
	if (text.Contains("z0Nominal")) yValues[val0] = val2;
	else if (text.Contains("z0Hi")) yErrorHi[val0] = val2;
	else if (text.Contains("z0Lo")) yErrorLo[val0] = val2;
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
		  if (muHypothesis == 0) mapPoIMu1[listPoI[i_p]] = 0.0;
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
	
	// Fill arrays for graphs:
	if (options.Contains("Only2016")) xValues[i_l] = lumiValues[i_l]/1000.0;
	else {
	  // Total uminosity - 2015 luminosity = 2016 luminosity:
	  xValues[i_l] = ((lumiValues[i_l] - 
			   config->getNum("AnalysisLuminosity")) / 1000.0);
	}
	xValuesSqrt[i_l] = sqrt(xValues[i_l]);
	if (i_e == 0) { // Nominal signal
	  yValues[i_l] = currZ0;
	  outputLumiText << "z0Nominal " << i_l << " " << xValues[i_l] << " " 
			 << yValues[i_l] << std::endl;
	}
	else if (i_e == 1) { // Signal + fit error
	  yErrorHi[i_l] = currZ0;
	  outputLumiText << "z0Hi " << i_l << " " << xValues[i_l] << " " 
			 << yErrorHi[i_l] << std::endl;
	}
	else if (i_e == -1) { // Signal - fit error
	  yErrorLo[i_l] = currZ0;
	  outputLumiText << "z0Lo " << i_l << " " << xValues[i_l] << " " 
			 << yErrorLo[i_l] << std::endl;
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
    yErrorHi[i_l] = fabs(yErrorHi[i_l] - yValues[i_l]);
    yErrorLo[i_l] = fabs(yValues[i_l] - yErrorLo[i_l]);
    if (std::isnan(yErrorLo[i_l]) && std::isnan(yErrorHi[i_l])) {
      yErrorLo[i_l] = 0.0;
      yErrorHi[i_l] = 0.0;
    }
    else if (std::isnan(yErrorLo[i_l])) yErrorLo[i_l] = yErrorHi[i_l];
    else if (std::isnan(yErrorHi[i_l])) yErrorHi[i_l] = yErrorLo[i_l];
  }
  
  // Create a graph with error bars:
  TGraphAsymmErrors *gZ0vsLumi_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValues, yValues, 0, 0,
			    yErrorLo, yErrorHi);
  TGraphAsymmErrors *gZ0vsLumiSqrt_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValuesSqrt, yValues, 0, 0,
			    yErrorLo, yErrorHi);
  gZ0vsLumi_err->SetNameTitle(Form("gZOvsLumi_err_Mu%d", muHypothesis),
			      Form("gZOvsLumi_err_Mu%d", muHypothesis));
  gZ0vsLumiSqrt_err->SetNameTitle(Form("gZOvsLumiSqrt_err_Mu%d", muHypothesis),
				  Form("gZOvsLumiSqrt_err_Mu%d", muHypothesis));
  if (muHypothesis == 1) {
    gZ0vsLumi_err->SetFillColor(kBlue+1);
    gZ0vsLumi_err->SetFillStyle(3245);
    gZ0vsLumiSqrt_err->SetFillColor(kBlue+1);
    gZ0vsLumiSqrt_err->SetFillStyle(3245);
  }
  else {
    gZ0vsLumi_err->SetFillColor(kRed+1);
    gZ0vsLumi_err->SetFillStyle(3254);
    gZ0vsLumiSqrt_err->SetFillColor(kRed+1);
    gZ0vsLumiSqrt_err->SetFillStyle(3254);
  }
  gZ0vsLumi_err->GetXaxis()->SetTitle("Luminosity at 13 TeV in 2016 [fb^{-1}]");
  gZ0vsLumiSqrt_err->GetXaxis()
    ->SetTitle("#sqrt{L} at 13 TeV in 2016 [#sqrt{fb^{-1}}]");
  if (options.Contains("Only2016")) {
    gZ0vsLumi_err->GetYaxis()->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gZ0vsLumiSqrt_err->GetYaxis()
      ->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
  }
  else {
    gZ0vsLumi_err->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gZ0vsLumiSqrt_err->GetYaxis()
      ->SetTitle("Combined Z_{0}^{Local} [#sigma]");
  }
  
  // Create a graph with the median value:
  TGraph *gZ0vsLumi_nom = new TGraph((int)lumiValues.size(), xValues, yValues);
  TGraph *gZ0vsLumiSqrt_nom
    = new TGraph((int)lumiValues.size(), xValues, yValues);
  gZ0vsLumi_nom->SetNameTitle(Form("gZOvsLumi_nom_Mu%d", muHypothesis),
			      Form("gZOvsLumi_nom_Mu%d", muHypothesis));
  gZ0vsLumiSqrt_nom->SetNameTitle(Form("gZOvsLumiSqrt_nom_Mu%d", muHypothesis),
				  Form("gZOvsLumiSqrt_nom_Mu%d", muHypothesis));
  if (muHypothesis == 1) {
    gZ0vsLumi_nom->SetLineColor(kBlue+1);
    gZ0vsLumiSqrt_nom->SetLineColor(kBlue+1);
  }
  else {
    gZ0vsLumi_nom->SetLineColor(kRed+1);
    gZ0vsLumiSqrt_nom->SetLineColor(kRed+1);
  }
  gZ0vsLumi_nom->SetLineWidth(2);
  gZ0vsLumi_nom->GetXaxis()->SetTitle("Luminosity at 13 TeV in 2016 [fb^{-1}]");
  gZ0vsLumiSqrt_nom->SetLineWidth(2);
  gZ0vsLumiSqrt_nom->GetXaxis()
    ->SetTitle("#sqrt{L} at 13 TeV in 2016 [#sqrt{fb^{-1}}]");
  if (options.Contains("Only2016")) {
    gZ0vsLumi_nom->GetYaxis()->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
    gZ0vsLumiSqrt_nom->GetYaxis()
      ->SetTitle("Z_{0}^{Local} for 2016 data [#sigma]");
  }
  else {
    gZ0vsLumi_nom->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
    gZ0vsLumiSqrt_nom->GetYaxis()->SetTitle("Combined Z_{0}^{Local} [#sigma]");
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
  
  // Then print canvas:
  can->Print(Form("%s/plot_Z0vsLumi_Mu%d_%s_%s.eps", outputDir.Data(),
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
  can->Write();
  outputFile->Close();
  
  // Clock the toys for final print summary:
  time = clock() - time;
  std::cout << "\nExtrapolateSig: procedure concluded." << std::endl;
  printf("\tProgram required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return 0;
}
