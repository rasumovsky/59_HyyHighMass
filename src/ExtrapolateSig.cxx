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
  double yValues[100] = {0.0};
  double yErrorHi[100] = {0.0};
  double yErrorLo[100] = {0.0};
  
  //----------------------------------------//
  // Loop over luminosities to get expectation.
  std::cout << "ExtrapolateSig: Looping over lumi for extraploation" << endl;
  for (int i_l = 0; i_l < (int)lumiValues.size(); i_l++) {
    std::cout << "ExtrapolateSig: Lumi point" << i_l << " of " 
	      << lumiValues.size() << std::endl;
    
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
      
      // For a maximum-likelihood fit, set the parameters of interest using mu=1
      // snapshot. The mu=1 fit params of interest are used even for the mu=0 
      // parameters of interest because the background-only fits and toys should
      // still have the spurious signal randomized and fitted at that mass.
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
	    else {
	      mapPoIMu0[listPoI[i_p]] = snapParMu1->getVal();
	      mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
	    }
	  }
	}
      }
      
      double lumiScaleFactor = (lumiValues[i_l] / 
				config->getNum("AnalysisLuminosity"));
      // If we are assuming no signal in 2016 data, then the amount of signal
      // should be constant (only what we saw in 2015):
      if (muHypothesis == 0) {
	mapPoIMu1[poiForNorm] = (mapPoIMu1[poiForNorm] / lumiScaleFactor);
      }
      testStat->scaleAsimovData(lumiScaleFactor, 
				config->getStrV("ParamsToScale"));
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
      // Total uminosity - 2015 luminosity = 2016 luminosity:
      xValues[i_l]
	= ((lumiValues[i_l] - config->getNum("AnalysisLuminosity")) / 1000.0);
      if (i_e == 0) { // Nominal signal
	yValues[i_l] = currZ0;
	std::cout << "ExtrapSig: Z0 [" << xValues[i_l] << " ] = " 
		  << yValues[i_l] << std::endl;
      }
      else if (i_e == 1) { // Signal + fit error
	yErrorHi[i_l] = currZ0;
	std::cout << "ExtrapSig: Z0_hi [" << xValues[i_l] << " ] = " 
		  << yErrorHi[i_l] << std::endl;
      }
      else if (i_e == -1) { // Signal - fit error
	yErrorLo[i_l] = currZ0;
	std::cout << "ExtrapSig: Z0_lo [" << xValues[i_l] << " ] = " 
		  << yErrorLo[i_l] << std::endl;
      }
      
      // Close the input file before the loop repeats:
      inputFile.Close();
      delete testStat;
      delete workspace; 
    }
  }
  
  //----------------------------------------//
  // Print TGraph of sensitivity:
  
  // First re-format the y-errors:
  for (int i_l = 0; i_l < (int)lumiValues.size(); i_l++) {
    yErrorHi[i_l] = fabs(yErrorHi[i_l] - yValues[i_l]);
    yErrorLo[i_l] = fabs(yValues[i_l] - yErrorLo[i_l]);
  }
  
  // Create a graph: double xValues[100] = {0.0};
  TGraphAsymmErrors *gZ0vsLumi_err
    = new TGraphAsymmErrors((int)lumiValues.size(), xValues, yValues, 0, 0,
			    yErrorLo, yErrorHi);
  gZ0vsLumi_err->SetNameTitle(Form("gZOvsLumi_err_Mu%d", muHypothesis),
			      Form("gZOvsLumi_err_Mu%d", muHypothesis));
  if (muHypothesis == 1) {
    gZ0vsLumi_err->SetFillColor(kBlue+1);
    gZ0vsLumi_err->SetFillStyle(3245);
  }
  else {
    gZ0vsLumi_err->SetFillColor(kRed+1);
    gZ0vsLumi_err->SetFillStyle(3254);
  }
  gZ0vsLumi_err->GetXaxis()->SetTitle("Luminosity at 13 TeV in 2016 [fb^{-1}]");
  gZ0vsLumi_err->GetYaxis()->SetTitle("z_{0}^{Local} [#sigma]");
  
  TGraph *gZ0vsLumi_nom = new TGraph((int)lumiValues.size(), xValues, yValues);
  gZ0vsLumi_nom->SetNameTitle(Form("gZOvsLumi_nom_Mu%d", muHypothesis),
			      Form("gZOvsLumi_nom_Mu%d", muHypothesis));
  if (muHypothesis == 1) gZ0vsLumi_nom->SetLineColor(kBlue+1);
  else gZ0vsLumi_nom->SetLineColor(kRed+1);
  gZ0vsLumi_nom->SetLineWidth(2);
  gZ0vsLumi_nom->GetXaxis()->SetTitle("Luminosity at 13 TeV in 2016 [fb^{-1}]");
  gZ0vsLumi_nom->GetYaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  
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
  can->Print(Form("%s/plot_Z0vsLumi_Mu%d_%s.eps", outputDir.Data(),
		  muHypothesis, (config->getStr("AnalysisType").Data())));
  
  // Also save a ROOT file:
  TFile *outputFile = new TFile(Form("%s/sigExtrap_Mu%d_%s.root",
				     outputDir.Data(), muHypothesis,
				     (config->getStr("AnalysisType").Data())),
				"RECREATE");
  gZ0vsLumi_err->Write();
  gZ0vsLumi_nom->Write();
  outputFile->Close();
  
  // Clock the toys for final print summary:
  time = clock() - time;
  std::cout << "\nExtrapolateSig: procedure concluded." << std::endl;
  printf("\tProgram required %d clock cycles (%f seconds).\n\n",
	 (int)time, ((float)time/CLOCKS_PER_SEC));
  return 0;
}
