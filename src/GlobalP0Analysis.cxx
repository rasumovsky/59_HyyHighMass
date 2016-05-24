////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: GlobalP0Analysis.cxx                                                //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/02/2016                                                          //
//                                                                            //
//  Calculates the global significance from a background-only toy MC ensemble.//
//                                                                            //
//  Macro options:                                                            //
//  - PlotAnalytic                                                            //
//  - StudyRetries
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "ToyAnalysis.h"
#include "TPolyLine.h"

/**
   -----------------------------------------------------------------------------
   Analytic function for the global significance Z global.
   @param zLocal - The local significance [in sigma] to convert to global Z.
   @param N - The trial factor for the analytic computation.
   @param alpha - Additional modification parameter (alpha=0 is nominal shape).
   @return - The global significance Z-global [in sigma].
*/
double analyticZGlobal(double zLocal, int N, double alpha) {
  return TMath::NormQuantile(TMath::Power(ROOT::Math::gaussian_cdf(zLocal*(1+alpha)), N));
}

/**
   -----------------------------------------------------------------------------
   Calculate p-value based on the Z (significance).
   @param zValue - The Z (significance).
   @return - The p-value (probability).
*/
double getPFromZ(double zValue) {
  return (1.0 - (ROOT::Math::gaussian_cdf(zValue)));
}

/**
   -----------------------------------------------------------------------------
   Calculate Z (significance) based on the p-value.
   @param p - The p-value.
   @return - The significance (# standard deviations).
*/
double getZFromP(double p) {
  return TMath::NormQuantile(1.0 - p);
}

/**
   Figure out a trial factor N from the zValue of the quantile corresponding to
   probability p.
   @param zValue - The local z value with global probability pValue.
   @param pValue - Global p-value corresponding to local significance zValue.
   @return - N, the trial factor.
*/
double getNfromQuantile(double zValue, double pValue) {
  return TMath::Log(1.0-pValue) / TMath::Log(ROOT::Math::gaussian_cdf(zValue));
}

/**
   -----------------------------------------------------------------------------
   Convert local significance to global significance using toy MC result.
   @param mappingGraph - The graph of z0_local -> z0_global.
   @param zLocal - The local significance of interest Z-local [sigma].
   @return - The global significance Z-global [sigma]
*/
double convertZLocalToGlobal(TGraphErrors *mappingGraph, double zLocal) {
  // Calculate the Z0 global value with errors:
  double xValue = 0.0; double yValue = 0.0;
  double xError = 0.0; double yError = 0.0;
  for (int i_p = 0; i_p < mappingGraph->GetN(); i_p++) {
    mappingGraph->GetPoint(i_p, xValue, yValue);
    xError = mappingGraph->GetErrorX(i_p);
    yError = mappingGraph->GetErrorY(i_p);
    if (((xValue + xError) >= zLocal) &&
	((xValue - xError) <= zLocal)) {
      break;
    }
  }
  return yValue;
}

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
  TString outputDir = Form("%s/%s/GlobalP0Analysis", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));

  // Load the toy analysis class, which does the analysis of toy MC jobs:
  TString toyAnaOptions = "";
  if (options.Contains("StudyRetries")) toyAnaOptions += "StudyRetries";
  ToyAnalysis *toyAna = new ToyAnalysis(configFile, toyAnaOptions);
  toyAna->setOutputDir(outputDir);
  toyAna->setStatHistRanges(500, 0, 20);

  std::vector<TString> fitTypes;
  fitTypes.clear(); fitTypes.push_back("0"); fitTypes.push_back("Free");
  toyAna->setFitTypes(fitTypes);
  toyAna->loadToy(0, Form("%s/%s/GlobalP0Toys/single_files/toy_mu0*", 
			  (config->getStr("MasterOutput")).Data(),
			  (config->getStr("JobName")).Data()));
  
  if (!(toyAna->areInputFilesOK())) {
    std::cout << "GlobalP0Analysis: ERROR loading toys." << std::endl;
    exit(0);
  }
    
  // Plot the toy MC nuisance parameter, global observables, and PoI:
  std::cout << "GlobalP0Analysis: Plotting parameters from toy MC" << std::endl;
  std::vector<TString> namesGlobs = toyAna->getNamesGlobalObservables();
  std::vector<TString> namesNuis = toyAna->getNamesNuisanceParameters();
  std::vector<TString> namesPars = toyAna->getNamesPoI();
  for (int i_g = 0; i_g < (int)namesGlobs.size(); i_g++) {
    if (!(namesGlobs[i_g]).Contains("gamma_stat_channel_bin")) {
      toyAna->plotHist(namesGlobs[i_g], 0);// Mu=0 toy data only
    }
  }
  for (int i_n = 0; i_n < (int)namesNuis.size(); i_n++) {
    if (!(namesNuis[i_n]).Contains("gamma_stat_channel_bin")) {
      toyAna->plotHist(namesNuis[i_n], 0);
    }
  }
  for (int i_p = 0; i_p < (int)namesPars.size(); i_p++) {
    toyAna->plotHist(namesPars[i_p], 0);
  }
  // Also plot the retries:
  if (options.Contains("StudyRetries")) toyAna->plotRetries(0);
  
  double observedZ0 = config->getNum("GlobalP0AnalysisSigma");
  
  //----------------------------------------//
  // Start the plotting!
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can",800,800);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.5, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.5);
  pad1->SetBottomMargin(0.00001);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad1->SetBorderMode(0);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();

  //----------//
  // Pad 1:
  pad1->cd();
  
  // Plot the distribution of maximum p0 values:
  TH1F *hMaxZ0 = toyAna->getStatHist("Z0", 0);
  hMaxZ0->SetLineWidth(2);
  hMaxZ0->SetLineColor(kBlue-1);
  hMaxZ0->SetFillColor(kBlue-10);
  hMaxZ0->GetYaxis()->SetTitle("Fraction of toy MC");
  hMaxZ0->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  hMaxZ0->GetXaxis()->SetTitleSize(0.07);
  hMaxZ0->GetXaxis()->SetLabelSize(0.06);
  hMaxZ0->GetYaxis()->SetTitleSize(0.07);
  hMaxZ0->GetYaxis()->SetTitleOffset(0.9);
  hMaxZ0->GetYaxis()->SetLabelSize(0.06);
  hMaxZ0->GetYaxis()->SetRangeUser(0.00001, 0.015);
  hMaxZ0->Draw("hist");
  //gPad->SetLogy();
  
  // Find median Z0:
  std::vector<double> valsZ0 = toyAna->getStatValues("Z0", 0);
  std::sort(valsZ0.begin(), valsZ0.end());
  double medianZ0 = valsZ0[(int)(((double)valsZ0.size())/2.0)];
  
  double pForMatching = getPFromZ(config->getNum("GlobalP0AnalysisZToMatch"));
  double zForMatching
    = valsZ0[(int)(((double)valsZ0.size())*(1.0-pForMatching))];
  
  // Also fit the histogram:
  
  // Single-parameter fit:
  //TF1 *fAnalytic = new TF1("dPdZ", "([0]*[1]*TMath::Power(ROOT::Math::gaussian_cdf(x), ([0]-1.0))*ROOT::Math::gaussian_pdf(x))", 

  // 2-parameter fit:
  TF1 *fAnalytic = new TF1("dPdZ", "([0]*[1]*TMath::Power(ROOT::Math::gaussian_cdf(x*(1+[2])), ([0]-1.0))*ROOT::Math::gaussian_pdf(x*(1+[2]))*(1+[2]))", 
			   hMaxZ0->GetXaxis()->GetXmin(),
			   hMaxZ0->GetXaxis()->GetXmax());

  fAnalytic->SetParameter(0, getNfromQuantile(medianZ0, 0.5));
  fAnalytic->SetParameter(1, 0.0088467);
  fAnalytic->SetParameter(2, 0.0);// NEW
  TString trialMethod = config->getStr("TrialMethod");
  if (trialMethod.Contains("FixTrial")) {
    fAnalytic->FixParameter(0, config->getNum("GlobalP0FixedN"));
  }
  else if (trialMethod.Contains("MatchTrial")) {
    fAnalytic->FixParameter(0, getNfromQuantile(zForMatching, pForMatching));
  }
  
  // Then fit and plot:
  if (options.Contains("PlotAnalytic")) {
    hMaxZ0->Fit(fAnalytic, "R");
    fAnalytic->SetLineWidth(2);
    fAnalytic->SetLineStyle(1);
    fAnalytic->SetLineColor(kBlue);
    fAnalytic->Draw("LSAME");
  }
  
  // Draw a line at the median:
  TLine *line1 = new TLine();
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->SetLineColor(kBlue-1);
  line1->DrawLine(medianZ0, hMaxZ0->GetYaxis()->GetXmin(),
		  medianZ0, hMaxZ0->GetMaximum());
  
  // Draw a line at the observed Z0Local:
  TLine *line2 = new TLine();
  line2->SetLineStyle(2);
  line2->SetLineWidth(2);
  line2->SetLineColor(kRed+1);
  line2->DrawLine(observedZ0, hMaxZ0->GetYaxis()->GetXmin(),
		  observedZ0, hMaxZ0->GetMaximum());
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.07);
  t.DrawLatex(0.20, 0.84, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.07);
  t.DrawLatex(0.32, 0.84, config->getStr("ATLASLabel"));
  t.DrawLatex(0.20, 0.76, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Legend:
  //TLegend leg1(0.2, 0.49, 0.49, 0.71);
  TLegend leg1(0.524, 0.60, 0.89, 0.91);
  leg1.SetTextFont(42); 
  leg1.SetTextSize(0.06);
  leg1.SetBorderSize(0);
  leg1.SetFillColor(0);
  if (options.Contains("PlotAnalytic")) {
    leg1.AddEntry(fAnalytic, 
		  "#Phi(z(1+#alpha))^{N-1}#Phi'(z(1+#alpha))", "l");
  }
  leg1.AddEntry(hMaxZ0, "Toy Monte Carlo", "F");
  leg1.AddEntry(line1, Form("Med. Z_{0}^{Local}=%2.1f#sigma",medianZ0), "l");
  leg1.AddEntry(line2, Form("Obs. Z_{0}^{Local}=%2.1f#sigma",observedZ0), "l");
  leg1.Draw("SAME");

  //----------//
  // Pad 2:
  
  // Create a new TGraph for global significance based on the histogram:
  pad2->cd();
  TGraphErrors *gZGlobal = new TGraphErrors();
  gZGlobal->SetNameTitle("GlobalZ0", "GlobalZ0");
  int pointIndex = 0;
  double normTotal = hMaxZ0->Integral();
  for (int i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
    double pValue = hMaxZ0->Integral(i_b, hMaxZ0->GetNbinsX());
    double zValue = getZFromP(pValue);
    if (pValue > 0 && i_b > 1) {
      double zLocalErr = hMaxZ0->GetBinWidth(i_b);
      int nToys = (int)valsZ0.size();
      double pGlobalErr = toyAna->calculateErrorPVal(pValue, nToys);
      //double pGlobalErr = toyAna->calculateErrorFromCounting(pValue, nToys);
      double zValueHi = getZFromP(pValue+pGlobalErr);
      double zGlobalErr = fabs(zValueHi - zValue);
      gZGlobal->SetPoint(pointIndex, hMaxZ0->GetBinCenter(i_b), zValue);
      gZGlobal->SetPointError(pointIndex, zLocalErr, zGlobalErr);
      pointIndex++;
    }
  }
  gZGlobal->SetLineWidth(2);
  gZGlobal->SetLineColor(kOrange+4);
  gZGlobal->SetFillColor(kOrange+1);

  gZGlobal->GetYaxis()->SetTitle("Z_{0}^{Global} [#sigma]");
  gZGlobal->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  gZGlobal->GetXaxis()->SetRangeUser(hMaxZ0->GetXaxis()->GetXmin(),
				     hMaxZ0->GetXaxis()->GetXmax());
  
  TH1F *hForAxis = new TH1F("hForAxis", "hForAxis", hMaxZ0->GetNbinsX(),
			    hMaxZ0->GetXaxis()->GetXmin(),
			    hMaxZ0->GetXaxis()->GetXmax());
  hForAxis->GetYaxis()->SetRangeUser(-2.0, 5.2);
  hForAxis->GetYaxis()->SetTitle("Z_{0}^{Global} [#sigma]");
  hForAxis->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  hForAxis->GetXaxis()->SetTitleSize(0.07);
  hForAxis->GetXaxis()->SetLabelSize(0.06);
  hForAxis->GetYaxis()->SetTitleSize(0.07);
  hForAxis->GetYaxis()->SetTitleOffset(0.9);
  hForAxis->GetYaxis()->SetLabelSize(0.06);
  hForAxis->SetLineColor(0);
  hForAxis->Draw();
  
  TLine *line3 = new TLine();
  line3->SetLineStyle(1);
  line3->SetLineWidth(2);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), 0,
		  hForAxis->GetXaxis()->GetXmax(), 0);
  line3->SetLineWidth(1);
  line3->SetLineStyle(3);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), 1,
		  hForAxis->GetXaxis()->GetXmax(), 1);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), -1,
		  hForAxis->GetXaxis()->GetXmax(), -1);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), 2,
		  hForAxis->GetXaxis()->GetXmax(), 2);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), 3,
		  hForAxis->GetXaxis()->GetXmax(), 3);
  line3->DrawLine(hForAxis->GetYaxis()->GetXmin(), 4,
		  hForAxis->GetXaxis()->GetXmax(), 4);
    
  gZGlobal->Draw("2SAME");
    
  // Now get a global significance just from the fit:
  TGraph *gAnalytic = new TGraph();
  for (int i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
    gAnalytic->SetPoint(i_b-1, hMaxZ0->GetBinCenter(i_b),
			analyticZGlobal(hMaxZ0->GetBinCenter(i_b),
					fAnalytic->GetParameter(0),
					fAnalytic->GetParameter(2)));
  }
  
  gAnalytic->SetLineWidth(2); 
  gAnalytic->SetLineColor(kBlue);
  gAnalytic->SetLineStyle(7);
  if (options.Contains("PlotAnalytic")) {
    gAnalytic->Draw("LSAME");
  }
      
  // Calculate the Z0 global value with errors:
  double xValue = 0.0; double yValue = 0.0;
  double xError = 0.0; double yError = 0.0;
  for (int i_p = 0; i_p < gZGlobal->GetN(); i_p++) {
    gZGlobal->GetPoint(i_p, xValue, yValue);
    xError = gZGlobal->GetErrorX(i_p);
    yError = gZGlobal->GetErrorY(i_p);
    if (((xValue + xError) >= observedZ0) &&
	((xValue - xError) <= observedZ0)) {
      break;
    }
  }
    
  // Create lines showin the Z0 local -> global conversion:
  Double_t xZ0Global[4] = {hMaxZ0->GetXaxis()->GetXmin(), 
			   hMaxZ0->GetXaxis()->GetXmin(),
			   xValue,
			   xValue};
  Double_t yZ0Global[4] = {yValue+yError,
			   yValue-yError,
			   yValue-yError,
			   yValue+yError};
  
  TPolyLine *lineZ0Global = new TPolyLine(4, xZ0Global, yZ0Global);
  lineZ0Global->SetFillColor(kRed+1);
  //lineZ0Global->SetFillStyle(3345);
  lineZ0Global->SetFillStyle(3245);
  lineZ0Global->SetLineColor(kRed+1);
  lineZ0Global->SetLineWidth(3);
  lineZ0Global->Draw("f");
  hForAxis->Draw("axisSAME");
  
  // Draw a line at the observed Z0Local:
  TLine *line4 = new TLine();
  line4->SetLineStyle(2);
  line4->SetLineWidth(2);
  line4->SetLineColor(kRed+1);
  line4->DrawLine(observedZ0, hMaxZ0->GetYaxis()->GetXmin(),
		  observedZ0, yValue);
  
  // Legend for Pad 2:
  TLegend leg2(0.18, 0.75, 0.59, 0.97);
  leg2.SetTextFont(42); 
  leg2.SetTextSize(0.06);
  leg2.SetBorderSize(0);
  leg2.SetFillColor(0);
  leg2.AddEntry(gZGlobal, "Toy MC #pm stat. error", "F");
  leg2.AddEntry(lineZ0Global, Form("Z_{0}^{Global}=%2.2f#sigma #pm %2.2f#sigma",
				   yValue, yError), "F");
  if (options.Contains("PlotAnalytic")) {
    leg2.AddEntry(gAnalytic, "#Phi^{-1}(#Phi(z(1+#alpha)^{N})", "l");
  }
  leg2.Draw("SAME");
  
  // Print the canvas:
  can->Print(Form("%s/plot_maxZ0_%s.eps", outputDir.Data(), anaType.Data()));
  can->Clear();
  
  // Delete pointers, close files, return:
  std::cout << "\nGlobalP0Analysis: Finished!" << std::endl;
  
  // Print the results:
  std::cout << "\tFrom toy: Z0Global( " << observedZ0 << " ) = " 
	    << yValue << " +/- " << yError << std::endl;

  std::cout << "\t From analytic function: Z0Global( " << observedZ0 << " ) = " 
	    << analyticZGlobal(observedZ0, fAnalytic->GetParameter(0), 
			       fAnalytic->GetParameter(2))
	    << std::endl;
  
  // Finally, save the TGraph containing the local -> global Z mapping:
  TFile *outZFile = new TFile(Form("%s/graph_maxZ0_%s.root", outputDir.Data(), 
				   anaType.Data()), "RECREATE");
  hMaxZ0->Write("hMaxZ0");
  gZGlobal->Write("gZGlobal");
  fAnalytic->Write("fAnalytic");
  gAnalytic->Write("gAnalytic");
  outZFile->Close();
  
  delete config;
  return 0;
}
