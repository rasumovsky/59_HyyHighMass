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
//  - PlotGauss                                                               //
//  - PlotLine                                                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "ToyAnalysis.h"
#include "TPolyLine.h"

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
   -----------------------------------------------------------------------------
   Compute the local significance to the global significance.
   @param mappingGraph - The graph of z0_local -> z0_global.
   @param zLocal - The local significance of interest.
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
  //hMaxZ0->GetYaxis()->SetRangeUser(0.0001, 0.014);
  hMaxZ0->GetYaxis()->SetRangeUser(0.0001, 0.014);
  hMaxZ0->Draw("hist");

  // Find median Z0:
  std::vector<double> valsZ0 = toyAna->getStatValues("Z0", 0);
  std::sort(valsZ0.begin(), valsZ0.end());
  double medianZ0 = valsZ0[(int)(((double)valsZ0.size())/2.0)];

  // Also fit the histogram with a Gaussian:
  TF1 *fGauss = new TF1("fGauss", "gaus", hMaxZ0->GetXaxis()->GetXmin(), 
  			hMaxZ0->GetXaxis()->GetXmax());
  TF1 *fGaussP1 = new TF1("fGaussP1", "gaus", hMaxZ0->GetXaxis()->GetXmin(), 
  			hMaxZ0->GetXaxis()->GetXmax());
  TF1 *fGaussN1 = new TF1("fGaussN1", "gaus", hMaxZ0->GetXaxis()->GetXmin(), 
  			hMaxZ0->GetXaxis()->GetXmax());
  
  if (options.Contains("PlotGauss")) {
    //fGauss->FixParameter(1, medianZ0);
    hMaxZ0->Fit(fGauss, "0");
    fGauss->SetLineWidth(2);
    fGauss->SetLineStyle(1);
    fGauss->SetLineColor(kBlue);
    fGauss->Draw("LSAME");
    
    // Then set +/-1 sigma values from fit:
    fGaussP1->SetParameter(1, fGauss->GetParameter(1));
    fGaussP1->SetParameter(2, 0.9*fGauss->GetParameter(2));
    fGaussN1->SetParameter(1, fGauss->GetParameter(1));
    fGaussN1->SetParameter(2, 1.1*fGauss->GetParameter(2));
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
  TLegend leg1(0.2, 0.49, 0.49, 0.71);
  leg1.SetTextFont(42); 
  leg1.SetTextSize(0.07);
  leg1.SetBorderSize(0);
  leg1.SetFillColor(0);
  leg1.AddEntry(line1, Form("Med. Z_{0}^{Local}=%2.1f#sigma",medianZ0), "l");
  leg1.AddEntry(line2, Form("Obs. Z_{0}^{Local}=%2.1f#sigma",observedZ0), "l");
  if (options.Contains("PlotGauss")) leg1.AddEntry(fGauss, "Gaussian", "l");
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
    
  // Also fit the graph:
  TF1 *fZGlobal = new TF1("fZGlobal", "pol1", hForAxis->GetXaxis()->GetXmin(),
  hForAxis->GetXaxis()->GetXmax());
  
  gZGlobal->Fit(fZGlobal, "0");
  fZGlobal->SetLineWidth(2);
  fZGlobal->SetLineStyle(2);
  fZGlobal->SetLineColor(1);
  
  gZGlobal->Draw("2SAME");
  if (options.Contains("PlotLine")) fZGlobal->Draw("LSAME");
  
  // Now get a global significance just from the Gaussian fit:
  TGraph *gFromGauss = new TGraph();
  TGraph *gFromGaussP1 = new TGraph();
  TGraph *gFromGaussN1 = new TGraph();
  
  double totalIntegral = fGauss->Integral(-5, hMaxZ0->GetXaxis()->GetXmax());
  double totalIntegralP1 = fGaussP1->Integral(-5,hMaxZ0->GetXaxis()->GetXmax());
  double totalIntegralN1 = fGaussN1->Integral(-5,hMaxZ0->GetXaxis()->GetXmax());
  for (int i_p = 0; i_p < gZGlobal->GetN(); i_p++) {
    double xCurr = 0.0; double yCurr = 0.0;
    gZGlobal->GetPoint(i_p, xCurr, yCurr);
    double gaussIntegral 
      = fGauss->Integral(xCurr, hMaxZ0->GetXaxis()->GetXmax());
    double gaussIntegralP1
      = fGaussP1->Integral(xCurr, hMaxZ0->GetXaxis()->GetXmax());
    double gaussIntegralN1
      = fGaussN1->Integral(xCurr, hMaxZ0->GetXaxis()->GetXmax());
    gFromGauss->SetPoint(i_p, xCurr, getZFromP(gaussIntegral/totalIntegral));
    gFromGaussP1->SetPoint(i_p, xCurr,
			   getZFromP(gaussIntegralP1/totalIntegralP1));
    gFromGaussN1->SetPoint(i_p, xCurr,
			   getZFromP(gaussIntegralN1/totalIntegralN1));
  }
  gFromGauss->SetLineWidth(2); gFromGauss->SetLineColor(kBlue);
  gFromGaussP1->SetLineWidth(1); gFromGaussP1->SetLineColor(kBlue);
  gFromGaussN1->SetLineWidth(1); gFromGaussN1->SetLineColor(kBlue);
  if (options.Contains("PlotGauss")) {
    gFromGauss->Draw("LSAME");
    gFromGaussP1->Draw("LSAME");
    gFromGaussN1->Draw("LSAME");
  }
  
  if (options.Contains("PlotLine")) {
    // Print functional form:
    TString fText = Form("Z_{0}^{Global} = %2.2f Z_{0}^{Local} + %2.2f",
			 fZGlobal->GetParameter(1), fZGlobal->GetParameter(0));
    fText.ReplaceAll("+ -","- ");
    t.DrawLatex(0.58, 0.34, fText);
    
    // Print the chi^2 probability:
    double chi2 = fZGlobal->GetChisquare();
    double probChi2 = TMath::Prob(fZGlobal->GetChisquare(), gZGlobal->GetN());
    t.DrawLatex(0.58, 0.24, Form("p(#chi^{2}) = %2.2f", probChi2));
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
  leg2.SetTextSize(0.07);
  leg2.SetBorderSize(0);
  leg2.SetFillColor(0);
  leg2.AddEntry(gZGlobal, "Toy MC #pm stat. error", "F");
  if (options.Contains("PlotLine")) {
    leg2.AddEntry(fZGlobal, "Fit to Toy MC", "l");
  }
  if (options.Contains("PlotGauss")) {
    leg2.AddEntry(gFromGauss, "Gaussian Fit", "l");
  }
  //leg2.AddEntry(lineExcl, "Upper bound on Z_{0}^{Global}", "F");
  leg2.AddEntry(lineZ0Global, Form("Z_{0}^{Global}=%2.2f#sigma #pm %2.2f#sigma",
				   yValue, yError), "F");
  leg2.Draw("SAME");
  
  // Print the canvas:
  can->Print(Form("%s/plot_maxZ0_%s.eps", outputDir.Data(), anaType.Data()));
  can->Clear();
  
  // Delete pointers, close files, return:
  std::cout << "\nGlobalP0Analysis: Finished!" << std::endl;
  
  // Print the results:
  std::cout << "\tFrom linear fit: Z0Global( " << observedZ0 << " ) = " 
	    << fZGlobal->Eval(observedZ0) << std::endl;
    
  std::cout << "\tFrom toy: Z0Global( " << observedZ0 << " ) = " 
	    << yValue << " +/- " << yError << std::endl;

  // Finally, save the TGraph containing the local -> global Z mapping:
  TFile *outZFile = new TFile(Form("%s/graph_maxZ0_%s.root", outputDir.Data(), 
				   anaType.Data()), "RECREATE");
  gZGlobal->Write();
  outZFile->Close();
  
  delete config;
  return 0;
}
