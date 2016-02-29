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
//  - "New"        Calculate everything from scratch.                         //
//  - "FromFile"   Load CL values from file.                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "ToyAnalysis.h"

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
  
  /*
  // Open the workspace:
  TFile wsFile(config->getStr("WorkspaceFile"), "read");
  RooWorkspace *workspace
    = (RooWorkspace*)wsFile.Get(config->getStr("WorkspaceName"));
  */

  ToyAnalysis *toyAna = new ToyAnalysis(configFile, "None");

  toyAna->setOutputDir(Form("%s/%s/ToyAnalysis", 
			  (config->getStr("MasterOutput")).Data(),
			  jobName.Data()));
  
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
    
  /*
  // Get the asymptotic test statistic distribution:
  getAsymptoticForm("QMu");// THIS SHOULD BE GENERALIZED!!!
  
    // Plot the results:
    plotProfiledMu();
    plotTestStat("QMu");
    plotTestStat("Q0");
    plotTestStatComparison("QMu");
    plotTestStatComparison("Q0");
    
    // Then plot the nuis, globs, and other parameters:
    for (int i_g = 0; i_g < (int)m_namesGlobs.size(); i_g++) {
      plotHist(m_namesGlobs[i_g], 0);
      plotHist(m_namesGlobs[i_g], 1);
    }
    for (int i_n = 0; i_n < (int)m_namesNuis.size(); i_n++) {
      plotHist(m_namesNuis[i_n], 0);
      plotHist(m_namesNuis[i_n], 1);
    }
    for (int i_p = 0; i_p < (int)m_namesPars.size(); i_p++) {
      plotHist(m_namesPars[i_p], 0);
      plotHist(m_namesPars[i_p], 1);
    }
  }
  */
  
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
  hMaxZ0->GetXaxis()->SetTitle("Z_{0}^{Local}");
  hMaxZ0->Draw("hist");
  
  // Draw a line at the me
  std::vector<double> valsZ0 = toyAna->getStatValues("Z0", 0);
  std::sort(valsZ0.begin(), valsZ0.end());
  double medianZ0 = valsZ0[(int)(((double)valsZ0.size())/2.0)];
  
  TLine *line1 = new TLine();
  line1->SetLineStyle(2);
  line1->SetLineWidth(3);
  line1->SetLineColor(kRed);
  line1->DrawLine(medianZ0, hMaxZ0->GetYaxis()->GetXmin(),
		  medianZ0, hMaxZ0->GetMaximum());
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.07);
  t.DrawLatex(0.62, 0.84, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.07);
  t.DrawLatex(0.74, 0.84, config->getStr("ATLASLabel"));
  t.DrawLatex(0.62, 0.76, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (config->getNum("AnalysisLuminosity")/1000.0)));
  
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
      double zValueHi = getZFromP(pValue+pGlobalErr);
      double zGlobalErr = fabs(zValueHi - zValue);
      gZGlobal->SetPoint(pointIndex, hMaxZ0->GetBinCenter(i_b), zValue);
      gZGlobal->SetPointError(pointIndex, zLocalErr, zGlobalErr);
      pointIndex++;
    }
  }
  gZGlobal->SetLineWidth(2);
  gZGlobal->SetLineColor(kRed+1);
  gZGlobal->SetFillColor(kRed+1);

  gZGlobal->GetYaxis()->SetTitle("Z_{0}^{Global}");
  gZGlobal->GetXaxis()->SetTitle("Z_{0}^{Local}");
  gZGlobal->GetXaxis()->SetRangeUser(hMaxZ0->GetXaxis()->GetXmin(),
				     hMaxZ0->GetXaxis()->GetXmax());
  
  TH1F *hForAxis = new TH1F("hForAxis", "hForAxis", hMaxZ0->GetNbinsX(),
			    hMaxZ0->GetXaxis()->GetXmin(),
			    hMaxZ0->GetXaxis()->GetXmax());
  hForAxis->GetYaxis()->SetRangeUser(-2.0, 5.2);
  hForAxis->GetYaxis()->SetTitle("Z_{0}^{Global}");
  hForAxis->GetXaxis()->SetTitle("Z_{0}^{Local}");
  hForAxis->SetLineColor(0);
  hForAxis->Draw();
  
  TLine *line2 = new TLine();
  line2->SetLineStyle(2);
  line2->SetLineWidth(3);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), 0,
		  hForAxis->GetXaxis()->GetXmax(), 0);
  line2->SetLineWidth(1);
  line2->SetLineStyle(3);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), 1,
		  hForAxis->GetXaxis()->GetXmax(), 1);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), -1,
		  hForAxis->GetXaxis()->GetXmax(), -1);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), 2,
		  hForAxis->GetXaxis()->GetXmax(), 2);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), 3,
		  hForAxis->GetXaxis()->GetXmax(), 3);
  line2->DrawLine(hForAxis->GetYaxis()->GetXmin(), 4,
		  hForAxis->GetXaxis()->GetXmax(), 4);
    
  // Also fit the graph:
  TF1 *fZGlobal = new TF1("fZGlobal", "pol1", hForAxis->GetXaxis()->GetXmin(), 
			  hForAxis->GetXaxis()->GetXmax());
  gZGlobal->Fit(fZGlobal);
  fZGlobal->SetLineWidth(2);
  fZGlobal->SetLineColor(1);
  fZGlobal->Draw("LSAME");
  
  //gZGlobal->Draw("LSAME");
  gZGlobal->Draw("2SAME");
  //gZGlobal->Draw("E1SAME");
  
  TString fText = Form("Z_{0}^{Global} = %2.2f Z_{0}^{Local} + %2.2f",
		       fZGlobal->GetParameter(1), fZGlobal->GetParameter(0));
  fText.ReplaceAll("+ -","- ");
  
  t.DrawLatex(0.2, 0.78, fText);
  
  double chi2 = fZGlobal->GetChisquare();
  double probChi2 = TMath::Prob(fZGlobal->GetChisquare(), gZGlobal->GetN());
  t.DrawLatex(0.2, 0.68, Form("p(#chi^{2}) = %2.2f",probChi2));
  
  // Legend:
  /*
  TLegend leg(0.54,0.20,0.79,0.38);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!config->getBool("DoBlind")) leg.AddEntry(gCLObs_toy,"Obs. limit","l");
  leg.AddEntry(gCLExp_toy,"Exp. limit","l");
  leg.AddEntry(gCLExp_toy_1s,"Exp. limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_toy_2s,"Exp. limit #pm2#sigma_{exp}","F");
  
  // Plotting options:
  if (options.Contains("toy")) {
    gCLExp_toy->Draw("AL");
    gCLExp_toy_2s->Draw("3same");
    gCLExp_toy_1s->Draw("3same");
    gCLExp_toy->Draw("LSAME");
    if (!config->getBool("DoBlind")) gCLObs_toy->Draw("LSAME");
  }
  gPad->RedrawAxis();
  leg.Draw("SAME");
  */
    
  // Print the canvas:
  can->Print(Form("%s/plot_maxZ0_%s.eps", outputDir.Data(), anaType.Data()));
  can->Clear();
   
  // Delete pointers, close files, return:
  std::cout << "GlobalP0Analysis: Finished!" << std::endl;
  //delete line;
  //delete pad1;
  //delete pad2;
  //delete can;
  //delete gZGlobal;
  //delete hMaxZ0;
  //delete workspace;
  delete config;
  //wsFile.Close();
  return 0;
}
