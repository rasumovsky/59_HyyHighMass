////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotGlobalP0.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 20/06/2016                                                          //
//                                                                            //
//  This program reads a list of local p0 values and plots the global p0      //
//  using the parameters stored in the config file.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "RooFitHead.h"
#include "statistics.h"

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
   Calculate p-value based on the Z (significance).
   @param zValue - The Z (significance).
   @return - The p-value (probability).
*/
double getPFromZ(double zValue) {
  return (1.0 - (ROOT::Math::gaussian_cdf(zValue)));
}

/**
   -----------------------------------------------------------------------------
   Analytic function for the global significance Z global.
   @param zLocal - The local significance [in sigma] to convert to global Z.
   @param N - The trial factor for the analytic computation.
   @param alpha - Additional modification parameter (alpha=0 is nominal shape).
   @return - The global significance Z-global [in sigma].
*/
float analyticZGlobal(double zLocal, int N, double alpha) {
  return TMath::NormQuantile(TMath::Power(ROOT::Math::gaussian_cdf(zLocal*(1+alpha)), N));
}

/**
   -----------------------------------------------------------------------------
   The main method loads a text file with local p0 values and converts it to
   global significance before plotting.
   @param configFile - The analysis configuration file.
   @param options - Job options. Can be "TestAlpha" "PlotP0"
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 4) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <p0file> <width>"
	      << std::endl;
    exit(0);
  }
  TString configFile = argv[1];
  TString inputFileP0 = argv[2];
  double width = atof(argv[3]);
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);
  TString outputDir = Form("%s/%s/PlotGlobalP0", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Declare the graph:
  TGraph *graphP0 = new TGraph();
  graphP0->SetLineColor(kBlack);
  graphP0->SetLineWidth(2);
  
  // Load p0 local values from file:
  int pointIndex = 0;
  double inMass=0.0; double inP0_1=0.0; double inP0_2=0.0;
  double massMin=0.0; double massMax=0.0;
  std::ifstream inputFile(inputFileP0);
  while (inputFile >> inMass >> inP0_1 >> inP0_2) {
    double zLocal = getZFromP(inP0_1);
    double zGlobal = analyticZGlobal(zLocal,
				     config->getNum("AnalyticZGlobal_N"), 
				     config->getNum("AnalyticZGlobal_alpha"));
    double pGlobal = getPFromZ(zGlobal);
    graphP0->SetPoint(pointIndex, inMass, pGlobal);
    if (pointIndex == 0) { massMin = inMass; massMax = inMass; }
    else if (inMass > massMax) massMax = inMass;
    else if (inMass < massMin) massMin = inMass;
    pointIndex++;
  }
  inputFile.close();
  
  // Create the canvases and pads for drawing later:
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  TH1F *hForAxis = new TH1F("hForAxis", "hForAxis", 2, massMin, massMax);
  hForAxis->GetYaxis()->SetTitle("p_{0}^{Global}");
  if ((config->getStr("AnalysisType")).Contains("Graviton")) {
    hForAxis->GetXaxis()->SetTitle("m_{G*} [GeV]");
  }
  else hForAxis->GetXaxis()->SetTitle("m_{X} [GeV]");
  hForAxis->GetYaxis()->SetRangeUser(0.0001, 1.0);
  hForAxis->Draw("axis");
  gPad->SetLogy();
  graphP0->Draw("LSAME");
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.2, 0.27, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.27, config->getStr("ATLASLabel"));
  t.DrawLatex(0.2, 0.21, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			      (config->getNum("AnalysisLuminosity")/1000.0)));
  if ((config->getStr("AnalysisType")).Contains("Scalar")) {
    t.DrawLatex(0.6, 0.27, "Spin-0 Selection");
    t.DrawLatex(0.6, 0.21,
		Form("X#rightarrow#gamma#gamma, #Gamma/m_{X}=%2.2f", width));
  }
  else if ((config->getStr("AnalysisType")).Contains("GravitonLoose")) {
    t.DrawLatex(0.6, 0.27, "Spin-2 Loose Iso.");
    t.DrawLatex(0.6, 0.21,
		Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{Pl}=%2.2f",
		     width));
  }
  else {
    t.DrawLatex(0.6, 0.27, "Spin-2 Selection");
    t.DrawLatex(0.6, 0.21,
		Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{Pl}=%2.2f",
		     width));
  }
  
  // Significance lines and text:
  TLatex sigma; sigma.SetTextColor(kRed+1);
  sigma.SetTextFont(42); sigma.SetTextSize(0.04);
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kRed+1);
  double sigmaVals[5] = {0.5,0.15865,0.02275,0.001349,0.000032};
  for (int i_s = 0; i_s < 4; i_s++) {
    double sigmaXPos = massMax - (0.07*(massMax - massMin));
    line->DrawLine(massMin, sigmaVals[i_s], massMax, sigmaVals[i_s]);
    sigma.DrawLatex(sigmaXPos, 1.1*sigmaVals[i_s], Form("%d#sigma",i_s));
  }
  
  // Save the canvas:
  can->Print(Form("%s/plot_global_p0_width%2.2f.eps", outputDir.Data(), width));
  
  
  delete config;
  return 0;
}
