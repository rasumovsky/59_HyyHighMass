////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: SimpleSignificance.cxx                                              //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 12/05/2016                                                          //
//                                                                            //
//  Try a simple calculation of the global significance by assuming that the  //
//  grid has independent asymptotic tests.                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "TPolyLine.h"
#include "RooFitHead.h"
#include "statistics.h"

/**
   -----------------------------------------------------------------------------
   Calculate the error from toy MC statistics on calculated p-value, using 
   binomial errors.
   @param pValue - The p-value for which we are computing the error.
   @param nToys - The number of toy MC used to compute the p-value.
   @return - The binomial error on the p-value.
*/
double calculateErrorPVal(double pValue, int nToys) {
  return sqrt(pValue * (1.0 - pValue) / ((double)nToys));
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
   Implements the functional form of q0 (same as qMu).
   @param x - The value of the test statistic.
   @return - The value of the asymptotic q0 test statistic distribution.
*/
double functionQ0(double x) {
  // This corresponds to the "special case" of mu'=0
  double result = (TMath::Exp(-1*fabs(x)/2.0) / 
		   (2.0*sqrt(2.0*TMath::Pi()*fabs(x))));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Evaluate an analytic expression for the global significance using the local
   p0 probability and the number of independent trials.
   @param pLocal - The local p probability
   @param N - The number of trials.
   @return - The global p0 probabiltiy.
*/
double analyticPGlobal(double pLocal, int N) {
  return (1.0 - TMath::Power((1.0 - pLocal), N));
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void printProgressBar(int index, int total) {
  if (index%100000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

/**
   -----------------------------------------------------------------------------
   The main method randomly samples a chi^2 distribution and computes a dummy
   local and global significance. 
   @param configFile - The analysis configuration file.
   @param options - Job options. Can be "toy" or "asymptotic" or "both"
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
  TString outputDir = Form("%s/%s/SimpleSignificance", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Color palette:
  Color_t lineColors[5] = {kGreen+1, kCyan+1, kBlue+1, kMagenta+1, kRed+1};
  Color_t fillColors[5] = {kGreen-10, kCyan-10, kBlue-10, kMagenta-10, kRed-10};
  
  // Histograms and graphs:
  int const testPoints = 5;
  TH1F *hMaxZ0[testPoints];
  TH1F *hMaxp0[testPoints];
  TGraph *gZAnalytic[testPoints];
  TGraph *gPAnalytic[testPoints];
  TH1F *hForAxis;
  
  // Create the canvases and pads for drawing later:
  TCanvas *can = new TCanvas("can", "can", 800, 800);
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
  pad1->cd();
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.07);
  t.DrawLatex(0.20, 0.84, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.07);
  t.DrawLatex(0.32, 0.84, config->getStr("ATLASLabel"));
  
  // Legend:
  TLegend leg(0.2, 0.56, 0.51, 0.78);
  leg.SetTextFont(42); 
  leg.SetTextSize(0.06);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  
  // Loop over numbers of trials:
  for (int i_e = 0; i_e < 5; i_e++) {
    int nTrials = (int)(TMath::Power(10, i_e));
        
    // Create a histogram to store max significance values from each experiment:
    hMaxZ0[i_e] = new TH1F(Form("hMaxZ0_%d",i_e), Form("hMaxZ0_%d",i_e),
			   1000, 0.0, 6.0);
    hMaxp0[i_e] = new TH1F(Form("hMaxp0_%d",i_e), Form("hMaxp0_%d",i_e),
			   1000, 0.0, 1.0);
    for (int i_b = 1; i_b <= hMaxZ0[i_e]->GetNbinsX(); i_b++) {
      double zA = (hMaxZ0[i_e]->GetBinCenter(i_b) - 
		   (0.5*hMaxZ0[i_e]->GetBinWidth(i_b)));
      double zB = (hMaxZ0[i_e]->GetBinCenter(i_b) + 
		   (0.5*hMaxZ0[i_e]->GetBinWidth(i_b)));
      hMaxZ0[i_e]->SetBinContent(i_b, (analyticPGlobal(getPFromZ(zB),nTrials) -
				       analyticPGlobal(getPFromZ(zA),nTrials)));
    }
    for (int i_b = 1; i_b <= hMaxp0[i_e]->GetNbinsX(); i_b++) {
      double pA = (hMaxp0[i_e]->GetBinCenter(i_b) - 
		   (0.5*hMaxp0[i_e]->GetBinWidth(i_b)));
      double pB = (hMaxp0[i_e]->GetBinCenter(i_b) + 
		   (0.5*hMaxp0[i_e]->GetBinWidth(i_b)));
      hMaxp0[i_e]->SetBinContent(i_b, (analyticPGlobal(pB,nTrials) -
				       analyticPGlobal(pA,nTrials)));
    }
    
    // Also create axis histogram:
    if (i_e == 0) {
      hForAxis = new TH1F("hForAxis", "hForAxis", hMaxZ0[i_e]->GetNbinsX(),
			  hMaxZ0[i_e]->GetXaxis()->GetXmin(),
			  hMaxZ0[i_e]->GetXaxis()->GetXmax());
      hForAxis->GetYaxis()->SetRangeUser(-2.0, 5.2);
      hForAxis->GetYaxis()->SetTitle("Z_{0}^{Global} [#sigma]");
      hForAxis->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
      hForAxis->GetXaxis()->SetTitleSize(0.07);
      hForAxis->GetXaxis()->SetLabelSize(0.06);
      hForAxis->GetYaxis()->SetTitleSize(0.07);
      hForAxis->GetYaxis()->SetTitleOffset(0.9);
      hForAxis->GetYaxis()->SetLabelSize(0.06);
      hForAxis->SetLineColor(0);
    }

    // Pad 1: Plot the distribution of maximum p0 values:
    pad1->cd();
    
    hMaxZ0[i_e]->Scale(1.0 / hMaxZ0[i_e]->Integral());
    hMaxZ0[i_e]->SetLineWidth(2);
    hMaxZ0[i_e]->SetLineColor(lineColors[i_e]);
    hMaxZ0[i_e]->SetFillColor(fillColors[i_e]);
    hMaxZ0[i_e]->GetYaxis()->SetTitle("Fraction of toy MC");
    hMaxZ0[i_e]->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
    hMaxZ0[i_e]->GetXaxis()->SetTitleSize(0.07);
    hMaxZ0[i_e]->GetXaxis()->SetLabelSize(0.06);
    hMaxZ0[i_e]->GetYaxis()->SetTitleSize(0.07);
    hMaxZ0[i_e]->GetYaxis()->SetTitleOffset(0.9);
    hMaxZ0[i_e]->GetYaxis()->SetLabelSize(0.06);
    if (i_e == 0) hMaxZ0[i_e]->Draw("hist");
    else hMaxZ0[i_e]->Draw("histSAME");
    leg.AddEntry(hMaxZ0[i_e], Form("N = %d", nTrials), "F");
    
    // Pad 2:
    pad2->cd();
    
    // Analytic Z 
    gZAnalytic[i_e] = new TGraph();
    gZAnalytic[i_e]->SetNameTitle(Form("zAnalytic_%dtrials", nTrials),
				  Form("zAnalytic_%dtrials", nTrials));
    int pointIndex = 0;
    for (int i_b = 1; i_b <= hMaxZ0[i_e]->GetNbinsX(); i_b++) {
      double pValue 
	= analyticPGlobal(getPFromZ(hMaxZ0[i_e]->GetBinCenter(i_b)), nTrials);
      double zValue = getZFromP(pValue);
      gZAnalytic[i_e]->SetPoint(pointIndex, hMaxZ0[i_e]->GetBinCenter(i_b),
				zValue);
      pointIndex++;
    }
    gZAnalytic[i_e]->SetLineWidth(2);
    gZAnalytic[i_e]->SetLineColor(lineColors[i_e]);
    gZAnalytic[i_e]->GetXaxis()
      ->SetRangeUser(hMaxZ0[i_e]->GetXaxis()->GetXmin(),
		     hMaxZ0[i_e]->GetXaxis()->GetXmax());
    
    // Analytic P
    gPAnalytic[i_e] = new TGraph();
    gZAnalytic[i_e]->SetNameTitle(Form("pAnalytic_%dtrials", nTrials),
				  Form("pAnalytic_%dtrials", nTrials));
    pointIndex = 0;
    for (int i_b = 1; i_b <= hMaxp0[i_e]->GetNbinsX(); i_b++) {
      double pValue = analyticPGlobal(hMaxp0[i_e]->GetBinCenter(i_b), nTrials);
      gPAnalytic[i_e]->SetPoint(pointIndex, hMaxp0[i_e]->GetBinCenter(i_b), 
				pValue);
      pointIndex++;
    }
    gPAnalytic[i_e]->SetLineWidth(2);
    gPAnalytic[i_e]->SetLineColor(lineColors[i_e]);
    gPAnalytic[i_e]->GetXaxis()
      ->SetRangeUser(hMaxp0[i_e]->GetXaxis()->GetXmin(),
		     hMaxp0[i_e]->GetXaxis()->GetXmax());
    
    if (i_e == 0) hForAxis->Draw();
    gZAnalytic[i_e]->Draw("LSAME");
  }

  // Print the canvas:
  pad2->cd();
  hForAxis->Draw("axisSAME");
  pad1->cd();
  leg.Draw("SAME");
  
  can->Print(Form("%s/plot_analyticZ0.eps", outputDir.Data()));
  
  delete config;
  return 0;
}
