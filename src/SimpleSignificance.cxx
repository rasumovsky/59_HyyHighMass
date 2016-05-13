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
  
  // Choose the number of experiments and independent trials in each experiment:
  int nExperiments = config->getInt("SimpleZNumExperiments");
  int nTrialsPerExperiment = config->getInt("SimpleZTrialsPerExp");
  
  // Create asymptotic 1/2 chi^2 distribution histogram for sampling:
  int nBinsQ = 1000;
  int qMin = -50;
  int qMax = 50;
  TH1F *qAsymptotic = new TH1F("qAsymptotic","qAsymptotic",nBinsQ,qMin,qMax);
  for (int i_b = 1; i_b <= qAsymptotic->GetNbinsX(); i_b++) {
    double qProbability = functionQ0(qAsymptotic->GetBinCenter(i_b));
    qAsymptotic->SetBinContent(i_b, qProbability);
  }
  qAsymptotic->Scale(1.0 / qAsymptotic->Integral());
  
  // Create a histogram to store max. significance values from each experiment:
  TH1F *hMaxZ0 = new TH1F(Form("hMaxZ0_%dtrials",nTrialsPerExperiment),
			  Form("hMaxZ0_%dtrials",nTrialsPerExperiment),
			  1000, 0.0, 6.0);
  
  // Analytic form of distribution:
  if (config->getBool("SimpleZDoAnalytic")) {
    for (int i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
      double zA = hMaxZ0->GetBinCenter(i_b) - (0.5*hMaxZ0->GetBinWidth(i_b));
      double zB = hMaxZ0->GetBinCenter(i_b) + (0.5*hMaxZ0->GetBinWidth(i_b));
      double currentDP = (analyticPGlobal(getPFromZ(zB),nTrialsPerExperiment) -
			  analyticPGlobal(getPFromZ(zA),nTrialsPerExperiment));
      hMaxZ0->SetBinContent(i_b, currentDP);
    }
  }
  
  // Loop over the number of "experiments":
  else {
    std::cout << "SimpleSignificance: Looping over " << nExperiments 
	      << " experiments." << std::endl;
    for (int i_e = 0; i_e < nExperiments; i_e++) {
      
      // For each experiment, keep track of the max Z0 from all the trials:
      double experimentMaxZ0 = 0.0;
      
      // Loop over the number of independent trials:
      for (int i_t = 0; i_t < nTrialsPerExperiment; i_t++) {
	printProgressBar(i_e*nTrialsPerExperiment+i_t,
			 nExperiments*nTrialsPerExperiment);
	
	// Generate a single value of q:
	double randQ0 = qAsymptotic->GetRandom();
	
	// Convert to significance:
	double currZ0 = (randQ0 > 0) ? 
	  sqrt(randQ0) : (-1.0 * sqrt(fabs(randQ0)));
	
	// Store if it is the largest of all the trials:
	if (currZ0 > experimentMaxZ0) experimentMaxZ0 = currZ0;
      }
      hMaxZ0->Fill(experimentMaxZ0);
    }
  }

  // Then extract the global significance:
  double zLocalObs = config->getNum("GlobalP0AnalysisSigma");
  int i_b = 0; 
  for (i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
    double binMin = hMaxZ0->GetBinCenter(i_b) - (0.5*hMaxZ0->GetBinWidth(i_b));
    double binMax = hMaxZ0->GetBinCenter(i_b) + (0.5*hMaxZ0->GetBinWidth(i_b));
    if (binMax >= zLocalObs && binMin < zLocalObs) {
      break;
    }
  }
  double pGlobalObs = (hMaxZ0->Integral(i_b, hMaxZ0->GetNbinsX()) / 
		       hMaxZ0->Integral());
  double zGlobalObs = getZFromP(pGlobalObs);

  // Draw asymptotic q0 (really r0) distribution.
  TCanvas *can1 = new TCanvas("can", "can");
  can1->cd();
  gPad->SetLogy();
  qAsymptotic->GetXaxis()->SetTitle("r_{0}");
  double GeVPerBin = ((double)(qMax - qMin) / (double)nBinsQ);
  qAsymptotic->GetYaxis()->SetTitle(Form("Entries / %2.3f",GeVPerBin));
  qAsymptotic->SetLineColor(kRed+1);
  qAsymptotic->Draw();
  can1->Print(Form("%s/q0Asymptotic.eps", outputDir.Data()));
  can1->Clear();

  //----------------------------------------//
  // Start plotting:
  TCanvas *can = new TCanvas(Form("can_%dtrials",nTrialsPerExperiment),
			     Form("can_%dtrials",nTrialsPerExperiment),
			     800,800);
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
  // Pad 1: Plot the distribution of maximum p0 values:
  pad1->cd();
  
  hMaxZ0->Scale(1.0 / hMaxZ0->Integral());
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
  hMaxZ0->Draw("hist");
  
  // Find median Z0:
  double medianZ0 = 0.0;
  for (i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
    if ((hMaxZ0->Integral(i_b,hMaxZ0->GetNbinsX()) / hMaxZ0->Integral()) < 0.5){
      medianZ0 = hMaxZ0->GetBinCenter(i_b);
      break;
    }
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
  line2->DrawLine(zLocalObs, hMaxZ0->GetYaxis()->GetXmin(),
		  zLocalObs, hMaxZ0->GetMaximum());
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.07);
  t.DrawLatex(0.20, 0.84, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.07);
  t.DrawLatex(0.32, 0.84, config->getStr("ATLASLabel"));
    
  // Legend:
  TLegend leg1(0.2, 0.56, 0.51, 0.78);
  leg1.SetTextFont(42); 
  leg1.SetTextSize(0.06);
  leg1.SetBorderSize(0);
  leg1.SetFillColor(0);
  leg1.AddEntry(line1, Form("Med. Z_{0}^{Local}=%2.1f#sigma",medianZ0), "l");
  leg1.AddEntry(line2, Form("Obs. Z_{0}^{Local}=%2.1f#sigma",zLocalObs), "l");
  leg1.Draw("SAME");
  
  //----------//
  // Pad 2:
  
  // Create a new TGraph for global significance based on the histogram:
  pad2->cd();
  TGraphErrors *gZGlobal = new TGraphErrors();
  gZGlobal->SetNameTitle(Form("GlobalZ0_%dtrials",nTrialsPerExperiment),
			 Form("GlobalZ0_%dtrials",nTrialsPerExperiment));
  TGraph *gZAnalytic = new TGraph();
  gZAnalytic->SetNameTitle(Form("GlobalZ0Analytic_%dtrials", 
				nTrialsPerExperiment),
			   Form("GlobalZ0Analytic_%dtrials",
				nTrialsPerExperiment));
  
  int pointIndex = 0;
  double normTotal = hMaxZ0->Integral();
  for (int i_b = 1; i_b <= hMaxZ0->GetNbinsX(); i_b++) {
    
    // Analytic form or toy form of distribution:
    double pValue = 0.0;
    if (config->getBool("SimpleZDoAnalytic")) {
      pValue = analyticPGlobal(getPFromZ(hMaxZ0->GetBinCenter(i_b)), 
			       nTrialsPerExperiment);
    }
    else pValue = hMaxZ0->Integral(i_b, hMaxZ0->GetNbinsX());
    double zValue = getZFromP(pValue);
    
    if (pValue > 0 && i_b > 1) {
      double zLocalErr = hMaxZ0->GetBinWidth(i_b);
      double pGlobalErr = calculateErrorPVal(pGlobalObs, nExperiments);
      double zValueHi = getZFromP(pValue+pGlobalErr);
      double zGlobalErr = fabs(zValueHi - zValue);
      gZGlobal->SetPoint(pointIndex, hMaxZ0->GetBinCenter(i_b), zValue);
      gZGlobal->SetPointError(pointIndex, zLocalErr, zGlobalErr);
      gZAnalytic->SetPoint(pointIndex, hMaxZ0->GetBinCenter(i_b), zValue);
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

  gZAnalytic->SetLineWidth(2);
  gZAnalytic->SetLineColor(kBlue-1);
  gZAnalytic->GetYaxis()->SetTitle("Z_{0}^{Global} [#sigma]");
  gZAnalytic->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  gZAnalytic->GetXaxis()->SetRangeUser(hMaxZ0->GetXaxis()->GetXmin(),
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

  double xValue = 0.0; double yValue = 0.0;
  double xError = 0.0; double yError = 0.0;
  
  if (config->getBool("SimpleZDoAnalytic")) {
    gZAnalytic->Draw("LSAME");
    xValue = zLocalObs;
    yValue = gZAnalytic->Eval(xValue);
  }
  else {
    gZGlobal->Draw("2SAME");
    // Calculate the Z0 global value with errors:
    for (int i_p = 0; i_p < gZGlobal->GetN(); i_p++) {
      gZGlobal->GetPoint(i_p, xValue, yValue);
      xError = gZGlobal->GetErrorX(i_p);
      yError = gZGlobal->GetErrorY(i_p);
      if (((xValue + xError) >= zLocalObs) &&
	  ((xValue - xError) <= zLocalObs)) {
	break;
      }
    }
  }
  
  hForAxis->Draw("axisSAME");
  
  // Draw a line at the observed Z0Local:
  TLine *line4 = new TLine();
  line4->SetLineStyle(2);
  line4->SetLineWidth(2);
  line4->SetLineColor(kRed+1);
  line4->DrawLine(zLocalObs, hMaxZ0->GetYaxis()->GetXmin(),
		  zLocalObs, yValue);
 
  // Draw a line at the observed Z0Global:
  TLine *line5 = new TLine();
  line5->SetLineStyle(2);
  line5->SetLineWidth(2);
  line5->SetLineColor(kRed+1);
  line5->DrawLine(hMaxZ0->GetXaxis()->GetXmin(), yValue, zLocalObs, yValue);
  
  // Legend for Pad 2:
  TLegend leg2(0.18, 0.77, 0.59, 0.97);
  leg2.SetTextFont(42); 
  leg2.SetTextSize(0.06);
  leg2.SetBorderSize(0);
  leg2.SetFillColor(0);
  if (config->getBool("SimpleZDoAnalytic")) {
    leg2.AddEntry(line5, Form("Z_{0}^{Global}=%2.2f#sigma", yValue), "l");
  }
  else {
    leg2.AddEntry(gZGlobal, "Toy MC #pm stat. error", "F");
    leg2.AddEntry(line5, Form("Z_{0}^{Global}=%2.2f#sigma #pm %2.2f#sigma",
			      yValue, yError), "l");
  }
  leg2.Draw("SAME");
  
  // Print the canvas:
  can->Print(Form("%s/plot_maxZ0_%s_%dtrials.eps", outputDir.Data(), 
		  (config->getStr("AnalysisType")).Data(), 
		  nTrialsPerExperiment));
  can->Clear();
  
  // Print the results:
  std::cout << "\nSimpleSignificance: " << std::endl;
  std::cout << "\tmedian expected z0 = " << medianZ0 << std::endl;
  std::cout << "\tp0 global = " << pGlobalObs << std::endl;
  std::cout << "\tz0 global = " << yValue << " +/- " << yError << std::endl;
  
  // Finally, save the TGraph containing the local -> global Z mapping:
  TFile *outZFile = new TFile(Form("%s/graph_maxZ0_%dtrials.root",
				   outputDir.Data(), nTrialsPerExperiment),
			      "RECREATE");
  can->Write();
  if (config->getBool("SimpleZDoAnalytic")) gZAnalytic->Write();
  else gZGlobal->Write();
  hMaxZ0->Write();
  outZFile->Close();
  
  delete config;
  return 0;
}
