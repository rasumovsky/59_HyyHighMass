////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CompareExtrap.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/04/2016                                                          //
//                                                                            //
//  This program compares the histograms of StudyData.cxx for two datasets.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"

TString m_outputDir;

/**
   -----------------------------------------------------------------------------
*/
void compareHistograms(TString name, TH1F *h1, TH1F *h2) {
    
  // Create a canvas with two pads (one main plot, one subtraction plot)
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
  pad1->SetLogy();
  
  h1->SetLineColor(kBlue+1);
  h1->SetMarkerColor(kBlue+1);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(0.9);
  h1->GetYaxis()->SetLabelSize(0.06);
  h1->Draw("E1");

  h2->SetLineColor(kRed+1);
  h2->SetMarkerColor(kRed+1);
  h2->Draw("E1SAME");
  h1->Draw("E1SAME");
  
  // Then make ratio:
  pad2->cd();
  
  TH1F *hRatio = new TH1F(Form("hRatio_%s",name.Data()),
			  Form("hRatio_%s",name.Data()), h1->GetNbinsX(), 
			  h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
  for (int i_b = 1; i_b <= h1->GetNbinsX(); i_b++) {
    if (h1->GetBinContent(i_b) > 0 && h2->GetBinContent(i_b) > 0) {
      hRatio->SetBinContent(i_b, 
			    (h1->GetBinContent(i_b)/h2->GetBinContent(i_b)));
      double error = (hRatio->GetBinContent(i_b) *
		      sqrt((h1->GetBinError(i_b) / h1->GetBinContent(i_b)) *
			   (h1->GetBinError(i_b) / h1->GetBinContent(i_b)) +
			   (h2->GetBinError(i_b) / h2->GetBinContent(i_b)) * 
			   (h2->GetBinError(i_b) / h2->GetBinContent(i_b))));
      hRatio->SetBinError(i_b, error);
    }
    else {
      hRatio->SetBinContent(i_b, 0);
      hRatio->SetBinError(i_b, 0);
    }
    
  }
  hRatio->SetLineColor(kBlue+1);
  hRatio->SetMarkerColor(kBlue+1);
  hRatio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  hRatio->GetYaxis()->SetTitle("Ratio");
  hRatio->GetXaxis()->SetTitleSize(0.07);
  hRatio->GetXaxis()->SetLabelSize(0.06);
  hRatio->GetYaxis()->SetTitleSize(0.07);
  hRatio->GetYaxis()->SetTitleOffset(0.9);
  hRatio->GetYaxis()->SetLabelSize(0.06);
  hRatio->Draw("E1");

  // Draw a line at the median:
  TLine *line1 = new TLine();
  line1->SetLineStyle(2);
  line1->SetLineWidth(2);
  line1->SetLineColor(kRed+1);
  line1->DrawLine(hRatio->GetXaxis()->GetXmin(), 1.0, 
		  hRatio->GetXaxis()->GetXmax(), 1.0);
  
  // Print the canvas:
  TString printName = Form("%s/comparison_%s.eps", m_outputDir.Data(),
			   name.Data());
  can->Print(printName);

  //delete can;
  //delete pad1;
  //delete pad2;
  //delete hRatio;
  //delete line1;
}

/**
   -----------------------------------------------------------------------------
   The main method tosses toys and saves data in a TTree.
   @param configFile - The name of the analysis config file.
   @param options - The options (see header note).
*/
int main(int argc, char **argv) {
  if (argc < 5) {
    std::cout << "Usage: " << argv[0]
	      << " <configFile> <file1> <file2> <options>"
	      << std::endl;
    exit(0);
  }
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString fileName1 = argv[2];
  TString fileName2 = argv[3];
  TString options = argv[4];
  
  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  
  // Construct the output directory:
  m_outputDir = Form("%s/%s/CompareData",
		     (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/", m_outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Open files, access histograms:
  TFile file1(fileName1);
  TFile file2(fileName2);
  
  // Loop over contents of file to find histogram names:
  TIter next(file1.GetListOfKeys());
  TObject *currObj;
  while ((currObj = (TObject*)next())) {
    TString currName = currObj->GetName();
    //std::cout << "current object name = " << currName << std::endl;
    if (currName.Contains("hist_")) {
      TH1F *h1 = (TH1F*)file1.Get(currName);
      TH1F *h2 = (TH1F*)file2.Get(currName);
      compareHistograms(currName, h1, h2);
      delete h1;
      delete h2;
    }
  }
  
  // Close files before exiting:
  file1.Close();
  file2.Close();
  
  delete currObj;
  delete config;
  
  return 0;
}
