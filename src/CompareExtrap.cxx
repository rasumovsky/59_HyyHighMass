////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CompareExtrap.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/04/2016                                                          //
//                                                                            //
//  This program plots the results of ExtrapolateSig.cxx for the scalar and   //
//  graviton searches with and without signal in 2016.                        //
//                                                                            //
//  Options:                                                                  //
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
  if (argc < 5) {
    std::cout << "Usage: " << argv[0]
	      << " <configFile> <jobNameGraviton> <jobNameScalar> <options>"
	      << std::endl;
    exit(0);
  }
  
  // Assign input parameters:
  TString configFile = argv[1];
  TString jobNameGraviton = argv[2];
  TString jobNameScalar = argv[3];
  TString options = argv[4];

  // Load the analysis configurations from file:
  Config *config = new Config(configFile);
  
  // Construct the output directory:
  TString outputDir = Form("%s/%s/CompareExtrap", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Open files, access histograms:
  TString fileTag = options.Contains("Only2016") ? "Only2016" : "Combined";
  TFile fGravitonMu0(Form("%s/%s/ExtrapolateSig/sigExtrap_Mu0_Graviton_%s.root",
			  (config->getStr("MasterOutput")).Data(),
			  jobNameGraviton.Data(), fileTag.Data()));
  TFile fGravitonMu1(Form("%s/%s/ExtrapolateSig/sigExtrap_Mu1_Graviton_%s.root",
			  (config->getStr("MasterOutput")).Data(),
			  jobNameGraviton.Data(), fileTag.Data()));
  TFile fScalarMu0(Form("%s/%s/ExtrapolateSig/sigExtrap_Mu0_Scalar_%s.root",
			(config->getStr("MasterOutput")).Data(),
			jobNameScalar.Data(), fileTag.Data()));
  TFile fScalarMu1(Form("%s/%s/ExtrapolateSig/sigExtrap_Mu1_Scalar_%s.root",
			(config->getStr("MasterOutput")).Data(),
			jobNameScalar.Data(), fileTag.Data()));
  
  TGraphAsymmErrors *gErr_Graviton_Mu1 = options.Contains("SqrtL") ? 
    (TGraphAsymmErrors*)fGravitonMu1.Get("gZOvsLumiSqrt_err_Mu1") :
    (TGraphAsymmErrors*)fGravitonMu1.Get("gZOvsLumi_err_Mu1");
  gErr_Graviton_Mu1->SetFillStyle(3254);
  gErr_Graviton_Mu1->SetFillColor(kRed+1);
  gErr_Graviton_Mu1->SetLineColor(kRed+1);
  
  TGraphAsymmErrors *gErr_Scalar_Mu1 = options.Contains("SqrtL") ? 
    (TGraphAsymmErrors*)fScalarMu1.Get("gZOvsLumiSqrt_err_Mu1") :
    (TGraphAsymmErrors*)fScalarMu1.Get("gZOvsLumi_err_Mu1");
  gErr_Scalar_Mu1->SetFillStyle(3245);
  gErr_Scalar_Mu1->SetFillColor(kBlue+1);
  gErr_Scalar_Mu1->SetLineColor(kBlue+1);
    
  TGraph *gNom_Graviton_Mu1 = options.Contains("SqrtL") ? 
    (TGraph*)fGravitonMu1.Get("gZOvsLumiSqrt_nom_Mu1") :
    (TGraph*)fGravitonMu1.Get("gZOvsLumi_nom_Mu1");
  gNom_Graviton_Mu1->SetLineColor(kRed+1);
  
  TGraph *gNom_Scalar_Mu1 = options.Contains("SqrtL") ? 
    (TGraph*)fScalarMu1.Get("gZOvsLumiSqrt_nom_Mu1") :
    (TGraph*)fScalarMu1.Get("gZOvsLumi_nom_Mu1");
  gNom_Scalar_Mu1->SetLineColor(kBlue+1);
  
  TGraph *gNom_Graviton_Mu0 = options.Contains("SqrtL") ? 
    (TGraph*)fGravitonMu0.Get("gZOvsLumiSqrt_nom_Mu0") :
    (TGraph*)fGravitonMu0.Get("gZOvsLumi_nom_Mu0");
  gNom_Graviton_Mu0->SetLineColor(kRed+1);
  gNom_Graviton_Mu0->SetLineStyle(2);
  
  TGraph *gNom_Scalar_Mu0 = options.Contains("SqrtL") ? 
    (TGraph*)fScalarMu0.Get("gZOvsLumiSqrt_nom_Mu0") :
    (TGraph*)fScalarMu0.Get("gZOvsLumi_nom_Mu0");
  gNom_Scalar_Mu0->SetLineColor(kBlue+1);
  gNom_Scalar_Mu0->SetLineStyle(2);
  
  double xMin = 0.0; double xMax = 8.0;
  //double yMin = 2.5; double yMax = 8.5;
  double yMin = 0.0; double yMax = 8.5;
  
  TH1F *hForRange = new TH1F("hRange", "hRange", 1, xMin, xMax);
  hForRange->Fill(1);
  hForRange->SetLineColor(0);
  hForRange->GetYaxis()->SetRangeUser(yMin, yMax);
  hForRange->GetXaxis()->SetRangeUser(xMin, xMax);
  hForRange->GetXaxis()->SetTitle(gErr_Scalar_Mu1->GetXaxis()->GetTitle());
  hForRange->GetYaxis()->SetTitle(gErr_Scalar_Mu1->GetYaxis()->GetTitle());

  // Create the TCanvas and begin plotting:
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  hForRange->Draw("axis");
  gErr_Graviton_Mu1->Draw("3SAME");
  gErr_Scalar_Mu1->Draw("3SAME");
  gNom_Graviton_Mu1->Draw("LSAME");
  gNom_Scalar_Mu1->Draw("LSAME");
  if (!options.Contains("Only2016")) {
    gNom_Graviton_Mu0->Draw("LSAME");
    gNom_Scalar_Mu0->Draw("LSAME");
  }
  
  // 5 sigma line:
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(1);
  line->DrawLine(xMin, 5.0, xMax, 5.0);
  
  // Create a TLegend:
  TLegend leg(0.20, 0.75, 0.60, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);

  leg.AddEntry(gErr_Graviton_Mu1,
	       "Profiled #sigma_{G*}#timesBR(G*#rightarrow#gamma#gamma)", "F");
  leg.AddEntry(gErr_Scalar_Mu1,
	       "Profiled #sigma_{X}#timesBR(X#rightarrow#gamma#gamma)", "F");
  if (!options.Contains("Only2016")) {
    leg.AddEntry(gNom_Graviton_Mu0, "Absence of signal", "L");
    leg.AddEntry(gNom_Scalar_Mu0, "Absence of signal", "L");
  }
  leg.Draw("SAME");
  
  // Then print canvas:
  if (options.Contains("SqrtL")) {
    can->Print(Form("%s/plot_Z0vsSqrtLumi_comparison_%s.eps",
		    outputDir.Data(), fileTag.Data()));
  }
  else {
    can->Print(Form("%s/plot_Z0vsLumi_comparison_%s.eps",
		    outputDir.Data(), fileTag.Data()));
  }
  
  return 0;
}
