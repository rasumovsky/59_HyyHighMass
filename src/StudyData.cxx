////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  StudyData.cxx                                                             //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 29/04/2016                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for plotting various quantities from     //
//  data MxAODs.                                                              //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "HGammaMxAOD.h"

// Globally scoped variables:
TString m_outputDir;
std::map<TString,TH1F*> m_histograms;
std::map<TString,TH2D*> m_histograms2D;
HGammaMxAOD *m_treeMxAOD;
TString m_categoryName;

int m_cutFlowCounter_Hist[100];
int m_cutFlowCounter_Flag[100];
std::vector<TString> m_cutNames;

std::vector<int> m_runList;
std::map<int,int> m_eventsPerRun;

TString m_options;

/**
   -----------------------------------------------------------------------------
   Retrieve the total number of events in the file for event normalization:
   @param file - The current MxAOD file in the TChain.
   @param maxCutFlowIndex - An int to store the number of cutflow hist bins.
   @return - The total number of weighted events in the file. Also, the max
   cutFlowIndex passed by reference.
*/
double getNTotEvtFromHist(TFile *file, int& maxCutFlowIndex) {
  // Find the cutflow histograms from the file based on limited name info:
  TIter next(file->GetListOfKeys());
  TObject *currObj;
  while ((currObj = (TObject*)next())) {
    TString currName = currObj->GetName();
    if (currName.Contains("CutFlow") && currName.Contains("weighted")
	&& currName.Contains("noDalitz")) {
      maxCutFlowIndex = ((TH1F*)file->Get(currName))->GetNbinsX();
      return (((TH1F*)file->Get(currName))->GetBinContent(3));
      //return (((TH1F*)file->Get(currName))->GetBinContent(3) * 
      //      ((TH1F*)file->Get(currName))->GetBinContent(2) /
      //      ((TH1F*)file->Get(currName))->GetBinContent(1));
    }
  }
  std::cout << "createSignalParameterization: ERROR! MxAOD doesn't have cutflow"
	    << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Get names and values of cuts from MxAOD directly. Also look at DMxAOD...
   @param file - The current MxAOD file in the TChain.
   @param maximumBin - The last bin for which cutflow data will be filled.
*/
void fillCutFlowFromMxAOD(TFile *file, int maximumBin) {
  // Find the cutflow histograms from the file based on limited name info:
  TH1F *cutFlowHist = NULL;
  TIter next(file->GetListOfKeys());
  TObject *currObj;
  while ((currObj = (TObject*)next())) {
    TString currName = currObj->GetName();
    if (currName.Contains("CutFlow_Run") || 
	(currName.Contains("CutFlow") && currName.Contains("noDalitz") && 
	 !currName.Contains("weighted"))) {
      cutFlowHist = (TH1F*)file->Get(currName);
      break;
    }
  }
  
  bool cutsDefined = ((int)m_cutNames.size() > 0);
  
  // Loop over bins to fill cutflow table and also get names of cuts (1st pass):
  for (int i_b = 1; i_b <= maximumBin; i_b++) {
    m_cutFlowCounter_Hist[i_b-1] += (int)(cutFlowHist->GetBinContent(i_b));
    // Get names of cuts:
    if (!cutsDefined) {
      m_cutNames.push_back(cutFlowHist->GetXaxis()->GetBinLabel(i_b));
    }
  }
  if (!cutsDefined) {
    m_cutNames.push_back("Preselection");
    m_cutNames.push_back("PID");
    m_cutNames.push_back("m_yy");
    m_cutNames.push_back("isolation");
    m_cutNames.push_back("pT_cuts");
  }
}

/**
   -----------------------------------------------------------------------------
   Select the proper category for the current event in the HGammaMxAOD. The 
   categorization to use is based on the global variable m_categoryName.
   @return - The category index.
*/
int chooseCategory() {

  // Eta categorization:
  if (m_categoryName.EqualTo("EtaCate")) {
    double barrelEnd = 1.37;
    double eta1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0];
    double eta2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1];
    if (fabs(eta1) < barrelEnd && fabs(eta2) < barrelEnd) return 0;
    else if (fabs(eta1) < barrelEnd) return 1;
    else if (fabs(eta2) < barrelEnd) return 2;
    else return 3;
  }
  
  // Eta categorization:
  if (m_categoryName.EqualTo("OSEndcap")) {
    double barrelEnd = 1.37;
    double eta1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0];
    double eta2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1];
    bool oppositeSign = ((eta1 >= 0 && eta2 < 0) || (eta1 < 0 && eta2 >= 0));
    if (fabs(eta1)>barrelEnd && fabs(eta2)>barrelEnd && oppositeSign) return 0;
    else return 1;
  }
  
  // Conversion categorization:
  else if (m_categoryName.EqualTo("ConversionCate")) {
    int conv1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[0];
    int conv2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[1];
    if (conv1 == 0 && conv2 == 0) return 0;
    else if (conv1 == 0) return 1;
    else if (conv2 == 0) return 2;
    else return 3;
  }
 
  // Mass categorization:
  else if (m_categoryName.EqualTo("MassCate")) {
    if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy < 700000) return 0;
    else if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy >= 700000 &&
	     m_treeMxAOD->HGamEventInfoAuxDyn_m_yy < 800000) return 1;
    else if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy >= 800000) return 2;
    else return 3;
  }
  
  // Gain categorization:
  else if (m_categoryName.EqualTo("GainCate")) {
    int gain1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_gain)[0];
    int gain2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_gain)[1];
    if (gain1 == 0 && gain2 == 0) return 0;
    else if (gain1 == 0) return 1; 
    else if (gain2 == 0) return 2;
    else return 3;
  }
  
  // Exit because inputs are unknown:
  else {
    std::cout << "StudyData: ERROR! Categorization not found" << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the name for the current category. The categorization to use is based on 
   the global variable m_categoryName.
   @param categoryIndex - The index of the current category.
   @return - The category name.
*/
TString nameCategory(int categoryIndex) {
  if (m_categoryName.EqualTo("EtaCate")) {
    if (categoryIndex == 0) return "BarrelBarrel";
    else if (categoryIndex == 1) return "BarrelEndcap";
    else if (categoryIndex == 2) return "EndcapBarrel";
    else return "EndcapEndcap";
  }
  else if (m_categoryName.EqualTo("OSEndcap")) {
    if (categoryIndex == 0) return "OppositeSideEndcap";
    else return "Rest";
  }
  else if (m_categoryName.EqualTo("ConversionCate")) {
    if (categoryIndex == 0) return "UnconvUnconv";
    else if (categoryIndex == 1) return "UnconvConv";
    else if (categoryIndex == 2) return "ConvUnconv";
    else return "ConvConv";
  }
  else if (m_categoryName.EqualTo("MassCate")) {
    if (categoryIndex == 0) return "mgg_600_700";
    else if (categoryIndex == 1) return "mgg_700_800";
    else if (categoryIndex == 2) return "mgg_800_inf";
    else return "mgg_rest";
  }
  else if (m_categoryName.EqualTo("GainCate")) {
    if (categoryIndex == 0) return "gain_0_0";
    else if (categoryIndex == 1) return "gain_0_1";
    else if (categoryIndex == 2) return "gain_1_0";
    else return "gain_1_1";
  }
  
  // Exit because inputs are unknown:
  else {
    std::cout << "StudyData: ERROR! Categorization not found" << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the name for the current category. The categorization to use is based on 
   the global variable m_categoryName.
   @param categoryIndex - The index of the current category.
   @return - The category name.
*/
TString histToCateName(TString histName) {
  if (histName.Contains("OppositeSideEndcap")) return "Opposite-side endcap";
  else if (histName.Contains("Rest")) return "Rest";
  else if (histName.Contains("BarrelBarrel")) return "2 Barrel #gamma's";
  else if (histName.Contains("BarrelEndcap")) return "Barrel-Endcap";
  else if (histName.Contains("EndcapBarrel")) return "Endcap-Barrel";
  else if (histName.Contains("EndcapEndcap")) return "2 Endcap #gamma's";
  else if (histName.Contains("UnconvUnconv")) return "2 Unconverted #gamma's";
  else if (histName.Contains("UnconvConv")) return "Unconverted-Converted";
  else if (histName.Contains("ConvUnconv")) return "Converted-Unconverted";
  else if (histName.Contains("ConvConv")) return "2 Converted #gamma's";
  else if (histName.Contains("mgg_600_700")) {
    return "m_{#gamma#gamma} < 700 GeV";
  }
  else if (histName.Contains("mgg_700_800")) {
    return "700 GeV < m_{#gamma#gamma} < 800 GeV";
  }
  else if (histName.Contains("mgg_800_inf")) {
    return "m_{#gamma#gamma} > 800 GeV";
  }
  else if (histName.Contains("mgg_rest")) return "mgg_rest";
  else if (histName.Contains("gain_0_0")) return "gain(#gamma1,#gamma2)=0";
  else if (histName.Contains("gain_0_1")) return "gain(#gamma1)=0";
  else if (histName.Contains("gain_1_0")) return "gain(#gamma2)=0";
  else if (histName.Contains("gain_1_1")) return "gain(#gamma1,#gamma2)>0";
  else return "";
}

/**
   -----------------------------------------------------------------------------
   Instantiate an inclusive and categorized histogram.
   @param histName - The name of the histogram.
   @param nCategories - The number of analysis categories for plotting.
   @param nBins - The number of histogram bins.
   @param xMin - The minimum of the histogram x-axis.
   @param xMax - The maximum of the histogram x-axis.
*/
void defineHistograms(TString histName, int nCategories, int nBins, double xMin,
		      double xMax) {
  double binning = (xMax - xMin) / ((double)nBins);
  
  // Create inclusive histogram:
  m_histograms[histName] = new TH1F(histName, histName, nBins, xMin, xMax);
  m_histograms[histName]->GetXaxis()->SetTitle(histName);
  m_histograms[histName]->Sumw2();
  if (histName.Contains("GeV")) {
    m_histograms[histName]
      ->GetYaxis()->SetTitle(Form("Entries / %2.1f GeV",binning));
  }
  else {
    m_histograms[histName]
      ->GetYaxis()->SetTitle(Form("Entries / %2.1f",binning));
  }

  // Create categorized histograms:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    TString cateName = nameCategory(i_c);
    TString currName = histName + "_" + cateName;
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      = new TH1F(currName, currName, nBins, xMin, xMax);
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetXaxis()->SetTitle(histName);
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]->Sumw2();
    if (histName.Contains("GeV")) {
      m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
	->GetYaxis()->SetTitle(Form("Entries / %2.1f GeV",binning));
    }
    else {
      m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
	->GetYaxis()->SetTitle(Form("Entries / %2.1f",binning));
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Instantiate an inclusive and categorized histogram.
   @param histName - The name of the histogram.
   @param nCategories - The number of analysis categories for plotting.
   @param xName - The name of the x-axis variable.
   @param xBins - The number of x-axis histogram bins.
   @param xMin - The minimum of the histogram x-axis.
   @param xMax - The maximum of the histogram x-axis.
   @param yName - The name of the y-axis variable.
   @param yBins - The number of y-axis histogram bins.
   @param yMin - The minimum of the histogram y-axis.
   @param yMax - The maximum of the histogram y-axis.
*/
void define2DHistograms(TString histName, int nCategories, TString xName, 
			int xBins, double xMin, double xMax, TString yName,
			int yBins, double yMin, double yMax) {
  
  double binningX = (xMax - xMin) / ((double)xBins);
  double binningY = (yMax - yMin) / ((double)yBins);
  
  // Create inclusive histogram:
  m_histograms2D[histName]
    = new TH2D(histName, histName, xBins, xMin, xMax, yBins, yMin, yMax);
  m_histograms2D[histName]->GetXaxis()->SetTitle(xName);
  m_histograms2D[histName]->GetYaxis()->SetTitle(yName);
  m_histograms2D[histName]->Sumw2();
  m_histograms2D[histName]->GetZaxis()->SetTitle("Entries");
  
  // Create categorized histograms:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    TString cateName = nameCategory(i_c);
    TString currName = histName + "_" + cateName;
    m_histograms2D[Form("%s_%s",histName.Data(), cateName.Data())]
      = new TH2D(currName, currName, xBins, xMin, xMax, yBins, yMin, yMax);
    m_histograms2D[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetXaxis()->SetTitle(xName);
    m_histograms2D[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetYaxis()->SetTitle(yName);
    m_histograms2D[Form("%s_%s",histName.Data(), cateName.Data())]->Sumw2();
    m_histograms2D[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetZaxis()->SetTitle("Entries");
  }
}

/**
   -----------------------------------------------------------------------------
   Fill inclusive and categorized histograms.
   @param histName - The name of the histogram.
   @param category -The index of the category into which this event falls.
   @param value - The value to fill into the histogram for this event.
   @param weight - The event weight.
*/
void fillHistograms(TString histName, int category, double value,
		    double weight) {
    
  m_histograms[histName]->Fill(value, weight);
  TString cateName = nameCategory(category);
  m_histograms[Form("%s_%s", histName.Data(), cateName.Data())]
    ->Fill(value, weight);
}

/**
   -----------------------------------------------------------------------------
   Fill inclusive and categorized histograms.
   @param histName - The name of the histogram.
   @param category -The index of the category into which this event falls.
   @param xValue - The x value to fill into the histogram for this event.
   @param yValue - The y value to fill into the histogram for this event.
   @param weight - The event weight.
*/
void fill2DHistograms(TString histName, int category, double xValue, 
		      double yValue, double weight) {
    
  m_histograms2D[histName]->Fill(xValue, yValue, weight);
  TString cateName = nameCategory(category);
  m_histograms2D[Form("%s_%s", histName.Data(), cateName.Data())]
    ->Fill(xValue, yValue, weight);
}

/**
   -----------------------------------------------------------------------------
   Format the histogram name so that it can be used for file names.
   @param histName - The name of the histogram.
   @return - A formatted histogram name.
*/
TString formatHistName(TString histName) {
  histName.ReplaceAll("/","");
  histName.ReplaceAll("(",""); 
  histName.ReplaceAll(")","");
  histName.ReplaceAll("[","");
  histName.ReplaceAll("]","");
  histName.ReplaceAll("{","");
  histName.ReplaceAll("}","");
  histName.ReplaceAll(">","");
  histName.ReplaceAll("<","");
  histName.ReplaceAll("#","");
  histName.ReplaceAll("^","");
  histName.ReplaceAll("*","");
  histName.ReplaceAll(" ","");
  return histName;
}

/**
   -----------------------------------------------------------------------------
   Plot the named histogram.
   @param histName - The name of the histogram.
   @param ana - The analysis type.
   @param label - The ATLAS label
   @param lumi - The dataset luminosity.
*/
void plotHistogram(TString histName, TString ana, TString label, double lumi) {
 
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  gPad->SetLogy();
  m_histograms[histName]->Draw("E1");
  histName = formatHistName(histName);
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.6, 0.87, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.72, 0.87, label);
  t.DrawLatex(0.6, 0.81, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",lumi));
  if (ana.Contains("Scalar")) t.DrawLatex(0.6, 0.75, "Spin-0 Selection");
  else if (ana.Contains("GravitonLoose")) {
    t.DrawLatex(0.6, 0.75, "Spin-2 Loose Iso.");
  }
  else t.DrawLatex(0.6, 0.75, "Spin-2 Selection");
  t.DrawLatex(0.6, 0.69, histToCateName(histName));

  TString printName = Form("%s/plot_%s.eps",m_outputDir.Data(),histName.Data());
  can->Print(printName);
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Plot the named histogram.
   @param histName - The name of the histogram.
   @param ana - The analysis type.
   @param label - The ATLAS label
   @param lumi - The dataset luminosity.
*/
void plot2DHistogram(TString histName, TString ana, TString label, double lumi){
  std::cout << "plot2DHistogram(" << histName << ", " << ana << ", " << label
	    << ", " << lumi << ")" << std::endl;
  
  TCanvas *can = new TCanvas("can", "can", 800, 1000);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.51, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.49);
  pad1->SetBottomMargin(0.00001);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.2);
  pad1->SetBorderMode(0);
  pad2->SetBorderMode(0);
  pad1->SetRightMargin(0.15);
  pad2->SetRightMargin(0.15);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  
  pad1->SetLogz();
  m_histograms2D[histName]->GetZaxis()
    ->SetRangeUser(0.1, m_histograms2D[histName]->GetMaximum());
  m_histograms2D[histName]->SetContour(256);
  m_histograms2D[histName]->Draw("colz");
  //m_histograms2D[histName]->Draw("scat=1.0");
  
  histName = formatHistName(histName);
  /*
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(kRed+1);
  line->DrawLine(700, m_histograms2D[histName]->GetYaxis()->GetXmin(),
  		 700, m_histograms2D[histName]->GetYaxis()->GetXmax());
  line->DrawLine(800, m_histograms2D[histName]->GetYaxis()->GetXmin(),
  		 800, m_histograms2D[histName]->GetYaxis()->GetXmax());
  */
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.55, 0.86, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.67, 0.86, label);
  t.DrawLatex(0.55, 0.80, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",lumi));
  if (ana.Contains("Scalar")) t.DrawLatex(0.55, 0.74, "Spin-0 Selection");
  else if (ana.Contains("GravitonLoose")) {
    t.DrawLatex(0.55, 0.74, "Spin-2 Loose Iso.");
  }
  else t.DrawLatex(0.55, 0.74, "Spin-2 Selection");
  t.DrawLatex(0.55, 0.68, histToCateName(histName));
  
  pad2->cd();
  
  TString printName
    = Form("%s/plot2D_%s.eps", m_outputDir.Data(), histName.Data());
  
  // Load the MC profile:
  if (m_options.Contains("CompareMC")) {
        
    TFile file("MCFile.root", "READ");
    TH2D *temp = (TH2D*)file.Get(Form("hist_%s",
				      (formatHistName(histName)).Data()));
    TProfile *profMC = temp->ProfileX();
    profMC->SetNameTitle("profMC", "profMC");
    profMC->GetYaxis()
      ->SetTitle(m_histograms2D[histName]->GetYaxis()->GetTitle());
    profMC->SetLineColor(kRed);
    profMC->SetMarkerColor(kRed);
    
    TProfile *profile = m_histograms2D[histName]->ProfileX();
    profile->GetYaxis()
      ->SetTitle(m_histograms2D[histName]->GetYaxis()->GetTitle());
    profile->Draw("E1");
    profMC->Draw("E1SAME");
    
    /*    
    line->DrawLine(700, profMC->GetMinimum(),
		   700, profMC->GetMaximum());
    line->DrawLine(800, profMC->GetMinimum(),
		   800, profMC->GetMaximum());
    */
    
    TLegend leg(0.62, 0.23, 0.81, 0.33);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetTextSize(0.04);
    leg.AddEntry(profile, "Data","EP");
    leg.AddEntry(profMC, "Sherpa LO MC","EP");
    leg.Draw("SAME");
    
    can->Print(printName);
    file.Close();
  }
  else {
    TProfile *profile = m_histograms2D[histName]->ProfileX();
    profile->GetYaxis()
      ->SetTitle(m_histograms2D[histName]->GetYaxis()->GetTitle());
    profile->Draw("E1");
    /*
    line->DrawLine(700, profile->GetMinimum(),
		   700, profile->GetMaximum());
    line->DrawLine(800, profile->GetMinimum(),
		   800, profile->GetMaximum());
    */
    can->Print(printName);
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Plot the named histogram compared with hists in other categories.
   @param histName - The name of the histogram.
   @param nCategories - The number of analysis categories for plotting.
   @param ana - The analysis type.
   @param label - The ATLAS label
   @param lumi - The dataset luminosity.
*/
void plotComparisonHist(TString histName, int nCategories, TString ana,
			TString label, double lumi) {
  std::cout << "plotComparisonHist(" << histName << ", " << nCategories
	    << ", " << ana << ", " << label << ", " << lumi << ")" << std::endl;
  
  // Get hist name without categorization information:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    TString currCateName = nameCategory(i_c);
    if (histName.Contains(Form("_%s",currCateName.Data()))) {
      return;
    }
  }
  
  // Create the canvas:
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  
  // Create a legend:
  TLegend leg(0.6, 0.72, 0.9, 0.91);
  leg.SetTextFont(42); 
  leg.SetTextSize(0.03);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  
  // Plot each histogram and add legend entry:
  Color_t lineColors[5] = {kBlue+1, kMagenta+1, kGreen+1, kCyan+1, kRed+1};
  for (int i_c = nCategories-1; i_c >= 0; i_c--) {
    if (i_c > 2) continue;
    TString currCateName = nameCategory(i_c);
    TString currHistName = Form("%s_%s", histName.Data(), 
				currCateName.Data());
    m_histograms[currHistName]
      ->Scale(1.0 / m_histograms[currHistName]->Integral());
    m_histograms[currHistName]->SetLineColor(lineColors[i_c]);
    m_histograms[currHistName]->SetMarkerColor(lineColors[i_c]);
    m_histograms[currHistName]->SetMarkerStyle(21 + i_c);
    if (i_c == nCategories-1) {
      m_histograms[currHistName]->GetYaxis()
	->SetRangeUser(0, 2.0 * m_histograms[currHistName]->GetMaximum());
      m_histograms[currHistName]->Draw("E1");
    }
    else m_histograms[currHistName]->Draw("E1SAME");
    leg.AddEntry(m_histograms[currHistName],
		 histToCateName(nameCategory(i_c)), "LEP");
  }
  histName = formatHistName(histName);
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.04);
  t.DrawLatex(0.35, 0.87, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.04);
  t.DrawLatex(0.45, 0.87, label);
  t.DrawLatex(0.35, 0.82, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",lumi));
  if (ana.Contains("Scalar")) t.DrawLatex(0.35, 0.77, "Spin-0 Selection");
  else if (ana.Contains("GravitonLoose")) {
    t.DrawLatex(0.35, 0.77, "Spin-2 Loose Iso.");
  }
  else t.DrawLatex(0.6, 0.77, "Spin-2 Selection");
  t.DrawLatex(0.35, 0.72, histToCateName(histName));
  leg.Draw("SAME");

  TString printName = Form("%s/plot_comparison_%s.eps", m_outputDir.Data(), 
			   histName.Data());
  can->Print(printName);
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @return - An updated list of file names.
*/
std::vector<TString> makeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "StudyData: Making local copies of inputs."
	    << std::endl;
  std::vector<TString> result; result.clear();
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    TString newName = Form("tempFile%d.root", i_f);
    if (fileNames[i_f].Contains("root://eosatlas/")) {
      system(Form("xrdcp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else if (fileNames[i_f].Contains("/eos/atlas/")) {
      system(Form("eos cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else {
      system(Form("cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    result.push_back(newName);
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Remove any files that were copied over for speed.
   @param fileNames - The original file names.
*/
void removeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "StudyData: Removing local copies of inputs."
	    << std::endl;
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    system(Form("rm %s", fileNames[i_f].Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Returns the maximum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @return - The maximum entry in the vector.
*/
double maxEntry(std::vector<double> currList) {
  double maximum = 0.0;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter > maximum) maximum = *mIter;
  }
  return maximum;
}

/**
   -----------------------------------------------------------------------------
   Returns the minimum entry in a vector of doubles.
   @param currList - The vector of doubles.
   @return - The minimum entry in the vector.
*/
double minEntry(std::vector<double> currList) {
  double minimum = 100000.0;
  for (std::vector<double>::iterator mIter = currList.begin(); 
       mIter != currList.end(); mIter++) {
    if (*mIter < minimum) minimum = *mIter;
  }
  return minimum;
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
   Check if a vector contains a value, and add the value if it doesn't.
   @param value - The value to add to the vector.
*/
void vectorAdd(int value) {
  for (int i_v = 0; i_v < (int)m_runList.size(); i_v++) {
    if (m_runList[i_v] == value) return;
  }
  m_runList.push_back(value);
  m_eventsPerRun[value] = 0;
}

/**
   -----------------------------------------------------------------------------
   The main method for this utility. Provide 1 argument - the location of the 
   config (.cfg) file, which should be stored in the data/ directory. The main()
   method runs over the samples provided, performs the fits requests, and gives
   comparisons of parameterized and non-parameterized fits. 
*/
int main(int argc, char *argv[])
{
  // Check that arguments are provided.
  if (argc < 2) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <options>"
	      << std::endl;
    exit(0);
  }
  
  Config *config = new Config(TString(argv[1]));
  m_options = argv[2];

  // Check that output directory exists:
  m_outputDir = Form("%s/%s/StudyData", (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Define output file with variables:
  std::ofstream outputForGAN(Form("%s/variablesForGAN.txt",m_outputDir.Data()));
  
  // Define per-event histograms:
  m_histograms.clear();
  int nBins = 10;
  int nCategories = config->getInt("MxAODNCategories");
  m_categoryName = config->getStr("MxAODCategorization");
  
  defineHistograms("z_{vertex} [mm]", nCategories, nBins, -150, 150);
  if ((config->getStr("AnalysisType")).Contains("Scalar")) {
    if (config->getBool("DoBlind")) {
      defineHistograms("m_{#gamma#gamma} [GeV]", nCategories, 44, 150, 700);
    }
    else {
      defineHistograms("m_{#gamma#gamma} [GeV]", nCategories, 74, 150, 2000);
    }
  }
  else {
    if (config->getBool("DoBlind")) {
      defineHistograms("m_{#gamma#gamma} [GeV]", nCategories, 40, 200, 700);
    }
    else {
      defineHistograms("m_{#gamma#gamma} [GeV]", nCategories, 72, 200, 2000);
    }
  }
  defineHistograms("p_{T}^{#gamma#gamma} [GeV]", nCategories, nBins, 0, 1000);
  defineHistograms("cos(#theta*)", nCategories, nBins, 0, 1);
  defineHistograms("N_{Jets}", nCategories, 9, -0.5, 8.5);
  defineHistograms("N_{Photons}", nCategories, 3, 1.5, 4.5);
  defineHistograms("N_{Electrons}", nCategories, 5, -0.5, 4.5);
  defineHistograms("N_{Muons}", nCategories, 5, -0.5, 4.5);

  // Define histograms of per-photon variables:
  for (int i_p = 0; i_p < 2; i_p++) {
    defineHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0.0, 1000.0);
    defineHistograms(Form("#eta(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -2.5, 2.5);
    defineHistograms(Form("#eta_{S2}(#gamma_{%d})",i_p+1),
		     nCategories, 11, -3.0, 3.0);
    defineHistograms(Form("#phi(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -3.141, 3.141);
    defineHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1), 
		     nCategories, nBins, 0.0, 50.0);
    defineHistograms(Form("E_{T}^{cone40}(#gamma_{%d}) [GeV]",i_p+1), 
		     nCategories, nBins, -5.0, 25.0);
    defineHistograms(Form("Conversion type #gamma_{%d}",i_p+1),
		     nCategories, 6, -0.5, 5.5);
    defineHistograms(Form("gain_{max-E} #gamma_{%d}",i_p+1),
		     nCategories, 3, -0.5, 2.5);
    
    // Also define some 2D histograms:
    define2DHistograms(Form("timing_maxEcell%d",i_p+1),
		       nCategories, "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("timing_{max-E cell} #gamma_{%d} [ns]", i_p+1),
		       50, -2, 2);
    define2DHistograms(Form("myy_vs_PID%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("PID_{#gamma%d}",i_p+1), 2, -0.5, 1.5);    
    define2DHistograms(Form("myy_vs_pt%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1), 25, 0, 1000);
    define2DHistograms(Form("myy_vs_fabseta%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("#eta(#gamma_{%d})",i_p+1), 25, 0, 2.5);
    define2DHistograms(Form("myy_vs_phi%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("#phi(#gamma_{%d})",i_p+1), 25,
		       -1.0*TMath::Pi()/2.0, TMath::Pi()/2.0);
    define2DHistograms(Form("myy_vs_conv%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("conv. bit #gamma_{%d}",i_p+1), 6, -0.5, 5.5);
    define2DHistograms(Form("myy_vs_gain%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("gain #gamma_{%d}",i_p+1), 3, -0.5, 2.5);
    define2DHistograms(Form("myy_vs_trkiso%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("trk iso. #gamma_{%d}) [GeV]",i_p+1), 50, 0, 10);
    define2DHistograms(Form("myy_vs_caloiso%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("calo. iso. #gamma_{%d}) [GeV]",i_p+1), 50, -10,10);
    
    // Shower-shape variables:
    define2DHistograms(Form("myy_vs_weta1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("weta1_{#gamma%d}",i_p+1), 40, 0.4, 0.9);
    
    define2DHistograms(Form("myy_vs_weta2_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("weta2_{#gamma%d}",i_p+1), 40, 0.4, 0.9);
    
    define2DHistograms(Form("myy_vs_wtots1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("wstots1_{#gamma%d}",i_p+1), 40, 0.0, 5.0);
    
    define2DHistograms(Form("myy_vs_e277_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("e277_{#gamma%d}",i_p+1), 40, 0.0, 1500.0);
    
    define2DHistograms(Form("myy_vs_relEreso_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("relEreso_{#gamma%d}",i_p+1), 40, 0.0, 1.0);
    
    define2DHistograms(Form("myy_vs_Eratio%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("Eratio_{#gamma%d}",i_p+1), 40, 0.7, 1.1);
    
    define2DHistograms(Form("myy_vs_Reta%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("Reta_{#gamma%d}",i_p+1), 40, 0.8, 1.0);
    
    define2DHistograms(Form("myy_vs_f1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("f1_{#gamma%d}",i_p+1), 40, 0.0, 1.0);
    
    define2DHistograms(Form("myy_vs_Rhad_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("Rhad_{#gamma%d}",i_p+1), 40, -0.04, 0.04);
    
    define2DHistograms(Form("myy_vs_Rhad1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("Rhad1_{#gamma%d}",i_p+1), 40, -0.04, 0.04);
    
    define2DHistograms(Form("myy_vs_Rphi_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("Rphi_{#gamma%d}",i_p+1), 40, 0.7, 1.0);
    
    define2DHistograms(Form("myy_vs_DeltaE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("DeltaE_{#gamma%d}",i_p+1), 50, -5.0, 5.0);
    
    define2DHistograms(Form("myy_vs_fracs1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("fracs1_{#gamma%d}",i_p+1), 50, 0.0, 1.0);
    
    // Ratio of sampling layers:
    define2DHistograms(Form("myy_vs_rawcl_ratioEs1Es2_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("ratio E_{S1}/E_{S2} #gamma_{%d} [GeV]",i_p+1),
		       50, 0, 2);
    
    // Cluster deposits:
    define2DHistograms(Form("myy_vs_clES0_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E0_{cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 200);
    define2DHistograms(Form("myy_vs_clES1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E1_{cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 500);
    define2DHistograms(Form("myy_vs_clES2_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E2_{cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 1000);
    define2DHistograms(Form("myy_vs_clES3_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E3_{cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 100);
    
    // Raw cluster deposits:
    define2DHistograms(Form("myy_vs_rawclES0_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E0_{raw cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 200);
    define2DHistograms(Form("myy_vs_rawclES1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E1_{raw cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 500);
    define2DHistograms(Form("myy_vs_rawclES2_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E2_{raw cl} #gamma_{%d} [GeV]",i_p+1), 50, 0,1000);
    define2DHistograms(Form("myy_vs_rawclES3_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E3_{raw cl} #gamma_{%d} [GeV]",i_p+1), 50, 0, 100);
        
    // ratio of cluster E to raw cluster E:
    define2DHistograms(Form("myy_vs_ratioES0_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E0_{cl}/E0_{raw cl} #gamma_{%d}",i_p+1),
		       50, 0.9, 1.2);
    define2DHistograms(Form("myy_vs_ratioES1_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E1_{cl}/E1_{raw cl} #gamma_{%d}",i_p+1),
		       50, 0.9, 1.2);
    define2DHistograms(Form("myy_vs_ratioES2_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E2_{cl}/E2_{raw cl} #gamma_{%d}",i_p+1),
		       50, 0.9, 1.2);
    define2DHistograms(Form("myy_vs_ratioES3_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E3_{cl}/E3_{raw cl} #gamma_{%d}",i_p+1),
		       50, 0.9, 1.2);
    
    // Ratio of cluster E to photon E
    define2DHistograms(Form("myy_vs_clES0overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E0_{cl}/E #gamma_{%d}",i_p+1), 50, 0, 0.5);
    define2DHistograms(Form("myy_vs_clES1overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E1_{cl}/E #gamma_{%d}",i_p+1), 50, 0, 0.5);
    define2DHistograms(Form("myy_vs_clES2overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E2_{cl}/E #gamma_{%d}",i_p+1), 50, 0, 1);
    define2DHistograms(Form("myy_vs_clES3overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E3_{cl}/E #gamma_{%d}",i_p+1), 50, 0, 0.3);

    // Ratio of raw cluster E to photon E
    define2DHistograms(Form("myy_vs_rawES0overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E0_{raw cl}/E #gamma_{%d}",i_p+1), 50, 0, 1);
    define2DHistograms(Form("myy_vs_rawES1overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E1_{raw cl}/E #gamma_{%d}",i_p+1), 50, 0, 1);
    define2DHistograms(Form("myy_vs_rawES2overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E2_{raw cl}/E #gamma_{%d}",i_p+1), 50, 0, 1);
    define2DHistograms(Form("myy_vs_rawES3overE_%d",i_p+1), nCategories, 
		       "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		       Form("E3_{raw cl}/E #gamma_{%d}",i_p+1), 50, 0, 1);
  }
  define2DHistograms("myy_vs_deta", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#Delta#eta_{#gamma#gamma}", 25, 0, 6);
  define2DHistograms("myy_vs_phi", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#Delta#phi_{#gamma#gamma}", 25, 0, 2.0*TMath::Pi());
  define2DHistograms("myy_vs_zvtx", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "z_{pointing} [mm]", 30, -400, 400);
  define2DHistograms("myy_vs_npv", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "N_{PV}", 40, -0.5, 40);
  

  define2DHistograms("myy_vs_actualIPX", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "Actual Interactions Per Crossing", 40, 0, 40);
  define2DHistograms("myy_vs_avgIPX", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "Avg. Interactions Per Crossing", 40, 0, 40);
  define2DHistograms("myy_vs_beamPosX", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "x_{beam} [mm]", 50, -25, 25);
  define2DHistograms("myy_vs_beamPosY", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "y_{beam} [mm]", 50, -25, 25);
  define2DHistograms("myy_vs_beamPosZ", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "z_{beam} [mm]", 50, -25, 25);
  define2DHistograms("myy_vs_beamPosSigmaX", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#sigma(x_{beam}) [mm]", 50, 0, 50);
  define2DHistograms("myy_vs_beamPosSigmaY", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#sigma(y_{beam}) [mm]", 50, 0, 50);
  define2DHistograms("myy_vs_beamPosSigmaZ", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#sigma(z_{beam}) [mm]", 50, 0, 50);
  define2DHistograms("myy_vs_beamPosSigmaXY", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#sigma(xy_{beam}) [mm]", 50, 0, 50);
  define2DHistograms("myy_vs_bunchDistanceFromFront", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "Distance From Front", 50, 0, 200);
  define2DHistograms("myy_vs_bunchGapBeforeTrain", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "Bunch Gap Before Train", 50, 0, 200);
  define2DHistograms("myy_vs_zSelected", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "z_{selected} [mm]", 50, -400, 400);
  define2DHistograms("myy_vs_zHardest", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "z_{hardest} [mm]", 50, -400, 400);
  define2DHistograms("myy_vs_zDiff", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "z_{hardest}-z_{selected} [mm]", 50, -400, 400);
  define2DHistograms("myy_vs_mu", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#mu", 40, 0, 40);
  define2DHistograms("myy_vs_pT_hard", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "p_{T}^{Hard} [GeV]", 50, 0, 500);
  define2DHistograms("myy_vs_met_TST", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#slash{E}_{T}^{TST} [GeV]", 50, 0, 2000);
  define2DHistograms("myy_vs_sumet_TST", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#sum E_{T}^{TST} [GeV]", 50, 0, 2000);
  define2DHistograms("myy_vs_met_phi", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "#phi_{#slash{E}_{T}^{TST}}", 40, -1*TMath::Pi(),
		     TMath::Pi());
  define2DHistograms("myy_vs_costheta", nCategories, 
		     "m_{#gamma#gamma} [GeV]", 36, 200, 2000, 
		     "cos(#theta*)", 50, 0, 1);
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = config->getStrV("MxAODsForData");
  // Make local copies of files if requested, to improve speed:
  if (config->getBool("MakeLocalMxAODCopies")) {
    fileNames = makeLocalFileCopies(fileNames);
  }
  // Create TChain of input files:
  TChain *chain = new TChain(config->getStr("MxAODTreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  m_treeMxAOD = new HGammaMxAOD(chain);
    
  //chain->MakeClass("test");
  
  // Count events:
  for (int i_c = 0; i_c < 50; i_c++) {
    m_cutFlowCounter_Hist[i_c] = 0;
    m_cutFlowCounter_Flag[i_c] = 0;
  }
  m_cutNames.clear();
  m_runList.clear();
  int countPass = 0;
  int countCate[50] = {0};
  m_eventsPerRun.clear();
  int prevRun = 0;
  
  // Variables for weighted events:
  double luminosity = config->getNum("AnalysisLuminosity");
  double nTotEvt = 1000.0;
  TString currFileName = "";
  int cutFlowIndex = 0;
  
  // Also create a text file of all passing events:
  std::ofstream eventList(Form("%s/eventList.txt", m_outputDir.Data()));
  
  //--------------------------------------//
  // Loop over events:
  int nEvents = m_treeMxAOD->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    m_treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    // Add each new file to the cutflow:
    if (!currFileName.EqualTo(chain->GetFile()->GetName())) {
      currFileName = chain->GetFile()->GetName();
      if (config->getBool("MxAODIsMC")) {
	nTotEvt = getNTotEvtFromHist(chain->GetFile(), cutFlowIndex);
      }
      fillCutFlowFromMxAOD(chain->GetFile(),
			   config->getInt("MxAODCutFlowIndex"));
    }
    
    if (prevRun != (int)(m_treeMxAOD->EventInfoAux_runNumber)) {
      prevRun = m_treeMxAOD->EventInfoAux_runNumber;
    }
    
    // Store a list of runs:
    vectorAdd(m_treeMxAOD->EventInfoAux_runNumber);
    
    // Add events to the cut-flow counter:
    for (int i_c = 0; i_c < m_treeMxAOD->HGamEventInfoAuxDyn_cutFlow; 
    	 i_c++) {
      if (i_c < config->getInt("MxAODCutFlowIndex")) {
	m_cutFlowCounter_Flag[i_c]++;
      }
    }
    // Also cut on events that don't make it to the desired cut:
    if (m_treeMxAOD->HGamEventInfoAuxDyn_cutFlow <
	config->getInt("MxAODCutFlowIndex")) continue;
    
    //---------- Pre-selection Cut ----------//
    if (!m_treeMxAOD->HGamEventInfoAuxDyn_isPassedPreselection) continue;
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")]++;
        
    // Choose the category:
    int category = chooseCategory();
    
    // Calculate an event weight (=1 for data):
    double weight = 1.0;    
    if (config->getBool("MxAODIsMC")) {
      weight = (luminosity * 
		m_treeMxAOD->HGamEventInfoAuxDyn_crossSectionBRfilterEff *
		m_treeMxAOD->HGamEventInfoAuxDyn_weight / nTotEvt);
      //std::cout << "weight = lumi(" << luminosity << ") * xsbrfe("
      //	<< m_treeMxAOD->HGamEventInfoAuxDyn_crossSectionBRfilterEff
      //	<< ") * weight(" << m_treeMxAOD->HGamEventInfoAuxDyn_weight
      //	<< ") / ntotevt(" << nTotEvt << ") = " << weight << std::endl;
    }
    
    // Fill a 2D plot of the photon ID bit:
    for (int i_p = 0; i_p < 2; i_p++) {
      fill2DHistograms(Form("myy_vs_PID%d",i_p+1), category, 
		  (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		  (int)((bool)(*m_treeMxAOD->HGamPhotonsAuxDyn_isTight)[i_p]),
		  weight);
    }
    
    //---------- PID Cut ----------//
    if (!m_treeMxAOD->HGamEventInfoAuxDyn_isPassedPID) continue;
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")+1]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")+1]++;
    
    //---------- Myy Cut ----------//
    if (((config->getStr("AnalysisType")).EqualTo("Scalar") && 
	 m_treeMxAOD->HGamEventInfoAuxDyn_m_yy <= 150000) ||
	((config->getStr("AnalysisType")).Contains("Graviton") && 
	 m_treeMxAOD->HGamEventInfoAuxDyn_m_yy <= 150000)) {
	 //m_treeMxAOD->HGamEventInfoAuxDyn_m_yy <= 200000)) {
      continue;
    }
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")+2]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")+2]++;
    
    //---------- Isolation selection ----------//
    //if (!m_treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy) {
    //continue;
    //}
    
    double isoConstant = 2450.00;
    if ((config->getStr("AnalysisType")).EqualTo("GravitonLoose")) {
      isoConstant = 7000.00;
    }
    
    double pT1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[0];
    double pT2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[1];
    bool isCaloIso1 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[0] <
			    ((0.022 * pT1) + isoConstant));
    bool isCaloIso2 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[1] <
			    ((0.022 * pT2) + isoConstant));
    bool isCaloIsoTight1 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[0] <
			    ((0.022 * pT1) + 2450.00));
    bool isCaloIsoTight2 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[1] <
			    ((0.022 * pT2) + 2450.00));
    
    bool isTrackIso1
      = ((*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[0] < (0.05 * pT1));
    bool isTrackIso2
      = ((*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[1] < (0.05 * pT2));
    
    if (((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	 !(isCaloIso1 && isCaloIso2 && isTrackIso1 && isTrackIso2)) ||
	((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	 !(isCaloIso1 && isCaloIso2 && isTrackIso1 && isTrackIso2)) ||
	((config->getStr("AnalysisType")).EqualTo("GravitonLoose") &&
	 !(isCaloIso1 && isCaloIso2))) {
      continue;
    }

    // Loose-not-tight selection:
    if ((config->isDefined("LooseNotTight") && 
	 config->getBool("LooseNotTight")) && 
	(config->getStr("AnalysisType")).EqualTo("GravitonLoose") && 
	isCaloIsoTight1 && isCaloIsoTight2 && isTrackIso1 && isTrackIso2) {
      continue;
    }
    
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")+3]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")+3]++;
    
    //---------- pT cut ----------//
    double pTRatio1 = (pT1 / m_treeMxAOD->HGamEventInfoAuxDyn_m_yy);
    double pTRatio2 = (pT2 / m_treeMxAOD->HGamEventInfoAuxDyn_m_yy);
    if (((config->getStr("AnalysisType")).EqualTo("Scalar") && 
	 (pTRatio1 <= 0.4 || pTRatio2 <= 0.3)) ||
	((config->getStr("AnalysisType")).Contains("Graviton") && 
	 (pT1 <= 55000 || pT2 <= 55000))) {
      continue;
    }
    /*
    if (((config->getStr("AnalysisType")).EqualTo("Scalar") && 
	 !m_treeMxAOD->HGamEventInfoAuxDyn_isPassedRelPtCutsLowHighMyy) ||
    	((config->getStr("AnalysisType")).Contains("Graviton") && 
    	 !m_treeMxAOD->HGamEventInfoAuxDyn_isPassedlPtCutsExotic)) {
      continue;
    }
    */
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")+4]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")+4]++;
    
    //---------- Select the event ----------//
    if ((config->getStr("AnalysisType")).Contains("Graviton") &&
	!(config->getStr("AnalysisType")).Contains("GravitonLoose") &&
	!(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedExotic && 
	  m_treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      //m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
        
    // Add to event counts:
    countPass++;
    countCate[category]++;
    m_eventsPerRun[m_treeMxAOD->EventInfoAux_runNumber]++;
    eventList << m_treeMxAOD->EventInfoAux_runNumber << " " 
	      << m_treeMxAOD->EventInfoAux_eventNumber << std::endl;
    
    
    //----------------------------------------//
    // Fill the event variable histograms:
    
    
    fillHistograms("z_{vertex} [mm]", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ, weight);
    fillHistograms("m_{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0, weight);
    fillHistograms("p_{T}^{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_pT_yy / 1000.0, weight);
    fillHistograms("cos(#theta*)", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_cosTS_yy, weight);
    
    fillHistograms("N_{Jets}", category,
		   (*m_treeMxAOD->HGamAntiKt4EMTopoJetsAuxDyn_pt).size(),
		   weight);
    fillHistograms("N_{Photons}", category,
		   (*m_treeMxAOD->HGamPhotonsAuxDyn_pt).size(), weight);
    fillHistograms("N_{Electrons}", category,
		   (*m_treeMxAOD->HGamElectronsAuxDyn_pt).size(), weight);
    fillHistograms("N_{Muons}", category,
		   (*m_treeMxAOD->HGamMuonsAuxDyn_pt).size(), weight);
    
    // Fill the photon variable histograms in loop over photons:
    for (int i_p = 0; i_p < 2; i_p++) {
      fillHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p] / 1000.0, 
		     weight);
      fillHistograms(Form("#eta(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p], weight);
      fillHistograms(Form("#eta_{S2}(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_eta_s2)[i_p], weight);
      fillHistograms(Form("#phi(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[i_p], weight);
      fillHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1),
		     category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[i_p] / 1000.0,
		     weight);
      fillHistograms(Form("E_{T}^{cone40}(#gamma_{%d}) [GeV]",i_p+1), 
		     category,
		     ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[i_p]
		      / 1000.0), weight);
      fillHistograms(Form("Conversion type #gamma_{%d}",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[i_p], 
		     weight);
      fillHistograms(Form("gain_{max-E} #gamma_{%d}",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_gain)[i_p],
		     weight);
      
      // Fill 2D plots: 
      fill2DHistograms(Form("timing_maxEcell%d",i_p+1), category,
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_time)[i_p],
		       weight);
      fill2DHistograms(Form("myy_vs_pt%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p] / 1000.0),
		       weight);
      fill2DHistograms(Form("myy_vs_fabseta%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       fabs((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]),
		       weight);
      fill2DHistograms(Form("myy_vs_phi%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[i_p]), weight);
      fill2DHistograms(Form("myy_vs_conv%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[i_p]),
		       weight);
      fill2DHistograms(Form("myy_vs_gain%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_gain)[i_p],
		       weight);
      fill2DHistograms(Form("myy_vs_trkiso%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[i_p]/1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_caloiso%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[i_p]
		       / 1000.0, weight);
      
      // Shower-shape variables:
      fill2DHistograms(Form("myy_vs_weta1_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_weta1)[i_p], weight);
      fill2DHistograms(Form("myy_vs_weta2_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_weta2)[i_p], weight);
      fill2DHistograms(Form("myy_vs_wtots1_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_wtots1)[i_p], weight);
      fill2DHistograms(Form("myy_vs_e277_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_e277)[i_p] / 1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_relEreso_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_relEreso)[i_p], weight);
      fill2DHistograms(Form("myy_vs_Eratio%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_Eratio)[i_p], weight);
      fill2DHistograms(Form("myy_vs_Reta%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_Reta)[i_p], weight);
      fill2DHistograms(Form("myy_vs_f1_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_f1)[i_p], weight);
      fill2DHistograms(Form("myy_vs_Rhad_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_Rhad)[i_p], weight);
      fill2DHistograms(Form("myy_vs_Rhad1_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_Rhad1)[i_p], weight);
      fill2DHistograms(Form("myy_vs_Rphi_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_Rphi)[i_p], weight); 
      fill2DHistograms(Form("myy_vs_DeltaE_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_DeltaE)[i_p] / 1000.0, 
		       weight);
      fill2DHistograms(Form("myy_vs_fracs1_%d",i_p+1), nCategories, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_fracs1)[i_p], weight);
                  
      /////
      
      // Ratio of sampling layers:
      fill2DHistograms(Form("myy_vs_rawcl_ratioEs1Es2_%d",i_p+1), category, 
		      (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		      (*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_ratioEs1Es2)[i_p], 
		       weight);
      
      // Cluster deposits:
      fill2DHistograms(Form("myy_vs_clES0_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es0)[i_p] / 1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_clES1_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es1)[i_p] / 1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_clES2_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es2)[i_p] / 1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_clES3_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es3)[i_p] / 1000.0,
		       weight);
      
      // Raw cluster deposits:
      fill2DHistograms(Form("myy_vs_rawclES0_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es0)[i_p]/1000.0,
		       weight);
      fill2DHistograms(Form("myy_vs_rawclES1_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es1)[i_p]/1000.0, 
		       weight);
      fill2DHistograms(Form("myy_vs_rawclES2_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es2)[i_p]/1000.0, 
		       weight);
      fill2DHistograms(Form("myy_vs_rawclES3_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       (*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es3)[i_p]/1000.0, 
		       weight);
      
      // ratio of cluster E to raw cluster E:
      fill2DHistograms(Form("myy_vs_ratioES0_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es0)[i_p] /
			(*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es0)[i_p]),
		       weight);
      fill2DHistograms(Form("myy_vs_ratioES1_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es1)[i_p] /
			(*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es1)[i_p]),
		       weight);
      fill2DHistograms(Form("myy_vs_ratioES2_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es2)[i_p] /
			(*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es2)[i_p]), 
		       weight);
      fill2DHistograms(Form("myy_vs_ratioES3_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es3)[i_p] /
			(*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es3)[i_p]),
		       weight);
      
      // Ratio of cluster E to photon E
      fill2DHistograms(Form("myy_vs_clES0overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es0)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_clES1overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es1)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_clES2overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es2)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_clES3overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_cl_Es3)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      
      // Ratio of raw cluster E to photon E
      fill2DHistograms(Form("myy_vs_rawES0overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es0)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_rawES1overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es1)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_rawES2overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es2)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
      fill2DHistograms(Form("myy_vs_rawES3overE_%d",i_p+1), category, 
		       (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		       ((*m_treeMxAOD->HGamPhotonsAuxDyn_rawcl_Es3)[i_p] /
			(cosh((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]) * 
			 (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p])), weight);
    }
    fill2DHistograms("myy_vs_deta", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		     fabs((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0] - 
			  (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1]), weight);
    fill2DHistograms("myy_vs_phi", category, 
    		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		     fabs((*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[0] - 
			  (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[1]), weight);
    fill2DHistograms("myy_vs_zvtx", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		     m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ, weight);
    fill2DHistograms("myy_vs_npv", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		     m_treeMxAOD->HGamEventInfoAuxDyn_numberOfPrimaryVertices,
		     weight);
    
    // Additional event-level variables:
    fill2DHistograms("myy_vs_actualIPX", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0),
		     m_treeMxAOD->EventInfoAux_actualInteractionsPerCrossing,
		     weight);
    fill2DHistograms("myy_vs_avgIPX", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->EventInfoAux_averageInteractionsPerCrossing,
		     weight);
    fill2DHistograms("myy_vs_beamPosX", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosX, weight);  
    fill2DHistograms("myy_vs_beamPosY", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosY, weight);
    fill2DHistograms("myy_vs_beamPosZ", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosZ, weight);
    fill2DHistograms("myy_vs_beamPosSigmaX", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosSigmaX, weight);
    fill2DHistograms("myy_vs_beamPosSigmaY", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosSigmaY, weight);
    fill2DHistograms("myy_vs_beamPosSigmaZ", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosSigmaZ, weight);
    fill2DHistograms("myy_vs_beamPosSigmaXY", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     1000.0 * m_treeMxAOD->EventInfoAux_beamPosSigmaXY, weight);
    fill2DHistograms("myy_vs_bunchDistanceFromFront", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->EventInfoAuxDyn_bunchDistanceFromFront,
		     weight);
    fill2DHistograms("myy_vs_bunchDistanceFromFront", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->EventInfoAuxDyn_bunchGapBeforeTrain, weight);
    fill2DHistograms("myy_vs_zSelected", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ, weight);
    fill2DHistograms("myy_vs_zHardest", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_hardestVertexZ, weight);
    fill2DHistograms("myy_vs_zDiff", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_hardestVertexZ -
		      m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ),weight);
    fill2DHistograms("myy_vs_mu", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_mu, weight);
    fill2DHistograms("myy_vs_pT_hard", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_pT_hard / 1000.0, weight);
    fill2DHistograms("myy_vs_met_TST", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_met_TST / 1000.0, weight);
    fill2DHistograms("myy_vs_sumet_TST", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_sumet_TST/1000.0, weight);
    fill2DHistograms("myy_vs_met_phi", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_phi_TST, weight);
    fill2DHistograms("myy_vs_costheta", category, 
		     (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0), 
		     m_treeMxAOD->HGamEventInfoAuxDyn_cosTS_yy, weight);
    
    // Output file to use in neural network study:
    outputForGAN << 0.25*TMath::Log(0.001*(*m_treeMxAOD
					   ->HGamPhotonsAuxDyn_pt)[0])-1.3
		 << " " 
		 << 0.25*TMath::Log(0.001*(*m_treeMxAOD
					   ->HGamPhotonsAuxDyn_pt)[1])-1.3
		 << " " 
		 << ((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0] / 4.0) << " " 
		 << ((*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1] / 4.0) << " " 
		 << ((*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[0])/(2*TMath::Pi())
		 << " " 
		 << ((*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[1])/(2*TMath::Pi())
		 << " " 
		 << 0.5*TMath::Log(0.001*m_treeMxAOD
				   ->HGamEventInfoAuxDyn_m_yy)-3
		 << std::endl;
    
  }// End of loop over events
  eventList.close();
  outputForGAN.close();
  std::cout << "End of event loop" << std::endl;

  // Create output file for histograms:
  TString histogramFileName = Form("%s/histogramFile_%s.root", 
				   m_outputDir.Data(), m_categoryName.Data());
  TFile *histogramFile = new TFile(histogramFileName, "RECREATE");
  
  // Plot and save all 1D histograms:
  std::cout << "Save 1D histograms" << std::endl;
  for (std::map<TString,TH1F*>::iterator histIter = m_histograms.begin(); 
       histIter != m_histograms.end(); histIter++) {
    // Plot individual histogram:
    plotHistogram(histIter->first, config->getStr("AnalysisType"),
		  config->getStr("ATLASLabel"), luminosity / 1000.0);
    // Plot comparison of histograms in each category:
    if (m_categoryName.EqualTo("MassCate")) {
      plotComparisonHist(histIter->first, nCategories, 
			 config->getStr("AnalysisType"),
			 config->getStr("ATLASLabel"), luminosity / 1000.0);
    }
    histIter->second
      ->Write(Form("hist_%s", (formatHistName(histIter->first)).Data()));
  }
  
  // Plot and save all 2D histograms:
  std::cout << "Save 2D histograms." << std::endl;
  for (std::map<TString,TH2D*>::iterator histIter2D = m_histograms2D.begin(); 
       histIter2D != m_histograms2D.end(); histIter2D++) {
    // Plot individual histogram:
    plot2DHistogram(histIter2D->first, config->getStr("AnalysisType"),
		    config->getStr("ATLASLabel"),
		    config->getNum("AnalysisLuminosity")/1000.0);
    histIter2D->second
      ->Write(Form("hist_%s", (formatHistName(histIter2D->first)).Data()));
  }
  
  // Write output file:
  std::cout << "Writing output file" << std::endl;
  histogramFile->Close();
  
  
  
  //----------------------------------------//
  // Remove files that were copied:
  if (config->getBool("MakeLocalMxAODCopies")) removeLocalFileCopies(fileNames);
  
  // Print a list of runs used:
  std::cout << "List of runs in the MxAOD:" << std::endl;
  for (int i_r = 0; i_r < (int)m_runList.size(); i_r++) {
    std::cout << "\t" << m_runList[i_r] << "\t" 
	      << m_eventsPerRun[m_runList[i_r]] << std::endl;
  }
  
  // Detailed summary:
  std::cout << "\nStudyData: Printing cut-flow from MxAOD and offline analysis."
	    << std::endl;
  for (int i_c = 0; i_c < (int)m_cutNames.size(); i_c++) {
    std::cout << "\t" << i_c << "\t" << m_cutNames[i_c] << " \t" 
	      << m_cutFlowCounter_Hist[i_c] << "\t" 
	      << m_cutFlowCounter_Flag[i_c] << std::endl;
  }
  
  // Category summary:
  std::ofstream cateFile(Form("%s/categoryCount_%s.txt", m_outputDir.Data(),
			      m_categoryName.Data()));
  std::cout << "\nThere were " << countPass << " in the final sample."
	    << std::endl;
  cateFile << "\nInclusive " << countPass << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::cout << "\t" << nameCategory(i_c)<< "\t" << countCate[i_c] 
	      << std::endl;
    cateFile << nameCategory(i_c)<< "\t" << countCate[i_c] << std::endl;
  }
  cateFile.close();
  
  // Print location of output ROOT file:
  std::cout << "\nCreated ROOT file with histograms: " << histogramFileName 
	    << std::endl;
  
  delete m_treeMxAOD;
  delete chain;
  delete config;
  delete histogramFile;  
  return 0;
}
