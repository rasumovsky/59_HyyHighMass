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
HGammaMxAOD *m_treeMxAOD;
TString m_categoryName;

int m_cutFlowCounter_Hist[100];
int m_cutFlowCounter_Flag[100];
std::vector<TString> m_cutNames;

std::vector<int> m_runList;
std::map<int,int> m_eventsPerRun;

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
    if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy >= 600000 &&
	m_treeMxAOD->HGamEventInfoAuxDyn_m_yy < 700000) return 0;
    else if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy >= 700000 &&
	     m_treeMxAOD->HGamEventInfoAuxDyn_m_yy < 800000) return 1;
    else if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy >= 800000 &&
	m_treeMxAOD->HGamEventInfoAuxDyn_m_yy < 900000) return 2;
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
  if (histName.Contains("BarrelBarrel")) return "2 Barrel #gamma's";
  else if (histName.Contains("BarrelEndcap")) return "Barrel-Endcap";
  else if (histName.Contains("EndcapBarrel")) return "Endcap-Barrel";
  else if (histName.Contains("EndcapEndcap")) return "2 Endcap #gamma's";
  else if (histName.Contains("UnconvUnconv")) return "2 Unconverted #gamma's";
  else if (histName.Contains("UnconvConv")) return "Unconverted-Converted";
  else if (histName.Contains("ConvUnconv")) return "Converted-Unconverted";
  else if (histName.Contains("ConvConv")) return "2 Converted #gamma's";
  else if (histName.Contains("mgg_600_700")) {
    return "600 GeV < m_{#gamma#gamma} < 700 GeV";
  }
  else if (histName.Contains("mgg_700_800")) {
    return "700 GeV < m_{#gamma#gamma} < 800 GeV";
  }
  else if (histName.Contains("mgg_800_inf")) {
    return "m_{#gamma#gamma} > 800 GeV";
  }
  else if (histName.Contains("mgg_rest")) return "mgg_rest";
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
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      = new TH1F(histName, histName, nBins, xMin, xMax);
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
   Fill inclusive and categorized histograms.
   @param histName - The name of the histogram.
   @param category -The index of the category into which this event falls.
   @param value - The value to fill into the histogram for this event.
*/
void fillHistograms(TString histName, int category, double value) {
  m_histograms[histName]->Fill(value);
  TString cateName = nameCategory(category);
  m_histograms[Form("%s_%s", histName.Data(), cateName.Data())]->Fill(value);
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
  //gPad->SetLogy();
  
  // Create a legend:
  TLegend leg(0.3, 0.75, 0.5, 0.93);
  leg.SetTextFont(42); 
  leg.SetTextSize(0.03);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  
  // Plot each histogram and add legend entry:
  Color_t lineColors[5] = {kBlue+1, kMagenta+1, kGreen+1, kCyan+1, kRed+1};
  for (int i_c = 0; i_c < nCategories; i_c++) {
    if (i_c > 2) continue;
    TString currCateName = nameCategory(i_c);
    TString currHistName = Form("%s_%s", histName.Data(), 
				currCateName.Data());
    m_histograms[currHistName]
      ->Scale(1.0 / m_histograms[currHistName]->Integral());
    m_histograms[currHistName]->SetLineColor(lineColors[i_c]);
    m_histograms[currHistName]->SetMarkerColor(lineColors[i_c]);
    m_histograms[currHistName]->SetMarkerStyle(21 + i_c);
    if (i_c == 0) m_histograms[currHistName]->Draw("E1");
    else m_histograms[currHistName]->Draw("E1SAME");
    leg.AddEntry(m_histograms[currHistName],
		 histToCateName(nameCategory(i_c)), "LEP");
  }
  histName = formatHistName(histName);
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.04);
  t.DrawLatex(0.6, 0.87, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.04);
  t.DrawLatex(0.72, 0.87, label);
  t.DrawLatex(0.6, 0.82, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",lumi));
  if (ana.Contains("Scalar")) t.DrawLatex(0.6, 0.77, "Spin-0 Selection");
  else if (ana.Contains("GravitonLoose")) {
    t.DrawLatex(0.6, 0.77, "Spin-2 Loose Iso.");
  }
  else t.DrawLatex(0.6, 0.77, "Spin-2 Selection");
  t.DrawLatex(0.6, 0.72, histToCateName(histName));
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
  TString options = argv[2];
  
  // Check that output directory exists:
  m_outputDir = Form("%s/%s/StudyData", (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
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
  defineHistograms("N_{Jets}", nCategories, 5, -0.5, 4.5);
  defineHistograms("N_{Photons}", nCategories, 5, -0.5, 4.5);
  defineHistograms("N_{Electrons}", nCategories, 5, -0.5, 4.5);
  defineHistograms("N_{Muons}", nCategories, 5, -0.5, 4.5);
  
  // Define histograms of per-photon variables:
  for (int i_p = 0; i_p < 2; i_p++) {
    defineHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0.0, 1000.0);
    defineHistograms(Form("#eta(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -2.5, 2.5);
    defineHistograms(Form("#eta_{S2}(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -2.5, 2.5);
    defineHistograms(Form("#phi(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -3.141, 3.141);
    defineHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1), 
		     nCategories, nBins, 0.0, 50.0);
    defineHistograms(Form("E_{T}^{cone40}(#gamma_{%d}) [GeV]",i_p+1), 
		     nCategories, nBins, -5.0, 95.0);
    defineHistograms(Form("Conversion type #gamma_{%d}",i_p+1),
		     nCategories, 6, -0.5, 5.5);
  
    defineHistograms(Form("Conversion radius #gamma_{%d}",i_p+1),
		     nCategories, nBins, 0.0, 1000.0);
    defineHistograms(Form("isEMTight #gamma_{%d}",i_p+1),
		     nCategories, 2, -0.5, 1.5);
    defineHistograms(Form("gain_{max-E} #gamma_{%d}",i_p+1),
		     nCategories, 3, -0.5, 2.5);
     
    defineHistograms(Form("E_{0}^{Raw}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0, 1000);
    defineHistograms(Form("E_{1}^{Raw}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0, 1000);
    defineHistograms(Form("E_{2}^{Raw}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0, 1000);
    defineHistograms(Form("E_{3}^{Raw}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0, 1000);
  }
  
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
  m_treeMxAOD = new HGammaMxAOD(chain, config->getStr("MxAODTag"));
  
  // Count events:
  for (int i_c = 0; i_c < 50; i_c++) {
    m_cutFlowCounter_Hist[i_c] = 0;
    m_cutFlowCounter_Flag[i_c] = 0;
  }
  m_cutNames.clear();
  m_runList.clear();
  TString currFileName = "";
  int countPass = 0;
  int countCate[50] = {0};
  m_eventsPerRun.clear();
  
  int prevRun = 0;
  
  // Also create a text file of all passing events:
  std::ofstream eventList(Form("%s/eventList.txt", m_outputDir.Data()));
  
  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  int nEvents = m_treeMxAOD->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    m_treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    // Add each new file to the cutflow:
    if (!currFileName.EqualTo(chain->GetFile()->GetName())) {
      currFileName = chain->GetFile()->GetName();
      fillCutFlowFromMxAOD(chain->GetFile(), 
			   config->getInt("MxAODCutFlowIndex"));
    }
    
    if (prevRun != (int)(m_treeMxAOD->EventInfoAux_runNumber)) {
      //std::cout << "NEW RUN! " << m_treeMxAOD->EventInfoAux_runNumber 
      //	<< std::endl;
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
    
    if (!m_treeMxAOD->HGamEventInfoAuxDyn_isPassedPreselection) continue;
    m_cutFlowCounter_Hist[config->getInt("MxAODCutFlowIndex")]++;
    m_cutFlowCounter_Flag[config->getInt("MxAODCutFlowIndex")]++;
    
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
    bool isCaloIsoTight1 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[0] <
			    ((0.022 * pT1) + isoConstant));
    bool isCaloIsoTight2 = ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[1] <
			    ((0.022 * pT2) + isoConstant));
    bool isTrackIso1
      = ((*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[0] < (0.05 * pT1));
    bool isTrackIso2
      = ((*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[1] < (0.05 * pT2));
    
    if (((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	 !(isCaloIsoTight1 && isCaloIsoTight2 && isTrackIso1 && isTrackIso2)) ||
	((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	 !(isCaloIsoTight1 && isCaloIsoTight2 && isTrackIso1 && isTrackIso2)) ||
	((config->getStr("AnalysisType")).EqualTo("GravitonLoose") &&
	 !(isCaloIsoTight1 && isCaloIsoTight2))) {
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
        
    /*
      // CHECK FOR PASSING EVENT COMPARISON
      std::cout << "PASSING EVENT: Run = " 
	      << m_treeMxAOD->EventInfoAux_runNumber << ", Event = " 
	      << m_treeMxAOD->EventInfoAux_eventNumber << std::endl;
    */    
    
    // Choose the category:
    int category = chooseCategory();
    
    // Add to event counts:
    countPass++;
    countCate[category]++;
    m_eventsPerRun[m_treeMxAOD->EventInfoAux_runNumber]++;
    eventList << m_treeMxAOD->EventInfoAux_runNumber << " " 
	      << m_treeMxAOD->EventInfoAux_eventNumber << std::endl;

    // Fill the event variable histograms:
    fillHistograms("z_{vertex} [mm]", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ);
    fillHistograms("m_{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0);
    fillHistograms("p_{T}^{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_pT_yy / 1000.0);
    fillHistograms("cos(#theta*)", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_cosTS_yy);
    //std::cout << "CHECK0" << std::endl;
    fillHistograms("N_{Jets}", category,
		   (*m_treeMxAOD->HGamAntiKt4EMTopoJetsAuxDyn_pt).size());
    fillHistograms("N_{Photons}", category,
		   (*m_treeMxAOD->HGamPhotonsAuxDyn_pt).size());
    fillHistograms("N_{Electrons}", category,
		   (*m_treeMxAOD->HGamElectronsAuxDyn_pt).size());
    fillHistograms("N_{Muons}", category,
		   (*m_treeMxAOD->HGamMuonsAuxDyn_pt).size());
    //std::cout << "CHECK1" << std::endl;
    
    // Fill the photon variable histograms in loop over photons:
    for (int i_p = 0; i_p < 2; i_p++) {
      fillHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p] / 1000.0);
      fillHistograms(Form("#eta(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]);
      fillHistograms(Form("#eta_{S2}(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_eta_s2)[i_p]);
      fillHistograms(Form("#phi(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[i_p]);
      fillHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1),
		     category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[i_p] / 1000.0);
      fillHistograms(Form("E_{T}^{cone40}(#gamma_{%d}) [GeV]",i_p+1), 
		     category,
		     ((*m_treeMxAOD->HGamPhotonsAuxDyn_topoetcone40)[i_p]
		      / 1000.0));
      fillHistograms(Form("Conversion type #gamma_{%d}",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[i_p]);
      //std::cout << "CHECK2" << std::endl;
      /*
      fillHistograms(Form("Conversion radius #gamma_{%d}",i_p+1), category, 
      	     (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionRadius)[i_p]);
      fillHistograms(Form("isEMTight #gamma_{%d}",i_p+1), category, 
      	     (*m_treeMxAOD->HGamPhotonsAuxDyn_isTight)[i_p]);
      fillHistograms(Form("gain_{max-E} #gamma_{%d}",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_maxEcell_gain)[i_p]);
      std::cout << "CHECK2.5" << std::endl;
      fillHistograms(Form("E_{0}^{Raw}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_E0_raw)[i_p]);
      fillHistograms(Form("E_{1}^{Raw}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_E1_raw)[i_p]);
      fillHistograms(Form("E_{2}^{Raw}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_E2_raw)[i_p]);
      fillHistograms(Form("E_{3}^{Raw}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_E3_raw)[i_p]);
      */
      //std::cout << "CHECK3" << std::endl;
    }
  }// End of loop over events
  eventList.close();
  
  // Then plot and save all histograms:
  TString histogramFileName = Form("%s/histogramFile_%s.root", 
				   m_outputDir.Data(), m_categoryName.Data());
  TFile *histogramFile = new TFile(histogramFileName, "RECREATE");
  for (std::map<TString,TH1F*>::iterator histIter = m_histograms.begin(); 
       histIter != m_histograms.end(); histIter++) {
    // Plot individual histogram:
    plotHistogram(histIter->first, config->getStr("AnalysisType"),
		  config->getStr("ATLASLabel"),
		  config->getNum("AnalysisLuminosity")/1000.0);
    // Plot comparison of histograms in each category:
    if (m_categoryName.EqualTo("MassCate")) {
      plotComparisonHist(histIter->first, nCategories, 
			 config->getStr("AnalysisType"),
			 config->getStr("ATLASLabel"),
			 config->getNum("AnalysisLuminosity")/1000.0);
    }
    histIter->second
      ->Write(Form("hist_%s", (formatHistName(histIter->first)).Data()));
  }
  histogramFile->Close();
  
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

/*
// DETAILED EVENT CHECKLIST:
std::cout << "\nRun = " << m_treeMxAOD->EventInfoAux_runNumber
<< ", Event = " << m_treeMxAOD->EventInfoAux_eventNumber
<< ", conv1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[0] 
<< ", conv2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[1] 
<< ", eta1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0]
<< ", eta2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1]
<< ", phi1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[0]
<< ", phi2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[1]
<< ", pT1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[0]
<< ", pT2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[1]
<< ", zvtx = " << m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ
<< ", E0,1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E0_raw)[0]
<< ", E0,2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E0_raw)[1]

<< ", E1,1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E1_raw)[0]
<< ", E1,2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E1_raw)[1]

<< ", E2,1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E2_raw)[0]
<< ", E2,2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E2_raw)[1]

<< ", E3,1 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E3_raw)[0]
<< ", E3,2 = " << (*m_treeMxAOD->HGamPhotonsAuxDyn_E3_raw)[1]
/<< ", myy = " << m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0
<< std::endl;
}
*/
