////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  StudyTwoDatasets.cxx                                                      //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 01/05/2016                                                          //
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
HGammaMxAOD *m_treeMxAOD1;
HGammaMxAOD *m_treeMxAOD2;
TString m_categorization;

/**
   -----------------------------------------------------------------------------
   Select the proper category for the current event in the HGammaMxAOD. The 
   categorization to use is based on the global variable m_categorization.
   @param tree - The MxAOD to use for the current categorization.
   @return - The category index.
*/
int chooseCategory(HGammaMxAOD *tree) {

  // Eta categorization:
  if (m_categorization.EqualTo("EtaCate")) {
    double barrelEnd = 1.37;
    double eta1 = (*tree->HGamPhotonsAuxDyn_eta)[0];
    double eta2 = (*tree->HGamPhotonsAuxDyn_eta)[1];
    if (fabs(eta1) < barrelEnd && fabs(eta2) < barrelEnd) return 0;
    else if (fabs(eta1) < barrelEnd) return 1;
    else if (fabs(eta2) < barrelEnd) return 2;
    else return 3;
  }

  // Conversion categorization:
  else if (m_categorization.EqualTo("ConversionCate")) {
    int conv1 = (*tree->HGamPhotonsAuxDyn_conversionType)[0];
    int conv2 = (*tree->HGamPhotonsAuxDyn_conversionType)[1];
    if (conv1 == 0 && conv2 == 0) return 0;
    else if (conv1 == 0) return 1;
    else if (conv2 == 0) return 2;
    else return 3;
  }
 
  // Exit because inputs are unknown:
  else {
    std::cout << "StudyTwoDatasets: ERROR! Categorization not found"
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the name for the current category. The categorization to use is based on 
   the global variable m_categorization.
   @param categoryIndex - The index of the current category.
   @return - The category name.
*/
TString nameCategory(int categoryIndex) {
  if (m_categorization.EqualTo("EtaCate")) {
    if (categoryIndex == 0) return "BarrelBarrel";
    else if (categoryIndex == 1) return "BarrelEndcap";
    else if (categoryIndex == 2) return "EndcapBarrel";
    else return "EndcapEndcap";
  }
  else if (m_categorization.EqualTo("ConversionCate")) {
    if (categoryIndex == 0) return "UnconvUnconv";
    else if (categoryIndex == 1) return "UnconvConv";
    else if (categoryIndex == 2) return "ConvUnconv";
    else return "ConvConv";
  }
  
  // Exit because inputs are unknown:
  else {
    std::cout << "StudyTwoDatasets: ERROR! Categorization not found" 
	      << std::endl;
    exit(0);
  }
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
  m_histograms[histName]
    ->GetYaxis()->SetTitle(Form("Entries / %2.2f",binning));
  
  // Create categorized histograms:
  for (int i_c = 0; i_c < nCategories; i_c++) {
    TString cateName = nameCategory(i_c);
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      = new TH1F(histName, histName, nBins, xMin, xMax);
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetXaxis()->SetTitle(histName);
    m_histograms[Form("%s_%s",histName.Data(), cateName.Data())]
      ->GetYaxis()->SetTitle(Form("Entries / %2.2f",binning));
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
*/
void plotHistogram(TString histName) {
  TCanvas *can = new TCanvas("can", "can");
  can->cd();
  gPad->SetLogy();
  m_histograms[histName]->Draw("E1");
  histName = formatHistName(histName);
  TString printName = Form("%s/plot_%s.eps", m_outputDir.Data(),
			   histName.Data());
  can->Print(printName);
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @param sample - The index of the sample (1 or 2).
   @return - An updated list of file names.
*/
std::vector<TString> makeLocalFileCopies(std::vector<TString> fileNames,
					 int sample) {
  std::cout << "StudyTwoDatasets: Making local copies of inputs."
	    << std::endl;
  std::vector<TString> result; result.clear();
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    TString newName = Form("tempFile%d_s%d.root", i_f, sample);
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
  std::cout << "StudyTwoDatasets: Removing local copies of inputs."
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
  if (index%10 == 0) {
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
  m_outputDir = Form("%s/%s/StudyTwoDatasets",
		     (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Define histograms:
  m_histograms.clear();
  int nBins = 50;
  int nCategories = config->getInt("MxAODNCategories");
  m_categorization = config->getStr("MxAODCategorization");
  
  defineHistograms("#Delta(z_{vertex}) [mm]", nCategories, nBins, -20, 20);
  defineHistograms("#Delta(m_{#gamma#gamma}) [GeV]",
		   nCategories, nBins, -10, 10);
  defineHistograms("#Delta(p_{T}^{#gamma#gamma}) [GeV]",
		   nCategories, nBins, -10, 10);
  defineHistograms("#Delta(cos(#theta*))", nCategories, nBins, -0.2, 0.2);
  
  // Loop over photons:
  for (int i_p = 0; i_p < 2; i_p++) {
    defineHistograms(Form("#Delta(p_{T}(#gamma_{%d})) [GeV]",i_p+1),
		     nCategories, nBins, -10.0, 10.0);
    defineHistograms(Form("#Delta(#eta(#gamma_{%d}))",i_p+1),
		     nCategories, nBins, -0.5, 0.5);
    defineHistograms(Form("#Delta(#phi(#gamma_{%d}))",i_p+1),
		     nCategories, nBins, -0.5, 0.5);
    defineHistograms(Form("#Delta(p_{T}^{cone20}(#gamma_{%d})) [GeV]",
			  i_p+1), nCategories, nBins, -5.0, 5.0);
  }
  
  // Make local copies of files if requested, to improve speed:
  std::vector<TString> fileNames1 = config->getStrV("MxAODsForData1");
  std::vector<TString> fileNames2 = config->getStrV("MxAODsForData2");
  if (config->getBool("MakeLocalMxAODCopies")) {
    fileNames1 = makeLocalFileCopies(fileNames1, 1);
    fileNames2 = makeLocalFileCopies(fileNames2, 2);
  }
  
  // Create TChain of input files for sample 1:
  TChain *chain1 = new TChain(config->getStr("MxAODTreeName"));
  for (int i_f1 = 0; i_f1 < (int)fileNames1.size(); i_f1++) {
    chain1->AddFile(fileNames1[i_f1]);
  }
  m_treeMxAOD1 = new HGammaMxAOD(chain1);
  
  // Create TChain of input files for sample 2:
  TChain *chain2 = new TChain(config->getStr("MxAODTreeName"));
  for (int i_f2 = 0; i_f2 < (int)fileNames2.size(); i_f2++) {
    chain2->AddFile(fileNames2[i_f2]);
  }
  m_treeMxAOD2 = new HGammaMxAOD(chain2);
  
  // Count common events:
  int countPass = 0;
  int countCate[100] = {0};

  //--------------------------------------//
  // Loop over events from dataset 1:  
  int nEvents1 = m_treeMxAOD1->fChain->GetEntries();
  std::cout << "There are " << nEvents1 << " events to process." << std::endl;
  for (int index1 = 0; index1 < nEvents1; index1++) {

    // Load event from MxAOD 1:
    m_treeMxAOD1->fChain->GetEntry(index1);
    printProgressBar(index1, nEvents1);
    
    // Select the event:
    if ((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	!(m_treeMxAOD1->HGamEventInfoAuxDyn_isPassedExotic && 
	  m_treeMxAOD1->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  m_treeMxAOD1->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(m_treeMxAOD1->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       m_treeMxAOD1->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
    
    //--------------------------------------//
    // Loop over events from dataset 2:
    int nEvents2 = m_treeMxAOD2->fChain->GetEntries();
    for (int index2 = 0; index2 < nEvents2; index2++) {
      
      m_treeMxAOD2->fChain->GetEntry(index2);
      
      if (m_treeMxAOD1->EventInfoAuxDyn_eventNumber == 
	  m_treeMxAOD2->EventInfoAuxDyn_eventNumber) {
	
	
	
	// Select the event:
	if ((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	    !(m_treeMxAOD2->HGamEventInfoAuxDyn_isPassedExotic && 
	      m_treeMxAOD2->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	      m_treeMxAOD2->HGamEventInfoAuxDyn_m_yy > 200000)) {
	  continue;
	}
	else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
		 !(m_treeMxAOD2->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
		   m_treeMxAOD2->HGamEventInfoAuxDyn_m_yy > 150000)) {
	  continue;
	}
	
	//// Choose the category:
	int category1 = chooseCategory(m_treeMxAOD2);
	int category2 = chooseCategory(m_treeMxAOD2);
	if (category1 != category2) continue;
	int category = category1;
	
	// Count events that pass both selections:
	countPass++;
	countCate[category]++;
	
	
	// Fill the event variable histograms:
	fillHistograms("#Delta(z_{vertex}) [mm]", category,
		       (m_treeMxAOD2->HGamEventInfoAuxDyn_selectedVertexZ - 
			m_treeMxAOD1->HGamEventInfoAuxDyn_selectedVertexZ));
	fillHistograms("#Delta(m_{#gamma#gamma}) [GeV]", category, 
		       (m_treeMxAOD2->HGamEventInfoAuxDyn_m_yy - 
			m_treeMxAOD1->HGamEventInfoAuxDyn_m_yy) / 1000.0);
	fillHistograms("#Delta(p_{T}^{#gamma#gamma}) [GeV]", category, 
		       (m_treeMxAOD2->HGamEventInfoAuxDyn_pT_yy - 
			m_treeMxAOD1->HGamEventInfoAuxDyn_pT_yy) / 1000.0);
	fillHistograms("#Delta(cos(#theta*))", category,
		       (m_treeMxAOD2->HGamEventInfoAuxDyn_cosTS_yy - 
			m_treeMxAOD1->HGamEventInfoAuxDyn_cosTS_yy));
	
	// Fill the photon variable histograms in loop over photons:
	for (int i_p = 0; i_p < 2; i_p++) {
	  fillHistograms(Form("#Delta(p_{T}(#gamma_{%d})) [GeV]",i_p+1),
			 category, 
			 ((*m_treeMxAOD2->HGamPhotonsAuxDyn_pt)[i_p] - 
			  (*m_treeMxAOD1->HGamPhotonsAuxDyn_pt)[i_p]) / 1000.0);
	  fillHistograms(Form("#Delta(#eta(#gamma_{%d}))",i_p+1), category, 
			 ((*m_treeMxAOD2->HGamPhotonsAuxDyn_eta)[i_p] - 
			  (*m_treeMxAOD1->HGamPhotonsAuxDyn_eta)[i_p]));
	  fillHistograms(Form("#Delta(#phi(#gamma_{%d}))",i_p+1), category, 
			 ((*m_treeMxAOD2->HGamPhotonsAuxDyn_phi)[i_p] - 
			  (*m_treeMxAOD1->HGamPhotonsAuxDyn_phi)[i_p]));
	  fillHistograms(Form("#Delta(p_{T}^{cone20}(#gamma_{%d})) [GeV]",
			      i_p+1), category, 
			 ((*m_treeMxAOD2->HGamPhotonsAuxDyn_ptcone20)[i_p] - 
			  (*m_treeMxAOD1->HGamPhotonsAuxDyn_ptcone20)[i_p]) /
			 1000.0);
	}
		
	break;
      }
    }
  }
  
  // Then plot and save all histograms:
  TString histogramFileName = Form("%s/histogramFile.root", m_outputDir.Data());
  TFile *histogramFile = new TFile(histogramFileName, "RECREATE");
  for (std::map<TString,TH1F*>::iterator histIter = m_histograms.begin(); 
       histIter != m_histograms.end(); histIter++) {
    plotHistogram(histIter->first);
    histIter->second
      ->Write(Form("hist_%s", (formatHistName(histIter->first)).Data()));
  }
  histogramFile->Close();
  
  // Remove files that were copied:
  if (config->getBool("MakeLocalMxAODCopies")) {
    removeLocalFileCopies(fileNames1);
    removeLocalFileCopies(fileNames2);
  }
  
  std::ofstream textFile(Form("%s/eventCount.txt", m_outputDir.Data()));
  
  // Quick summary:
  std::cout << "\nThere were " << countPass << " in the final sample."
	    << std::endl;
  textFile << "\nThere were " << countPass << " in the final sample." 
	   << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::cout << "\t" << nameCategory(i_c)<< "\t" << countCate[i_c] 
	      << std::endl;
    textFile << "\t" << nameCategory(i_c)<< "\t" << countCate[i_c] << std::endl;
  }
  textFile.close();
  
  std::cout << "\nCreated ROOT file with histograms: " << histogramFileName 
	    << std::endl;
  
  delete m_treeMxAOD1;
  delete m_treeMxAOD2;
  delete chain1;
  delete chain2;
  delete config;
  delete histogramFile;
  
  return 0;
}
