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
  
  // Exit because inputs are unknown:
  else {
    std::cout << "StudyData: ERROR! Categorization not found" << std::endl;
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
  
  // Define histograms:
  m_histograms.clear();
  int nBins = 50;
  int nCategories = config->getInt("MxAODNCategories");
  m_categoryName = config->getStr("MxAODCategorization");
  
  defineHistograms("z_{vertex} [mm]", nCategories, nBins, -150, 150);
  defineHistograms("m_{#gamma#gamma} [GeV]", nCategories, nBins, 0, 2000);
  defineHistograms("p_{T}^{#gamma#gamma} [GeV]", nCategories, nBins, 0, 1000);
  defineHistograms("cos(#theta*)", nCategories, nBins, 0, 1);
  
  // Loop over photons:
  for (int i_p = 0; i_p < 2; i_p++) {
    defineHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1),
		     nCategories, nBins, 0.0, 1000.0);
    defineHistograms(Form("#eta(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -2.5, 2.5);
    defineHistograms(Form("#phi(#gamma_{%d})",i_p+1),
		     nCategories, nBins, -3.141, 3.141);
    defineHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1), 
		     nCategories, nBins, 0.0, 50.0);
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
  m_treeMxAOD = new HGammaMxAOD(chain);
  
  // Count events:
  int countPass = 0;
  int countCate[100] = {0};

  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  
  int nEvents = m_treeMxAOD->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    m_treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    /*
    if ((m_treeMxAOD->EventInfoAux_runNumber == 280464 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 187603384) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 280862 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 128240362) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 281317 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 144044837) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 282631 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 444615930) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 282712 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 1400405523) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 283155 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 212630049) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 283429 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 1931043956) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 283608 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 458220275) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 284213 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 151967530) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 284484 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 1014170264)) {
      
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
	//<< ", myy = " << m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0
		<< std::endl;
    }
    */

    // Select the event:
    if ((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	!(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedExotic && 
	  m_treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
    
    // Choose the category:
    int category = chooseCategory();
    
    // Add to event counts:
    countPass++;
    countCate[category]++;

    //std::cout << "run = " << m_treeMxAOD->EventInfoAux_runNumber
    //	      << " event = " << m_treeMxAOD->EventInfoAux_eventNumber
    //	      << std::endl;
    // Print out information on specific events:    
    if ((m_treeMxAOD->EventInfoAux_runNumber == 280319 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 844110200) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 280319 &&
	 m_treeMxAOD->EventInfoAux_eventNumber == 1844702370) || 
	(m_treeMxAOD->EventInfoAux_runNumber == 280862 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 2428284929) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 280231 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 434101590) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 280319 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 844110200) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 280319 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 1844702370) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 280862 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 2428284929) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 280950 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 202310071) ||
	(m_treeMxAOD->EventInfoAux_runNumber == 284213 && 
	 m_treeMxAOD->EventInfoAux_eventNumber == 3096852034)) {
    
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

	//<< ", myy = " << m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0
		<< std::endl;
    }
    
    /*
    if (m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0 > 1600) {
      std::cout << "\nSPECIAL   "
		<< "Run = " << m_treeMxAOD->EventInfoAux_runNumber
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
		<< ", myy = " << m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0
		<< std::endl;
    }
    */

    // Fill the event variable histograms:
    fillHistograms("z_{vertex} [mm]", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ);
    fillHistograms("m_{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0);
    fillHistograms("p_{T}^{#gamma#gamma} [GeV]", category, 
		   m_treeMxAOD->HGamEventInfoAuxDyn_pT_yy / 1000.0);
    fillHistograms("cos(#theta*)", category,
		   m_treeMxAOD->HGamEventInfoAuxDyn_cosTS_yy);
    
    // Fill the photon variable histograms in loop over photons:
    for (int i_p = 0; i_p < 2; i_p++) {
      fillHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p] / 1000.0);
      fillHistograms(Form("#eta(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]);
      fillHistograms(Form("#phi(#gamma_{%d})",i_p+1), category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_phi)[i_p]);
      fillHistograms(Form("p_{T}^{cone20}(#gamma_{%d}) [GeV]",i_p+1),
		     category, 
		     (*m_treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[i_p] / 1000.0);
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
  if (config->getBool("MakeLocalMxAODCopies")) removeLocalFileCopies(fileNames);
  
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
  
  delete m_treeMxAOD;
  delete chain;
  delete config;
  delete histogramFile;
  
  return 0;
}
