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

/**
   -----------------------------------------------------------------------------
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
    m_histograms[Form("%s_c%d",histName.Data(),i_c)]
      = new TH1F(histName, histName, nBins, xMin, xMax);
    m_histograms[Form("%s_c%d",histName.Data(),i_c)]
      ->GetXaxis()->SetTitle(histName);
    m_histograms[Form("%s_c%d",histName.Data(),i_c)]
      ->GetYaxis()->SetTitle(Form("Entries / %2.2f",binning));
  }
}

/**
   -----------------------------------------------------------------------------
*/
void fillHistograms(TString histName, int category, double value) {
  m_histograms[histName]->Fill(value);
  m_histograms[Form("%s_c%d", histName.Data(), category)]->Fill(value);
}

/**
   -----------------------------------------------------------------------------
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
  
  // Check that output directory exists:
  m_outputDir = Form("%s/%s/StudyData", (config->getStr("MasterOutput")).Data(),
		     (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Define histograms:
  m_histograms.clear();
  int nBins = 50;
  int nCategories = 1;
 
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
    defineHistograms(Form("p_{T}^{#DeltaR<20 GeV}(#gamma_{%d}) [GeV]",i_p+1), 
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
  HGammaMxAOD *treeMxAOD = new HGammaMxAOD(chain);
    
  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  int passingEvents = 0;
  int nEvents = treeMxAOD->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    // Select the event:
    if ((config->getStr("AnalysisType")).EqualTo("Graviton") &&
	!(treeMxAOD->HGamEventInfoAuxDyn_isPassedExotic && 
	  treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  treeMxAOD->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(treeMxAOD->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
    
    passingEvents++;
    
    // Choose the category:
    int category = 0;
    
    // Fill the event variable histograms:
    fillHistograms("z_{vertex} [mm]", category,
		   treeMxAOD->HGamEventInfoAuxDyn_selectedVertexZ);
    fillHistograms("m_{#gamma#gamma} [GeV]", category, 
		   treeMxAOD->HGamEventInfoAuxDyn_m_yy / 1000.0);
    fillHistograms("p_{T}^{#gamma#gamma} [GeV]", category, 
		   treeMxAOD->HGamEventInfoAuxDyn_pT_yy / 1000.0);
    fillHistograms("cos(#theta*)", category,
		   treeMxAOD->HGamEventInfoAuxDyn_cosTS_yy);
    
    // Fill the photon variable histograms in loop over photons:
    for (int i_p = 0; i_p < 2; i_p++) {
      fillHistograms(Form("p_{T}(#gamma_{%d}) [GeV]",i_p+1), category, 
		     (*treeMxAOD->HGamPhotonsAuxDyn_pt)[i_p] / 1000.0);
      fillHistograms(Form("#eta(#gamma_{%d})",i_p+1), category, 
		     (*treeMxAOD->HGamPhotonsAuxDyn_eta)[i_p]);
      fillHistograms(Form("#phi(#gamma_{%d})",i_p+1), category, 
		     (*treeMxAOD->HGamPhotonsAuxDyn_phi)[i_p]);
      fillHistograms(Form("p_{T}^{#DeltaR<20 GeV}(#gamma_{%d}) [GeV]",i_p+1),
		     category, 
		     (*treeMxAOD->HGamPhotonsAuxDyn_ptcone20)[i_p] / 1000.0);
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
  
  // Quick summary:
  std::cout << "\nThere were " << passingEvents << " in the final sample."
	    << std::endl;
  std::cout << "Created ROOT file with histograms: " << histogramFileName 
	    << std::endl;

  delete treeMxAOD;
  delete chain;
  delete config;
  delete histogramFile;
  
  return 0;
}
