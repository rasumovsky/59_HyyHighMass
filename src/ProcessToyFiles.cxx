////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotCLScan.cxx                                                      //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 01/04/2016                                                          //
//                                                                            //
//  Plots the resonant analysis limits as a function of MX.                   //
//                                                                            //
//  Macro options:                                                            //
//  - "New" or "FromFile"                                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "CLScan.h"
#include "Config.h"

TString m_inputDir;
TString m_outputDir;

std::vector<int> m_massValues;
std::vector<int> m_widthValues;
std::vector<int> m_xsValues;

/**
   -----------------------------------------------------------------------------
   Check if a value belongs to a vector.
   @param theVector - The vector that might contain the value.
   @param theValue - The value to check for membership in the vector.
   @return - True iff the vector already contains the value.
*/
bool vectorContainsValue(std::vector<int> theVector, int theValue) {
  for (int i_v = 0; i_v < (int)theVector.size(); i_v++) {
    if (theVector[i_v] == theValue) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Look at the available toy MC files to detect which mass, width, and cross-
   section values have associated pseudo-experiment ensembles.
   @param toyDirectory - The directory in which the toy MC files reside.
*/
void detectMassWidthXSFiles(TString toyDirectory) {

  // First clear existing list:
  m_massValues.clear();
  m_widthValues.clear();
  m_xsValues.clear();
  
  // Make a list of toy MC files and loop over the list:
  system(Form("ls %s | tee tempToyList.txt", toyDirectory.Data()));
  std::ifstream toyFileList("tempToyList.txt");
  TString currToyName;
  while (toyFileList >> currToyName) {
    // Only look at toys made for scans:
    if (currToyName.Contains("ForScan")) {
      
      // Split up the file name into components:
      TObjArray *array = currToyName.Tokenize("_");
      for (int i_t = 0; i_t < array->GetEntries(); i_t++) {
	TString currElement = ((TObjString*)array->At(i_t))->GetString();
	// Token for mass value:
	if (currElement.Contains("mass")) {
	  TString massStr = currElement;
	  massStr.ReplaceAll("mass", "");
	  int massVal = massStr.Atoi();
	  if (!vectorContainsValue(m_massValues, massVal)) {
	    m_massValues.push_back(massVal);
	  }
	}
	// Token for width value:
	else if (currElement.Contains("width")) {
	  TString widthStr = currElement;
	  widthStr.ReplaceAll("width", "");
	  int widthVal = widthStr.Atoi();
	  if (!vectorContainsValue(m_widthValues, widthVal)) {
	    m_widthValues.push_back(widthVal);
	  }
	}
	// Token for cross-section value:
	else if (currElement.Contains("xs")) {
	  TString xsStr = currElement;
	  xsStr.ReplaceAll("xs", "");
	  xsStr.ReplaceAll(".root", "");
	  int xsVal = xsStr.Atoi();
	  if (!vectorContainsValue(m_xsValues, xsVal)) {
	    m_xsValues.push_back(xsVal);
	  }
	}
      }
    }
  }
  // Close and remove list of files:
  toyFileList.close();
  system("rm tempToyList.txt");
  
  // Then sort the vectors of unique mass, width, and cross-section values:
  std::sort(m_massValues.begin(), m_massValues.end());
  std::sort(m_widthValues.begin(), m_widthValues.end());
  std::sort(m_xsValues.begin(), m_xsValues.end());
}

/**
   -----------------------------------------------------------------------------
   Hadd several ROOT files together.
*/
void haddCurrentFiles() {
  std::cout << "ProcessToyFiles: hadd a batch of files..." << std::endl;
  // Detect cross-sections, masses, widths:
  //detectMassWidthXSFiles(m_inputDir);
  detectMassWidthXSFiles("Graviton_Mar30/GenericToys/single_files/");
  
  // Loop over cross-section, mass, width:
  for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
    for (int i_w = 0; i_w < (int)m_widthValues.size(); i_w++) {
      for (int i_x = 0; i_x < (int)m_xsValues.size(); i_x++) {
	for (int mu = 0; mu < 2; mu++) {
	  TString outputFile = Form("%s/toy_ALL_mu%d_ForScan_mass%d_width%d_xs%d.root", m_outputDir.Data(), mu, m_massValues[i_m], m_widthValues[i_w], m_xsValues[i_x]);
	  TString inputFile = Form("Graviton_Mar30/GenericToys/single_files/toy_mu%d_*_ForScan_mass%d_width%d_xs%d.root", mu, m_massValues[i_m], m_widthValues[i_w], m_xsValues[i_x]);
	  // This three-step is done to prevent the output from being deleted
	  // but also to allow previous outputs to be added together.
	  system(Form("hadd -dummy.root %s %s", inputFile.Data()));
	  system(Form("rm %s", inputFile.Data()));
	  system(Form("mv dummy.root %s", outputFile.Data()));

	}
      }
    }
  }
}

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param inputDir - The directory containing .tar files.
   @param outputDir - The directory where the output .root files belong.
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 3) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <options>"
	      << std::endl;
    exit(0);
  }
  
  // Input directory: location of downloaded files
  m_inputDir = argv[1];
  // Output directory: location for ROOT files to go
  m_outputDir = argv[2];
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  
  // Create a list of input tar files:
  int maxIndex = 0;
  system(Form("ls %s | tee tarFileList.txt", m_inputDir.Data()));
  std::ifstream tarFileList("tarFileList.txt");
  if (tarFileList.is_open()) {

    // Loop over .tar files:
    TString currTarFile = "";
    while (tarFileList >> currTarFile) {
      std::cout << "Tar file " << maxIndex << " = " << currTarFile << std::endl;
      
      // Extract the files:
      system(Form("tar zxvf %s/%s", m_inputDir.Data(), currTarFile.Data()));
     
      // Continue extraction if fewer than 50 files have been extracted:
      if (maxIndex < 50) maxIndex++;
      
      // Otherwise hadd the existing files together:
      else {
	haddCurrentFiles();
	maxIndex = 0;
      }
      //system(Form("rm %s/%s", inputDir.Data(), currTarFile.Data()));
    }
    
    // Then hadd the remaining files:
    haddCurrentFiles();
    tarFileList.close();
  }
  
  // Delete the list of tar files:
  system("rm tarFileList.txt");
  
  return 0;
}
