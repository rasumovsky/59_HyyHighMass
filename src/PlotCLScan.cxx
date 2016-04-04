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

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options. Can be "toy" or "asymptotic" or "both"
   @param resMass - The resonance mass.
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
  TString inputDir = Form("%s/%s/GenericToys/single_files",
			  (config->getStr("MasterOutput")).Data(),
			  (config->getStr("JobName")).Data());
  TString outputDir = Form("%s/%s/PlotCLScan", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Use the CLScan class to load toys and calculate 95% CL limits:
  CLScan *scan = new CLScan(configFile, config->getStr("PlotCLScanOptions"));
  scan->setInputDirectory(inputDir);
  scan->setOutputDirectory(outputDir);
  std::vector<int> widths = scan->listWidths();
  for (int i_w = 0; i_w < (int)widths.size(); i_w++) {
    scan->scanMass(widths[i_w], options.Contains("New"));
  }
  
  return 0;
}
