////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: PlotStatScan.cxx                                                    //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 11/04/2016                                                          //
//                                                                            //
//  Plots p0 values or CLs limits (observed, expected) as a function of mass. //
//                                                                            //
//  Macro options:                                                            //
//  - "New" or "FromFile"                                                     //
//  - "95CL" or "p0"                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "CLScan.h"
#include "Config.h"

/**
   -----------------------------------------------------------------------------
   The main method plots the p0 or 95CL values as a function of mass.
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
  TString inputDir = Form("%s/%s/GenericToys/single_files",
			  (config->getStr("MasterOutput")).Data(),
			  (config->getStr("JobName")).Data());
  TString outputDir = Form("%s/%s/PlotStatScan", 
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  // Use the CLScan class to load toys and calculate 95% CL limits:
  CLScan *scan = new CLScan(configFile, config->getStr("PlotStatScanOptions"));
  scan->setInputDirectory(inputDir);
  scan->setOutputDirectory(outputDir);
  std::vector<int> widths = scan->listWidths();
  for (int i_w = 0; i_w < (int)widths.size(); i_w++) {
    if (options.Contains("ScanLimit")) {
      scan->scanMassLimit(widths[i_w], options.Contains("New"));
    }
    else if (options.Contains("ScanP0")) {
      scan->scanMassP0(widths[i_w], options.Contains("New"));
    }
    else {
      std::cout << "PlotStatScan: ERROR! Option must specify p0 or 95CL"
		<< std::endl;
    }
  }
  
  return 0;
}
