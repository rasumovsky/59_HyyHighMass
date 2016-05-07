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
#include "StatScan.h"
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
  
  // Define the scan settings:
  bool makeNew = config->getBool("MakeNewScan");
  bool useAsymptotics = config->getBool("UseAsymptoticsForScan");
  
  // Use the StatScan class to load toys and calculate 95% CL limits:
  StatScan *scan = new StatScan(configFile, 
				config->getStr("StatScanOptions"));
  if (options.Contains("New")) scan->setInputDirectory(inputDir);
  else scan->setInputDirectory(outputDir);
  scan->setOutputDirectory(outputDir);
  
  // Get scan points from config file for asymptotics:
  if (useAsymptotics) {
    scan->useTheseMasses(config->getIntV("AsympStatScanMasses"));
    scan->useTheseWidths(config->getIntV("AsympStatScanWidths"));
    //scan->useTheseXS(config->getIntV("AsympStatScanXS"));
  }
  
  // Loop over widths:
  std::vector<int> widths = scan->listWidths();
  for (int i_w = 0; i_w < (int)widths.size(); i_w++) {
    std::cout << "PlotStatScan: Start scan for width " << widths[i_w]
	      << std::endl;
    
    if (options.Contains("ScanLimit")) {
      scan->scanMassLimit(widths[i_w], makeNew, useAsymptotics);
    }
    else if (options.Contains("ScanP0")) {
      scan->scanMassP0(widths[i_w], makeNew, useAsymptotics);
    }
    else {
      std::cout << "PlotStatScan: ERROR! Option must = ScanLimit or ScanP0"
		<< std::endl;
    }
    std::cout << "PlotStatScan: Finished scan for width " << widths[i_w]
	      << std::endl;
  }
  std::cout << "PlotStatScan: Finished looping over " << widths.size() 
	    << " widths." << std::endl;
  return 0;
}
