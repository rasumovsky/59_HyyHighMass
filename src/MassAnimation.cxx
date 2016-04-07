////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: MassAnimation.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/04/2016                                                          //
//                                                                            //
//  Creates an animated diphoton mass specturm .gif.                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "MassAnimation.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the MassAnimation class.
   @param configFileName - The name of the analysis config file.
   @param options - Options for the CL scan.
*/
MassAnimation::MassAnimation(TString configFileName, TString options) {
    
  // Load the config file:
  m_configFileName = configFileName;
  m_config = new Config(m_configFileName);
  m_options = options;

  printer(Form("MassAnimation::MassAnimation(%s, %s)", configFileName.Data(),
	       options.Data()), false);

  // Set output directory:
  setInputDirectory(Form("%s/%s/GenericToys/single_files",
			 (m_config->getStr("MasterOutput")).Data(),
			 (m_config->getStr("JobName")).Data()));
  setOutputDirectory(Form("%s/%s/MassAnimation",
			  (m_config->getStr("MasterOutput")).Data(),
			  (m_config->getStr("JobName")).Data()));
  
  // Clear data:
  clearData();
  
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
}

/**
   -----------------------------------------------------------------------------
   Clear the class data.
*/
void MassAnimation::clearData() {
  
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void MassAnimation::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Set the output file location.
   @param directory - The directory for storing output files.
*/
void MassAnimation::setOutputDirectory(TString directory) {
  m_outputDir = directory;
  // Create output directory if it doesn't already exist:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
}
