////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  DupSleuth.cxx                                                             //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 25/07/2016                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This main method provides a tool for performing individual fits to the    //
//  resonance Monte Carlo. Settings for the utility are provided in           //
//  singleSigFitExample.cfg.                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/// c++ headers
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <list>
#include <map>
#include <cstdlib>

// ROOT headers:
#include "TChain.h"
#include "TFile.h"

std::map<UInt_t,std::vector<ULong64_t> > r2e;

/**
   -----------------------------------------------------------------------------
   Check for duplciate run+event combinations. Use a map to reduce the size of
   the problem (2nd order loops over individual runs instead of all runs). 
   @param run - The run number of the event.
   @param event - The event number.
   @return - True iff. event is duplicated. 
*/
bool isDuplicate(UInt_t run, ULong64_t event) {
  
  // Check to see if map already exists for the run:
  if (r2e.count(run) > 0) {
    
    // Loop through events in the given run:
    for (int i_v = 0; i_v < (int)(r2e[run]).size(); i_v++) {
      if (event == r2e[run][i_v]) return true;
    }
    
    // Not duplicated:
    r2e[run].push_back(event);
    return false;    
  }
  
  // Map entry does not exist for this run. Create a new map.
  else {
    r2e[run].clear();
    r2e[run].push_back(event);
    return false;
  }
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void PrintProgressBar(int index, int total) {
  if (index%10000 == 0) {
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
  // Check that the config file location is provided.
  if (argc < 2) {
    std::cout << "ERROR! Improper number of args." << std::endl;
    exit(0);
  }
  
  // Create TChain of input files:
  TChain *chain = new TChain("CollectionTree");
  for (int i_a = 1; i_a < argc; i_a++) {
    TString currFileName = argv[i_a];
    std::cout << "\t Adding file " << currFileName << std::endl;
    chain->AddFile(currFileName);
  }
  
  // Track number of duplicate events and store list of run,event combinations:
  std::ofstream duplicateList("duplicateList.txt");
  int numDuplicates = 0;
  r2e.clear();
  
  // Assign the MxAOD/TTree branches to variables:
  UInt_t runNumber;
  ULong64_t eventNumber;
  chain->SetMakeClass(1);
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("EventInfoAux.eventNumber", 1);
  chain->SetBranchStatus("EventInfoAux.runNumber", 1);
  chain->SetBranchAddress("EventInfoAux.runNumber", &runNumber);
  chain->SetBranchAddress("EventInfoAux.eventNumber", &eventNumber);
    
  //--------------------------------------//
  // Loop over events to check for duplicates:
  int nEvents = chain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." 
	    << std::endl;
  for (int index = 0; index < nEvents; index++) {
    chain->GetEntry(index);
    PrintProgressBar(index, nEvents);
    if (isDuplicate(runNumber, eventNumber)) {
      numDuplicates++;
      duplicateList << "run=" << runNumber << " \tevent=" << eventNumber
		    << std::endl;
    }
  }
  
  //--------------------------------------//
  // Save and print summary statement:
  duplicateList.close();
  std::cout << "There were " << numDuplicates << " duplicated events.\n"
	    << "\tFor details: duplicateList.txt." << std::endl;
  delete chain;
  return 0;
}
