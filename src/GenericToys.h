////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: GenericToys.h                                                       //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 22/03/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "TestStat.h"
#include <time.h>
#include "RooFitHead.h"
#include "RooStatsHead.h"
#include "statistics.h"

Config *m_config;

TString m_configFile;
TString m_options;
int m_nToysPerJob;
int m_inputPoIVal;
  
// Name of the copied workspace file:
TString m_copiedFile;

// Output directory:
TString m_outputDir;

// Output file and TTree:
TFile *m_outputFile;
TTree *m_outputTree;

// Indices for tracking files and inputs:
int m_toyMassInt;
int m_toyWidthInt;
int m_toyXSectionInt;

// Variables to store in the TTree:
int m_seed;
double m_toyMass;
double m_toyWidth;
double m_toyXSection;
int m_bestFitUpdate;
double m_weight;
double m_numEvents;
bool m_convergedMu1;
bool m_convergedMu0;
bool m_convergedMuFree;
double m_profiledPOIVal;
double m_nllMu0;
double m_nllMu1;
double m_nllMuFree;
double m_llrL1L0;
double m_llrL0Lfree;
double m_llrL1Lfree;
std::vector<std::string> m_namesNP;
std::vector<double> m_valuesNPMu0;
std::vector<double> m_valuesNPMu1;
std::vector<double> m_valuesNPMuFree;
std::vector<string> m_namesGlobs;
std::vector<double> m_valuesGlobsMu0;
std::vector<double> m_valuesGlobsMu1;
std::vector<double> m_valuesGlobsMuFree;
std::vector<string> m_namesPoIs;
std::vector<double> m_valuesPoIsMu0;
std::vector<double> m_valuesPoIsMu1;
std::vector<double> m_valuesPoIsMuFree;
std::vector<double> m_numEventsPerCate;
std::vector<double> m_nllPerRetry;
