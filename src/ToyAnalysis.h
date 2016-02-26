////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  ToyAnalysis.h                                                             //
//  Class: ToyAnalysis.cxx                                                    //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 26/02/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ToyAnalysis_h
#define ToyAnalysis_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "ToyTree.h"
#include "TestStat.h"
#include "RooFitHead.h"
#include "statistics.h"

class ToyAnalysis {

 public:
  
  ToyAnalysis(TString newConfigFile, TString options, int resonanceMass);
  virtual ~ToyAnalysis() {};
  
  bool areInputFilesOK();
  double calculateCLsFromToy(double qMu);
  double calculateCLFromToy(double qMu);
  double calculateErrorPVal(double pValue, int nToys);
  double calculateErrorCLVal(double qMu);
  double calculateBkgQMuForN(double N);
  double calculatePBFromToy(double qMu);
  double calculatePMuFromToy(double qMu);
  double getPbFromN(double N);

  void fillToyHistograms(int muValue, ToyTree *toyTree);
  void getAsymptoticForm(TString statistic);
  TH1F* getAsymptoticHist();
  TH1F* getHist(TString paramName, TString fitType, int toyMu);
  TH1F* getMuHist(int toyMu);
  TH1F* getStatHist(TString statistic, int toyMu);
  void plotHist(TString paramName, int toyMu);
  void plotProfiledMu(); 
  void plotTestStat(TString statistic);
  void plotTestStatComparison(TString statistic);
    
 private:

  void printer(TString statement, bool isFatal);
  TString printStatName(TString statistic);
  
  // Private member variables:
  TString m_outputDir;
  Config *m_config;
  bool m_filesLoaded;
  int m_resonanceMass;

  // Classes for statistics access:
  TestStat *m_ts;
  RooWorkspace *m_workspace;
  
  // Test statistic binning:
  int m_nBins;
  int m_binMin;
  int m_binMax;
  
  // Fit types:
  std::vector<TString> m_fitTypes;
  
  // Storage of QMu for pMu calculation:
  std::vector<double> m_valuesQMu_Mu0;
  std::vector<double> m_valuesQMu_Mu1;

  // Histograms:
  TH1F *m_hAsymptotic;
  TH1F *m_hMuProfiled[2];
  TH1F *m_hQ0[2];
  TH1F *m_hQMu[2];
  //TH1F *m_hQMuTilde[2];
  
  // For the nuis, globs, and regular parameters:
  std::map<TString,TH1F*> m_histStorage;
  
  // Parameter data:
  std::vector<TString> m_namesGlobs;
  std::vector<TString> m_namesNuis;
  std::vector<TString> m_namesPars;
};

#endif

