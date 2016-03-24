////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  ToyAnalysis.h                                                             //
//  Class: ToyAnalysis.cxx                                                    //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/02/2016                                                          //
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
  
  ToyAnalysis(TString newConfigFile, TString options);
  virtual ~ToyAnalysis() {};
  
  bool areInputFilesOK();
  double calculateCLsFromToy(double qMu);
  double calculateCLFromToy(double qMu);
  double calculateErrorPVal(double pValue, int nToys);
  double calculateErrorFromCounting(double pValue, int nToys);
  double calculateErrorCLVal(double qMu);
  double calculateBkgQMuForN(double N);
  double calculatePBFromToy(double qMu);
  double calculatePMuFromToy(double qMu);
  bool doThisFit(TString fitType);
  double getPbFromN(double N);

  void fillToyHistograms(int muValue, ToyTree *toyTree);
  TH1F* getAsymptoticHist(TString statistic);
  std::vector<TString> getFitTypes(); 
  TH1F* getHist(TString paramName, TString fitType, int toyMu);
  TH1F* getMuHist(int toyMu);
  std::vector<TString> getNamesGlobalObservables();
  std::vector<TString> getNamesNuisanceParameters();
  std::vector<TString> getNamesPoI();
  TH1F* getStatHist(TString statistic, int toyMu);
  std::vector<double> getStatValues(TString statistic, int toyMu);
  void loadToy(int toyMu, TString toyFileForm);
  void plotHist(TString paramName, int toyMu);
  void plotProfiledMu(); 
  void plotRetries(int muValue);
  void plotTestStat(TString statistic);
  void plotTestStatComparison(TString statistic);
  void setFitTypes(std::vector<TString> fitTypes);
  void setOutputDir(TString outputDirectory);
  void setStatHistRanges(int nBins, int binMin, int binMax);
  void sortPairedVectors(std::vector<double> &vec1, std::vector<double> &vec2);
  
 private:

  void printer(TString statement, bool isFatal);
  TString printStatName(TString statistic);
  
  // Private member variables:
  TString m_options;
  TString m_outputDir;
  Config *m_config;
  bool m_filesLoaded;
  int m_resonanceMass;

  // Classes for statistics access:
  TestStat *m_ts;
  TFile *m_workspaceFile;
  RooWorkspace *m_workspace;
  ModelConfig *m_mc;
  
  // Test statistic binning:
  int m_nBins;
  int m_binMin;
  int m_binMax;
  
  // Fit types:
  std::vector<TString> m_fitTypes;
  
  // Storage of QMu for pMu calculation:
  std::vector<double> m_weightsIS_Mu0;
  std::vector<double> m_weightsIS_Mu1;
  std::vector<double> m_valuesQMu_Mu0;
  std::vector<double> m_valuesQMu_Mu1;
  std::vector<double> m_valuesMuHat_Mu0;
  std::vector<double> m_valuesMuHat_Mu1;

  // Storage of best-fit CL and p0 for calculation:
  std::vector<double> m_valuesBestFit_AsymCL_Mu0;
  std::vector<double> m_valuesBestFit_AsymCL_Mu1;
  std::vector<double> m_valuesBestFit_AsymZ0_Mu0;
  std::vector<double> m_valuesBestFit_AsymZ0_Mu1;
  
  // Histograms:
  TH1F *m_hAsymptoticQ0;
  TH1F *m_hAsymptoticQMu;
  TH1F *m_hMuProfiled[2];
  TH1F *m_hQ0[2];
  TH1F *m_hQMu[2];
  //TH1F *m_hQMuTilde[2];
  TH1F *m_hZ0[2];
  TH1F *m_hCL[2];

  // For investigations:
  TH1F *m_hRetries[2];
  TH1F *m_hImprovement[2];
  TH1F *m_hMedImprovement[2];
  TH1F *m_hCounter[2];
  TH2F *m_h2RetriesZ[2];
  TH2F *m_h2ZImprovement[2];
  TH1F *m_hImpAtThisStep[2];

  TH1F *m_hZ0Retries[2][50];

  // For the nuis, globs, and regular parameters:
  std::map<TString,TH1F*> m_histStorage;
  
  // Parameter data:
  std::vector<TString> m_namesGlobs;
  std::vector<TString> m_namesNuis;
  std::vector<TString> m_namesPoIs;
};

#endif

