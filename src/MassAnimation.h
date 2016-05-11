////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  MassAnimation.h                                                           //
//  Class: MassAnimation.cxx                                                  //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/04/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef MassAnimation_h
#define MassAnimation_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "TestStat.h"
#include "time.h"
#include "HGammaMxAOD.h"

class MassAnimation {

 public:
  
  MassAnimation(TString configFileName, TString options);
  virtual ~MassAnimation() {};
  
  void getDataForFrames();
  void makeGIF();
  void makeSingleFrame(int frame);
  void setNFrames(int nFrames);
  void setOutputDirectory(TString directory);
  
 private:
  
  std::vector<TString> makeLocalMxAODCopies(std::vector<TString> fileNames);
  TGraphErrors* plotSubtraction(RooAbsData *data, RooAbsPdf *pdf, 
				RooRealVar *observable, double xBins);
  void printer(TString statement, bool isFatal);
  void printProgressBar(int index, int total);
  void removeLocalMxAODCopies(std::vector<TString> fileNames);
  TString timeToString(UInt_t timeInt);
  
  // Private member variables:
  TString m_options;
  TString m_outputDir;
  TString m_configFileName;
  Config *m_config;
  
  // Workspace stuff:
  TFile *m_workspaceFile;
  RooWorkspace *m_workspace;
  ModelConfig *m_model;
  RooSimultaneous* m_combPdf;
  RooArgSet *m_observables;
  RooCategory *m_categories;

  int m_nFrames;
  std::map<int,TString> m_times;// use EventInfoAux.timeStamp
  std::map<int,RooDataSet*> m_data;
    
};

#endif

