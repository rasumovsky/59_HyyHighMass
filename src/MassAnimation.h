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
  int getNFrames();
  void makeGIF();
  void makeAllFrames();
  void makeSingleFrame(int frame);
  void setNFrames(int nFrames);
  void setOutputDirectory(TString directory);
  
 private:
  
  std::vector<TString> makeLocalMxAODCopies(std::vector<TString> fileNames);
  TGraphErrors* plotComparison(RooAbsData *data, RooAbsPdf *pdf, 
			       RooRealVar *observable, double xBins,
			       bool doRatio=true);
  void printer(TString statement, bool isFatal);
  void printProgressBar(int index, int total);
  void removeLocalMxAODCopies(std::vector<TString> fileNames);
  TString timeToString(time_t dateValue);
  
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
  
  TString m_obsName = "";
  TString m_pdfName = "";
  int m_geVPerBin;
  int m_nFrames;
  std::map<int,TString> m_times;// use EventInfoAux.timeStamp
  std::map<int,RooDataSet*> m_data;
    
};

#endif

