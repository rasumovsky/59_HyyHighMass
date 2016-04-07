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

class MassAnimation {

 public:
  
  MassAnimation(TString configFileName, TString options);
  virtual ~MassAnimation() {};
  void clearData();
  void setOutputDirectory(TString directory);
  
 private:
  
  void printer(TString statement, bool isFatal);
    
  // Private member variables:
  TString m_options;
  TString m_outputDir;
  TString m_configFileName;
  Config *m_config;
    
};

#endif

