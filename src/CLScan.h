////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  CLScan.h                                                                  //
//  Class: CLScan.cxx                                                         //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/03/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef CLScan_h
#define CLScan_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "TestStat.h"
#include "ToyAnalysis.h"

class CLScan {

 public:
  
  CLScan(TString configFileName, TString options);
  virtual ~CLScan() {};
  void clearData();
  double getLimit(int mass, int width, bool expected, int N);
  std::vector<int> listMasses();
  std::vector<int> listWidths();
  std::vector<int> listXS();
  void scanMass(int width, bool makeNew);
  void setInputDirectory(TString directory);
  void setLimit(int mass, int width, bool expected, int N, double limitValue);
  void setOutputDirectory(TString directory);
  bool singleCLScan(int mass, int width, bool makeNew);
  
 private:
  
  void detectMassWidthXSFiles(TString toyDirectory);
  double getIntercept(TGraph *graph, double valueToIntercept);
  void printer(TString statement, bool isFatal);
  bool vectorContainsValue(std::vector<int> theVector, int theValue);
  
  // Private member variables:
  TString m_options;
  TString m_inputDir;
  TString m_outputDir;
  TString m_configFileName;
  Config *m_config;
  bool m_filesLoaded;
  
  // Store mass, width, xs values for which toys have been generated:
  std::vector<int> m_massValues;
  std::vector<int> m_widthValues;
  std::vector<int> m_xsValues;
  
  // Store the limit values at a particular point:
  std::map<TString,double> m_values95CL;
  
};

#endif

