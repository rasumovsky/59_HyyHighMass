////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  StatScan.h                                                                //
//  Class: StatScan.cxx                                                       //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/03/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef StatScan_h
#define StatScan_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "TestStat.h"
#include "ToyAnalysis.h"

class StatScan {

 public:
  
  StatScan(TString configFileName, TString options);
  virtual ~StatScan() {};
  void clearData();
  double getCL(int mass, int width, int crossSection, bool expected,
	       bool asymptotic, int N);
  double getLimit(int mass, int width, bool expected, bool asymptotic, int N);
  double getP0(int mass, int width, bool expected, bool asymptotic);
  std::vector<int> listMasses();
  std::vector<int> listWidths();
  std::vector<int> listXS();
  void scanMassLimit(int width, bool makeNew, bool asymptotic);
  void scanMassP0(int width, bool makeNew, bool asymptotic);
  void setInputDirectory(TString directory);
  void setCL(int mass, int width, int crossSection, bool expected, 
	     bool asymptotic, int N, double CLValue);
  void setLimit(int mass, int width, bool expected, bool asymptotic, int N, 
		double limitValue);
  void setP0(int mass, int width, bool expected, bool asymptotic, 
	     double p0Value);
  void setOutputDirectory(TString directory);
  bool singleCLScan(int mass, int width, bool makeNew, bool asymptotic);
  bool singleCLTest(int mass, int width, int crossSection, bool makeNew,
		    bool asymptotic);
  bool singleP0Test(int mass, int width, int crossSection, bool makeNew, 
		    bool asymptotic);
  void useTheseMasses(std::vector<int> massValues);
  void useTheseWidths(std::vector<int> widthValues);
  void useTheseXS(std::vector<int> xsValues);
  
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
  
  // Store the statistical values at a particular point:
  std::map<TString,double> m_valuesCL;
  std::map<TString,double> m_valuesLimit;
  std::map<TString,double> m_valuesP0;
  
};

#endif

