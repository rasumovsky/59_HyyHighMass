////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: TestStat.h                                                          //
//  Class: TestStat.cxx                                                       //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 26/02/2016                                                          //
//                                                                            //
//  This class is used for statistics calculations.                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef TestStat_h
#define TestStat_h

// Package libraries:
#include "CommonHead.h"
#include "CommonFunc.h"
#include "Config.h"
#include "RooFitHead.h"
#include "statistics.h"

class TestStat {
  
 public:
  
  TestStat(TString newConfigFile, TString newOptions,
	   RooWorkspace *newWorkspace);
  virtual ~TestStat() {};
  
  void addGhostEvents(RooAbsData *dataset, RooRealVar *observable, 
		      RooRealVar *weight);
  std::vector<double> asymptoticCL(std::map<TString,double> mapPoI, 
				   TString datasetName, TString snapshotName,
				   TString poiForNorm, bool doTilde=false);
  //std::vector<double> asymptoticLimit(std::map<TString,double> mapPoI, 
  //				      TString datasetName, TString snapshotName,
  //				      TString poiForNorm);
  std::map<TString,double> asymptoticLimit(std::map<TString,double> mapPoI, 
					   TString datasetName, 
					   TString snapshotName,
					   TString poiForNorm, 
					   bool doTilde=false);
  //void asymptoticLimitBisector(std::map<TString,double> mapPoI, 
  //			       TString datasetName, TString poiForNorm, 
  //			       double nllMuHat, double profiledNorm,
  //			       double& muLimit, double& qMuLimit);
  std::vector<double> asymptoticP0(std::map<TString,double> mapPoI, 
				   TString datasetName, TString snapshotName,
				   TString poiForNorm);
  RooAbsData *binDataSet(TString binnedDataName, RooAbsData *unbinnedDataSet, 
			 RooRealVar* observable, RooRealVar* weight);
  void clearData();
  void clearFitParamSettings();
  RooAbsData* createAsimovData(int valPoI, TString snapshotName,
			       std::map<TString,double> namesAndValsPoI,
			       TString asimovDataName="");
  RooDataSet* createPseudoData(int seed, int valPoI, TString snapshotName,
			       std::map<TString,double> namesAndValsPoI, 
			       int toyIndex = -1);
  //double findMuUp(double muUpMed, double qMuAsimov, int N, double alpha=0.05);
  bool fitsAllConverged();
  double functionQ0(double x);
  double functionQMu(double x);
  double functionQMuTilde(double x, double asimovTestStat);
  double getCLFromCLs(double CLs);
  double getCLsFromCL(double CL);
  //double getCLFromQMu(double qMu, double N);
  //double getCLsFromQMu(double qMu, double N);
  double getCLFromQMu(double qMu, double sigma, double mu);
  double getCLsFromQMu(double qMu, double sigma, double mu);
  double getFitNLL(TString datasetName, int valPoI, bool fixPoI,
		   std::map<TString,double> namesAndValsPoI,
		   bool resetParams = true);
  std::map<std::string,double> getGlobalObservables();
  std::vector<double> getNEventsToys();
  std::map<std::string,double> getNuisanceParameters();
  std::map<std::string,double> getPoIs();
  double getP0FromQ0(double q0);
  double getP0FromR0(double r0);
  double getPFromN(double N);
  double getPbFromQMu(double qMu, double sigma, double mu);
  double getPbFromQMuTilde(double qMuTilde, double sigma, double mu);
  double getPMuFromQMu(double qMu);
  double getPMuFromQMuTilde(double qMuTilde, double sigma, double mu);
  double getPMuFromRMu(double rMu);
  double getQ0FromNLL(double nllMu0, double nllMuHat, double muHat);
  double getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
		       double muTest);
  double getQMuTildeFromNLL(double nllMu, double nllMu0, double nllMuHat,
			    double muHat, double muTest);
  double getR0FromNLL(double nllMu0, double nllMuHat, double muHat);
  double getRMuFromNLL(double nllMu, double nllMuHat, double muHat, 
		       double muTest);
  double getRMuTildeFromNLL(double nllMu, double nllMu0, double nllMuHat,
			    double muHat, double muTest);
  //double getSigma(double qMu, double mu, double muHat, bool doTilde, 
  //		  int direction);
  double getSigma(double qMu, double mu, double muHat, bool doTilde=false);
  double getZ0FromQ0(double q0);
  double getZFromP(double p);
  double getZMuFromQMu(double qMu);
  double graphIntercept(TGraph *graph, double valueToIntercept);
  TGraph *nllScanGraph();
  void resetParamsAfterFit(bool doResetParamsAfterFit);
  void saveSnapshots(bool doSaveSnapshot);
  void scaleAsimovData(double scaleFactor, std::vector<TString> varsToScale);
  std::map<int,double> scanNLL(TString scanName, TString datasetName,
			       TString varToScan, bool fixPoI,
			       std::map<TString,double> namesAndValsPoI,
			       int nScanPoints);
  void setFitOptions(TString fitOptions);
  void setNominalSnapshot(TString nominalSnapshot);
  void setPlotAxis(bool useLogScale, double yMin, double yMax, 
		   double GeVPerBin);
  void setPlotDirectory(TString directory);
  void setParam(TString paramName, double paramVal, bool doSetConstant);
  void setSubPlot(TString ratioOrSubtraction);
  RooWorkspace *theWorkspace();
  ModelConfig *theModelConfig();
  void useTwoSidedTestStat(bool useTwoSided);
  
 private:

  TString nameOfVar(TString varForm);
  TGraphErrors* plotDivision(TString dataName, TString pdfName, TString obsName,
			     double xMin, double xMax, double xBins);
  void plotFits(TString fitType, TString datasetName);
  TGraphErrors* plotSubtract(TString dataName, TString pdfName,
			     TString obsName, double xMin, double xMax,
			     double xBins);
  void printer(TString statement, bool isFatal);
  void printSet(TString setName, RooArgSet* set);
  void storeParams(RooArgSet *set, std::map<std::string,double>& map);
  double upperLimitFinder(std::map<TString,double> mapPoI, TString datasetName, 
			  TString poiForNorm, bool doTilde=false);
  
  // From the initialization:
  TString m_anaType;    // The analysis type ("Res", "NonRes").
  TString m_jobName;    // The name of the group of jobs (for I/O purposes).
  TString m_options;    // Job options.
  TString m_outputDir;  // The output directory for statistics.
  
  TString m_dataForObsQ0;  // The dataset for expected results.
  TString m_dataForObsQMu; // The dataset for observed results.
  TString m_dataForExpQ0;  // The dataset for expected results.
  TString m_dataForExpQMu; // The dataset for observed results.
  
  TString m_plotDir;    // The output directory for plots (not set by default).
  bool m_doSaveSnapshot;// Option to save snapshots of PoI and nuisance params.
  bool m_doPlot;        // Sets whether or not to plot fit results.

  // Option to use two-sided test-statistics:
  bool m_useTwoSided;

  // Pointer to the input file:
  TFile *inputFile;
  
  // The configuration of the analysis:
  Config *m_config;

  // Check whether all fits successful:
  bool m_allGoodFits;
  bool m_doResetParamsAfterFit;
  TString m_fitOptions; 
  
  // The workspace for the fits:
  RooWorkspace *m_workspace;
  ModelConfig *m_mc;
  TGraph *m_graphNLL;
  TString m_nominalSnapshot;
  
  // Store fit parameters from NLL calculation:
  std::map<std::string,double> m_mapGlobs;
  std::map<std::string,double> m_mapNP;
  std::map<std::string,double> m_mapPoI;
  
  // Store number of toy MC events generated:
  std::vector<double> m_numEventsPerCate;
  
  // In case special parameter settings are used for a fit:
  std::map<TString,double> m_paramValToSet;
  std::map<TString,bool> m_paramConstToSet;
  
  // Asimov data settings:
  double m_AsimovScaleFactor;
  std::vector<TString> m_AsimovVarsToScale;
  
  // Plot settings:
  bool m_useLogScale;
  double m_yMin;
  double m_yMax;
  double m_geVPerBin;
  TString m_ratioOrSubtraction;
};

#endif

