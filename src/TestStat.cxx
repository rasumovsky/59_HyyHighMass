////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: TestStat.cxx                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 26/02/2016                                                          //
//                                                                            //
//  This class is used for statistics calculations.                           //
//                                                                            //
//  Note, May 25: Disabling expected p0 fits, which double the run time.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "TestStat.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the TestStat class. 
   @param newConfigFile - The analysis config file.
   @param newOptions - The job options ("New", "FromFile"), etc.
   @param newWorkspace - The workspace with the model for the test stats. 
*/
TestStat::TestStat(TString newConfigFile, TString newOptions,
		   RooWorkspace *newWorkspace) {
  std::cout << "TestStat: Initializing...\n\t" << newConfigFile << "\n\t"
	    << newOptions << "\n\t" << std::endl;
  
  // Assign input variables:
  m_options = newOptions;
  
  // Start with a clean class:
  clearData();
  
  // Set the analysis configuration:
  m_config = new Config(newConfigFile);
  m_jobName = m_config->getStr("JobName");
  m_anaType = m_config->getStr("AnalysisType");
  
  // Use Asimov data if the analysis is blind.
  m_dataForExpQ0 = m_config->getStr("WorkspaceAsimovDataMu1");
  m_dataForExpQMu = m_config->getStr("WorkspaceAsimovDataMu0");
  m_dataForObsQ0 = m_config->getStr("WorkspaceObsData");
  m_dataForObsQMu = m_config->getStr("WorkspaceObsData");
  
  // Load from file if the pointer passed is NULL:
  if (newWorkspace == NULL) {
    inputFile = new TFile(m_config->getStr("WorkspaceFile"), "read");
    if (inputFile->IsOpen()) {
      printer("TestStat: Loading workspace.", false);
      m_workspace
	= (RooWorkspace*)inputFile->Get(m_config->getStr("WorkspaceName"));
    }
    else {
      printer(Form("TestStat: ERROR loading file %s.", 
		   (m_config->getStr("WorkspaceFile")).Data()), true);
    }
  }
  // Use the workspace passed to the class constructor:
  else m_workspace = newWorkspace;
  
  // Then get the model configuration:
  if (m_workspace->obj(m_config->getStr("WorkspaceModelConfig"))) {
    m_mc = (ModelConfig*)m_workspace
      ->obj(m_config->getStr("WorkspaceModelConfig"));
  }
  else {
    printer(Form("TestStat: ERROR, missing model config %s",
		 (m_config->getStr("WorkspaceModelConfig")).Data()), true);
  }
  
  // Reset nuisance parameters after fitting by default:
  resetParamsAfterFit(true);
    
  // Use linear plot y-axis by default:
  setPlotAxis(false, 0.0, 100.0, 1.0);
  setSubPlot("Subtraction");
  
  // Create output directories:
  m_outputDir = Form("%s/%s/TestStat", 
		     (m_config->getStr("MasterOutput")).Data(),
		     m_jobName.Data());
  system(Form("mkdir -vp %s", m_outputDir.Data()));
  system(Form("mkdir -vp %s/CL/", m_outputDir.Data()));
  system(Form("mkdir -vp %s/p0/", m_outputDir.Data()));
  
  // Finished instantiating the test statistic class. 
  std::cout << "TestStat: Initialized Successfully!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Add ghost events to a dataset.
   @param dataset - A pointer to the dataset w/o ghosts.
   @param observable - A pointer to the RooRealVar observable in the dataset.
   @param weight - A pointer to the RooRealVar weight variable in the dataset.
*/
void TestStat::addGhostEvents(RooAbsData *dataset, RooRealVar* observable,
			      RooRealVar* weight) {
  for(double obsVal = observable->getMin(); obsVal < observable->getMax(); 
      obsVal += (observable->getMax() - observable->getMin())/1000.0) {
    double obsWeight = 0.0000001;
    observable->setVal(obsVal);
    weight->setVal(obsWeight);
    dataset->add(RooArgSet(*observable, *weight), obsWeight);
  }
}

/**
   -----------------------------------------------------------------------------
   Calculate the asymptotic CL values using fits.
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset in the workspace.
   @param snapshotName - The name of the ML snapshot for generating mu=0 Asimov.
   @param poiForNorm - The name of the normalization poi.
   @param doTilde - True iff. using qMuTilde instead of qMu.
   @return - A vector with obs, exp+2s, exp+1s, exp, exp-1s, exp-2s CL values.
*/
std::vector<double> TestStat::asymptoticCL(std::map<TString,double> mapPoI, 
					   TString datasetName, 
					   TString snapshotName,
					   TString poiForNorm, bool doTilde) {
  printer(Form("TestStat::asymptoticCL(%s, %s)",
	       datasetName.Data(), poiForNorm.Data()), false);
  
  printer("REQUIRES REVISION, DEPRECATED!", true);

  // Create unique Asimov dataset name based on PoI:
  TString asimovDataMu0Name = "asimovDataMu0";

  // Need to set the non-norm PoI constant...
  for (std::map<TString,double>::iterator poiIter = mapPoI.begin(); 
       poiIter != mapPoI.end(); poiIter++) {
    TString poiName = poiIter->first;
    double poiValue = poiIter->second;
    if (!poiName.EqualTo(poiForNorm)) setParam(poiName, poiValue, true);
    asimovDataMu0Name.Append(Form("_%s%2.2f", poiName.Data(), poiValue));
  }
  
  // Generate mu0 asimov data, 
  if (!m_workspace->data(asimovDataMu0Name)) {
    double originNorm = mapPoI[poiForNorm];
    mapPoI[poiForNorm] = 0.0;//temporarily set the signal strength to zero.
    createAsimovData(0.0, snapshotName, mapPoI, asimovDataMu0Name);
    mapPoI[poiForNorm] = originNorm;//restore the non-zero signal strength.
  }
  
  // For  now, use mu = 1 (qmu for mu=1), since we are adjusting the cross-
  // section parameter and not the mu parameter.
  double muForQMu = 1.0;
  
  // Calculate observed qMu: 
  double nllMu1Obs = getFitNLL(datasetName, 1.0, true, mapPoI, false);
  double nllMuHatObs = getFitNLL(datasetName, 1.0, false, mapPoI, false);
  double profiledNormObs = (getPoIs())[(std::string)(poiForNorm)];
  double muHatObs = (profiledNormObs / mapPoI[poiForNorm]);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, muForQMu);
    
  // Calculate expected qMu:
  double nllMu1Exp = getFitNLL(asimovDataMu0Name, 1.0, true, mapPoI, false);
  double nllMuHatExp = getFitNLL(asimovDataMu0Name, 0.0, false, mapPoI, false);
  double profiledNormExp = (getPoIs())[(std::string)(poiForNorm)];
  double muHatExp = (profiledNormExp / mapPoI[poiForNorm]);
  double expQMu = getQMuFromNLL(nllMu1Exp, nllMuHatExp, muHatExp, muForQMu);
    
  // Replace qMu -> qMuTilde if requested (nllMu0 requires mu=0 fit):
  if (doTilde) {
    // Temporarily set normalization parameter to zero:
    double originNorm = mapPoI[poiForNorm];
    mapPoI[poiForNorm] = 0.0;
    
    // Calculate observed qMuTilde:
    double nllMu0Obs = getFitNLL(datasetName, 0.0, true, mapPoI, false);
    obsQMu = getQMuTildeFromNLL(nllMu1Obs, nllMu0Obs, nllMuHatObs, muHatObs,
				muForQMu);
    
    // Calculate expected qMuTilde:
    double nllMu0Exp = getFitNLL(asimovDataMu0Name, 0.0, true, mapPoI, false);
    expQMu = getQMuTildeFromNLL(nllMu1Exp, nllMu0Exp, nllMuHatExp, muHatExp,
				muForQMu);
    // Reset normalization parameter:
    mapPoI[poiForNorm] = originNorm;
  }

  // Get Sigma:
  //int direction = (N < 0) ? -1 : 1;
  double sigmaExp = getSigma(expQMu, muForQMu, muHatExp, doTilde);
  double sigmaObs = getSigma(obsQMu, muForQMu, muHatObs, doTilde);
  
  /*
    
  // Calculate nllMuHat only once:
  double nllMuHat = getFitNLL(datasetName, 0.0, false, mapPoI, false);
  double profiledNorm = (getPoIs())[(std::string)(poiForNorm)];
      
    // Set the signal strength:
    mapPoI[poiForNorm] = muLimit;
    
    // Calculate qMu:
    double nllMu1 = getFitNLL(datasetName, 1.0, true, mapPoI, false);
    double muHat = (profiledNorm / mapPoI[poiForNorm]);
    double qMu = getQMuFromNLL(nllMu1, nllMuHat, muHat, muForQMu);
    
    // Replace qMu with qMuTilde if requested (requires extra fit for nllMu0):
    if (doTilde) {
      mapPoI[poiForNorm] = 0.0;
      double nllMu0 = getFitNLL(datasetName, 0.0, true, mapPoI, false);
      mapPoI[poiForNorm] = muLimit;
      qMu = getQMuTildeFromNLL(nllMu1, nllMu0, nllMuHat, muHat, muForQMu);
    }
    
    // Get Sigma:
    //int direction = (N < 0) ? -1 : 1;
    double sigma = getSigma(qMu, muForQMu, muHat, doTilde);
    
    // Which Mu to use? Probably muHat, right?
    CLs = getCLsFromQMu(qMu, sigma, muForQMu);
  */
  
  // Calculate CL:
  double expCL = getCLsFromQMu(expQMu, sigmaExp, muForQMu);
  double obsCL = getCLsFromQMu(obsQMu, sigmaObs, muForQMu);
  
  // Print summary:
  std::cout << "\t expected CL nom = " << expCL << std::endl;
  std::cout << "\t observed CL = " << obsCL << "\n" << std::endl;
  
  // Return observed and expected CL values:
  std::vector<double> resultsCL;
  resultsCL.clear();
  resultsCL.push_back(obsCL);
  resultsCL.push_back(expCL);
  return resultsCL;
}

/*
   -----------------------------------------------------------------------------
   Calculate the asymptotic limit values using fits.
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset in the workspace.
   @param snapshotName - The name of the ML snapshot for generating mu=0 Asimov.
   @param poiForNorm - The name of the normalization poi.
   @return - A vector with obs, exp+2s, exp+1s, exp, exp-1s, exp-2s CL values.

std::vector<double> TestStat::asymptoticLimit(std::map<TString,double> mapPoI, 
					      TString datasetName, 
					      TString snapshotName,
					      TString poiForNorm) {
  printer(Form("TestStat::asymptoticLimit(%s, %s)",
	       datasetName.Data(), poiForNorm.Data()), false);
  
  // Create unique Asimov dataset name based on PoI:
  TString asimovDataMu0Name = "asimovDataMu0";
  
  // Need to set the non-norm PoI constant...
  for (std::map<TString,double>::iterator poiIter = mapPoI.begin(); 
       poiIter != mapPoI.end(); poiIter++) {
    TString poiName = poiIter->first;
    double poiValue = poiIter->second;
    if (!poiName.EqualTo(poiForNorm)) setParam(poiName, poiValue, true);
    asimovDataMu0Name.Append(Form("_%s%2.2f", poiName.Data(), poiValue));
  }
  
  // Generate mu0 asimov data, 
  if (!m_workspace->data(asimovDataMu0Name)) {
    double originNorm = mapPoI[poiForNorm];
    mapPoI[poiForNorm] = 0.0;//temporarily set the signal strength to zero.
    createAsimovData(0.0, snapshotName, mapPoI, asimovDataMu0Name);
    mapPoI[poiForNorm] = originNorm;//restore the non-zero signal strength.
  }
  
  // Use bisection method to find median observed 95% CL limit on mu:
  double observedLimit = 0.0; double qMuLimitObs = 0.0;
  double nllMuHatObs = getFitNLL(datasetName, 0.0, false, mapPoI, false);
  double profiledNormObs = (getPoIs())[(std::string)(poiForNorm)];
  asymptoticLimitBisector(mapPoI, datasetName, poiForNorm, nllMuHatObs,
			  profiledNormObs, observedLimit, qMuLimitObs);
  
  // Use bisection method to find median expected 95% CL limit on mu:
  double expectedLimit = 0.0; double qMuLimitExp = 0.0;
  double nllMuHatExp = getFitNLL(asimovDataMu0Name, 0.0, false, mapPoI, false);
  double profiledNormExp = (getPoIs())[(std::string)(poiForNorm)];
  asymptoticLimitBisector(mapPoI, asimovDataMu0Name, poiForNorm, nllMuHatExp,
			  profiledNormExp, expectedLimit, qMuLimitExp);
  
  // Then calculate expected error bands:
  double expectedLimit_n2 = findMuUp(expectedLimit, qMuLimitExp, -2);
  double expectedLimit_n1 = findMuUp(expectedLimit, qMuLimitExp, -1);
  double expectedLimit_p1 = findMuUp(expectedLimit, qMuLimitExp, 1);
  double expectedLimit_p2 = findMuUp(expectedLimit, qMuLimitExp, 2);
  
  // Print summary:
  std::cout << "\n\t observed 95% CL limit = " << observedLimit << "\n"
	    << "\t expected 95% CL limit -2s = " << expectedLimit_n2 << "\n"
	    << "\t expected 95% CL limit -1s = " << expectedLimit_n1 << "\n"
	    << "\t expected 95% CL limit nom = " << expectedLimit << "\n"
	    << "\t expected 95% CL limit +1s = " << expectedLimit_p1 << "\n"
	    << "\t expected 95% CL limit +2s = " << expectedLimit_p2 << "\n"
	    << std::endl;
  
  // Return observed and expected CL values:
  std::vector<double> resultsLimit;
  resultsLimit.clear();
  resultsLimit.push_back(observedLimit);
  resultsLimit.push_back(expectedLimit_n2);
  resultsLimit.push_back(expectedLimit_n1);
  resultsLimit.push_back(expectedLimit);
  resultsLimit.push_back(expectedLimit_p1);
  resultsLimit.push_back(expectedLimit_p2);
  return resultsLimit;
}
*/

/**
   -----------------------------------------------------------------------------
   Calculate the asymptotic limit values using fits.
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset in the workspace.
   @param snapshotName - The name of the ML snapshot for generating mu=0 Asimov.
   @param poiForNorm - The name of the normalization poi.
   @param doTilde - True iff. using qMuTilde instead of qMu.
   @return - A vector with obs, exp+2s, exp+1s, exp, exp-1s, exp-2s CL values.
*/
std::map<TString,double> TestStat::asymptoticLimit(
   std::map<TString,double> mapPoI, TString datasetName, TString snapshotName,
   TString poiForNorm, bool doTilde) {
  
  printer(Form("TestStat::asymptoticLimit(%s, %s)",
	       datasetName.Data(), poiForNorm.Data()), false);
  
  // Map to store limit results:
  std::map<TString,double> limitMap; limitMap.clear();
  
  // Create unique Asimov dataset name based on PoI:
  TString asimovDataMu0Name = "asimovDataMu0";
  
  // Need to set the non-norm PoI constant...
  for (std::map<TString,double>::iterator poiIter = mapPoI.begin(); 
       poiIter != mapPoI.end(); poiIter++) {
    TString poiName = poiIter->first;
    double poiValue = poiIter->second;
    if (!poiName.EqualTo(poiForNorm)) setParam(poiName, poiValue, true);
    asimovDataMu0Name.Append(Form("_%s%2.2f", poiName.Data(), poiValue));
  }
  
  // Generate mu0 asimov data:
  if (!m_workspace->data(asimovDataMu0Name)) {
    double originNorm = mapPoI[poiForNorm];
    mapPoI[poiForNorm] = 0.0;//temporarily set the signal strength to zero.
    createAsimovData(0.0, snapshotName, mapPoI, asimovDataMu0Name);
    mapPoI[poiForNorm] = originNorm;//restore the non-zero signal strength.
  }
  
  // Do NLL Scan to get N=-2,-1,0,+1,+2 values of mu:
  TString nllScanName = "scan_";
  for (std::map<TString,double>::iterator poiIter = mapPoI.begin();
       poiIter != mapPoI.end(); poiIter++) {
    nllScanName.Append(Form("%s_%f",(poiIter->first).Data(),poiIter->second));
  }
  std::map<int,double> mapNToNorm 
    = scanNLL(nllScanName, asimovDataMu0Name, poiForNorm, false, mapPoI, 10);
  
  // Set normalization variable range:
  double normMinOrigin = m_workspace->var(poiForNorm)->getMin();
  double normMaxOrigin = m_workspace->var(poiForNorm)->getMax();
  m_workspace->var(poiForNorm)
    ->setMin(mapNToNorm[0] - (5.0 * fabs(mapNToNorm[-1] - mapNToNorm[0])));
  m_workspace->var(poiForNorm)
    ->setMax(mapNToNorm[0] + (5.0 * fabs(mapNToNorm[1] - mapNToNorm[0])));
  
  //----------------------------------------//
  // Loop over values of N to define expected CLs exclusion limit and bands:
  for (int N = -2; N <= 2; N++) {
    
    // Create Asimov data:
    TString asimovCurr = asimovDataMu0Name;
    if (N != 0) {
      asimovCurr = Form("asimovDataBand%d",N);
      double originNorm = mapPoI[poiForNorm];
      mapPoI[poiForNorm] = mapNToNorm[N];
      createAsimovData(0.0, snapshotName, mapPoI, asimovCurr);
      mapPoI[poiForNorm] = originNorm;//restore the non-zero signal strength.
    }
    
    // Find mu that gives CLs = 0.05 via iteration:
    limitMap[Form("ExpN%d",N)]
      = upperLimitFinder(mapPoI, asimovCurr, poiForNorm, doTilde);
  }
  
  // Find the observed CLs exclusion limit:
  m_workspace->var(poiForNorm)
    ->setMin(mapNToNorm[0] - (10.0 * fabs(mapNToNorm[-1] - mapNToNorm[0])));
  m_workspace->var(poiForNorm)
    ->setMax(mapNToNorm[0] + (10.0 * fabs(mapNToNorm[1] - mapNToNorm[0])));
  limitMap["Obs"] = upperLimitFinder(mapPoI, datasetName, poiForNorm, doTilde);
  
  // Print summary:
  std::cout << "\n\t observed 95% CL limit = " << limitMap["Obs"] << "\n"
	    << "\t expected 95% CL limit -2s = " << limitMap["ExpN-2"] << "\n"
	    << "\t expected 95% CL limit -1s = " << limitMap["ExpN-1"] << "\n"
	    << "\t expected 95% CL limit nom = " << limitMap["ExpN0"] << "\n"
	    << "\t expected 95% CL limit +1s = " << limitMap["ExpN1"] << "\n"
	    << "\t expected 95% CL limit +2s = " << limitMap["ExpN2"] << "\n"
	    << std::endl;
  
  // Reset the normalization variable range:
  m_workspace->var(poiForNorm)->setMin(normMinOrigin);
  m_workspace->var(poiForNorm)->setMax(normMaxOrigin);
  
  // Return observed and expected CL values:
  return limitMap;
}

/*
   -----------------------------------------------------------------------------
   Bisection method to search for 95% limit intercept. Basically adjust mu for 
   the numerator in qMu (lMu/lfree) until sqrt(qMu)=1.64. muLimit and qMuLimit
   are passed by reference and modified to give the 95% CL limit values.
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset to fit.
   @param poiForNorm - The name of the normalization poi.
   @param nllMuHat - The nll value from the maximum likelihood fit.
   @param profiledNorm - The signal strength from the ML fit.
   @param muLimit - The variable to store the 95% CL limit on mu.
   @param qMuLimit - The variable to store the 95% CL limit on qMu.

void TestStat::asymptoticLimitBisector(std::map<TString,double> mapPoI, 
				       TString datasetName, TString poiForNorm,
				       double nllMuHat, double profiledNorm,
				       double& muLimit, double& qMuLimit) {
  printer(Form("TestStat::asymptoticLimitBisector(%s, %s)",
	       datasetName.Data(), poiForNorm.Data()), false);
  
  double precision = 0.0001;
  int nIterations = 0;
  int maxIterations = 30;
  double stepSize = (m_workspace->var(poiForNorm)->getMax() - 
		     m_workspace->var(poiForNorm)->getMin()) / 2.0;
  muLimit = (m_workspace->var(poiForNorm)->getMax() + 
	     m_workspace->var(poiForNorm)->getMin()) / 2.0;
  qMuLimit = 0.0;
  double qMuIntercept = 1.64; // 1-phi(1.64) = 0.05 -> 95% CL.
  while ((fabs(sqrt(qMuLimit) - qMuIntercept) / qMuIntercept) > precision && 
	 nIterations <= maxIterations) {
    std::cout << "TestStat:asymptoticLimitBisector: Iteration " << nIterations 
	      << " starting, muLimit=" << muLimit << ", qMuLimit=" << qMuLimit
	      << std::endl;
    
    // Set the signal strength:
    mapPoI[poiForNorm] = muLimit;
    
    // Calculate expected qmu:
    double nllMu1 = getFitNLL(datasetName, 1.0, true, mapPoI, false);
    double muHat = (profiledNorm / mapPoI[poiForNorm]);
    qMuLimit = getQMuFromNLL(nllMu1, nllMuHat, muHat, 1);
    
    // Increase the iteration count and reduce the step size:
    nIterations++;
    stepSize = 0.5 * stepSize;
    
    // Update the value of the mu limit:
    if (sqrt(qMuLimit) > qMuIntercept) muLimit -= stepSize;
    else muLimit += stepSize;
  } // At this point, muLimit should be the 95% CL limit on signal strength.
  
  if (nIterations == maxIterations) {
    printer("TestStat::asymptoticLimitBisector: Maximum Iterations!", true);
  }
  
}
*/

/**
   -----------------------------------------------------------------------------
   Calculate the p0 value using model fits.
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset in the workspace.
   @param snapshotName - The name of the ML snapshot for generating mu=1 Asimov.
   @param poiForNorm - The name of the normalization poi.
   @return - A vector with obs, exp.
*/
std::vector<double> TestStat::asymptoticP0(std::map<TString,double> mapPoI, 
					   TString datasetName,
					   TString snapshotName,
					   TString poiForNorm) {
  printer(Form("TestStat::asymptoticP0(%s, %s)", 
	       datasetName.Data(), poiForNorm.Data()), false);
  
  // Create unique Asimov dataset name based on PoI:
  TString asimovDataMu1Name = "asimovDataMu1";

  // Need to set the non-norm PoI constant...
  for (std::map<TString,double>::iterator poiIter = mapPoI.begin(); 
       poiIter != mapPoI.end(); poiIter++) {
    TString poiName = poiIter->first;
    double poiValue = poiIter->second;
    if (!poiName.EqualTo(poiForNorm)) setParam(poiName, poiValue, true);
    asimovDataMu1Name.Append(Form("_%s%2.2f", poiName.Data(), poiValue));
  }
  
  // Generate mu1 asimov data for expected:
  //if (!m_workspace->data(asimovDataMu1Name)) {
  //createAsimovData(1.0, snapshotName, mapPoI, asimovDataMu1Name);
  //}
  
  // Then set the signal strength to zero for fitting:
  double originNorm = mapPoI[poiForNorm];
  mapPoI[poiForNorm] = 0.0;
  
  // Calculate observed q0: 
  double nllMu0Obs = getFitNLL(datasetName, 0.0, true, mapPoI, false);
  double nllMuHatObs = getFitNLL(datasetName, 0.0, false, mapPoI, false);
  double profiledNormObs = (getPoIs())[(std::string)(poiForNorm)];
  double muHatObs = (profiledNormObs / originNorm);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  double obsR0 = getR0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  /*
  double nllMu0Exp = getFitNLL(asimovDataMu1Name, 0.0, true, mapPoI, false);
  double nllMuHatExp = getFitNLL(asimovDataMu1Name, 0.0, false, mapPoI, false);
  double profiledNormExp = (getPoIs())[(std::string)(poiForNorm)];
  double muHatExp = (profiledNormExp / originNorm);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  double expR0 = getR0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  */
  // Calculate p0 from q0:
  double expP0 = 0.5;//m_useTwoSided ? getP0FromR0(expR0) : getP0FromQ0(expQ0);
  double obsP0 = m_useTwoSided ? getP0FromR0(obsR0) : getP0FromQ0(obsQ0);
  
  // Print summary:
  std::cout << "\n\t Expected p0 = " << expP0 << std::endl;
  std::cout << "\t Observed p0 = " << obsP0 << "\n" << std::endl;
  
  // Return observed and expected p0 values:
  std::vector<double> resultsP0;
  resultsP0.clear();
  resultsP0.push_back(obsP0);
  resultsP0.push_back(expP0);
  return resultsP0;
}

/**
   -----------------------------------------------------------------------------
   Create a binned dataset from an unbinned dataset. 
   @param binnedDataName - A name for the new binned dataset.
   @param unbinnedDataSet - A pointer to the unbinned dataset.
   @param observable - A pointer to the RooRealVar observable in the dataset.
   @param weight - A pointer to the RooRealVar weight variable in the dataset.
   @return - A binned version of the unbinned dataset that was provided. 
*/
RooAbsData *TestStat::binDataSet(TString binnedDataName, 
				 RooAbsData *unbinnedDataSet, 
				 RooRealVar* observable, RooRealVar* weight) {
  
  // Convert dataset into TH1F:
  TH1F *hDataTemp = (TH1F*)unbinnedDataSet
    ->createHistogram("hDataTemp", *observable,
		      RooFit::Binning(observable->getBins(),
				      observable->getMin(),
				      observable->getMax()));
  
  // Create a dataset to fill with binned entries:
  RooDataSet* binnedData = new RooDataSet(binnedDataName, binnedDataName,
					  RooArgSet(*observable, *weight),
					  RooFit::WeightVar(*weight));
  for (int i_b = 1; i_b <= hDataTemp->GetNbinsX(); i_b++) {
    observable->setVal(hDataTemp->GetXaxis()->GetBinCenter(i_b));
    weight->setVal(hDataTemp->GetBinContent(i_b));
    binnedData->add(RooArgSet(*observable, *weight),
		    hDataTemp->GetBinContent(i_b));
  }
  
  return binnedData;
}

/**
   -----------------------------------------------------------------------------
   Clears all data stored by the class, but does not modify the workspace.
*/
void TestStat::clearData() {
  m_allGoodFits = true;
  m_mapGlobs.clear();
  m_mapNP.clear();
  m_mapPoI.clear();
  m_doSaveSnapshot = false;
  m_doPlot = false;
  m_plotDir = "";
  m_graphNLL = NULL;
  m_fitOptions = ""; 
  m_nominalSnapshot = "paramsOrigin";
  m_AsimovScaleFactor = 1.0;
  m_AsimovVarsToScale.clear();
  
  clearFitParamSettings();
  useTwoSidedTestStat(false);
}

/**
   -----------------------------------------------------------------------------
   Clears all specifications for parameter values during fits.
*/
void TestStat::clearFitParamSettings() {
  m_paramValToSet.clear();
  m_paramConstToSet.clear();
}

/**
   -----------------------------------------------------------------------------
   Create an Asimov dataset. WARNING! Method currently assumes inclusive
   categorization scheme in order to get category normalization.
   @param valPoI - The value of the parameter of interest.
   @param snapshotName - The name of the parameter snapshot to load.
   @param namesAndValsPoI - The names and values of the parameters of interest.
   @param asimovDataName - Specify a name for asimov data if default is no good.
   @return - An Asimov dataset.
*/
RooAbsData* TestStat::createAsimovData(int valPoI, TString snapshotName,
				       std::map<TString,double> namesAndValsPoI,
				       TString asimovDataName) {
  
  printer(Form("TestStat::createAsimovData(PoI=%d, snapshot=%s)", 
	       valPoI, snapshotName.Data()), false);
  
  // Load the parameters from profiling data:
  if (m_workspace->getSnapshot(snapshotName)) {
    m_workspace->loadSnapshot(snapshotName);
    printer(Form("TestStat: Loaded snapshot %s", snapshotName.Data()), false);
  }
  else {
    printer(Form("TestStat: ERROR! No snapshot %s", snapshotName.Data()), true);
  }
  
  // First check that PDF is of simultaneous type:
  if (!((TString)(m_mc->GetPdf()->ClassName())).Contains("Simultaneous")) {
    printer("TestStat::createPseudoData() ERROR. PDF not RooSimultaneous",true);
  }
    
  // Set the values of parameters of interest as specified in the input.
  for (std::map<TString,double>::iterator iterPoI = namesAndValsPoI.begin();
       iterPoI != namesAndValsPoI.end(); iterPoI++) {
    if (m_workspace->var(iterPoI->first)) {
      m_workspace->var(iterPoI->first)->setVal(iterPoI->second);
      m_workspace->var(iterPoI->first)->setConstant(true);
    }
    else {
      printer(Form("TestStat: ERROR! %s parameter missing",
		   (iterPoI->first).Data()), true);
    }
  }
  
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("TestStat: Error! Parameter %s not found in workspace! See printout above for clues...", (iterParam->first).Data()), true);
    }
  }

  // Scale Asimov data if set by user:
  for (int i_v = 0; i_v < (int)m_AsimovVarsToScale.size(); i_v++) {
    double currValue = m_workspace->var(m_AsimovVarsToScale[i_v])->getVal();
    m_workspace->var(m_AsimovVarsToScale[i_v])
      ->setVal(currValue * m_AsimovScaleFactor);
  }
  
  // Create and return the Asimov dataset:
  RooAbsData* asimovData
    = RooStats::AsymptoticCalculator::GenerateAsimovData(*m_mc->GetPdf(), 
						       *m_mc->GetObservables());
  if (asimovDataName.EqualTo("")) {
    asimovData->SetNameTitle(Form("asimovDataMu%d",valPoI),
			     Form("asimovDataMu%d",valPoI));
  }
  else asimovData->SetNameTitle(asimovDataName, asimovDataName);
  m_workspace->import(*asimovData);
  return asimovData;
}

/**
   -----------------------------------------------------------------------------
   Create a pseudo-dataset.
   @param seed - The random seed for dataset generation.
   @param valPoI - The value of the parameter of interest.
   @param snapshotName - The name of the parameter snapshot to load.
   @param namesAndValsPoI - The names and values of the parameters of interest.
   @param toyIndex - In case multiple toys are to be saved in this workspace.
   @return - A pseudo-dataset.
*/
RooDataSet* TestStat::createPseudoData(int seed, int valPoI, 
				       TString snapshotName,
				       std::map<TString,double> namesAndValsPoI,
				       int toyIndex) {
  printer(Form("TestStat::createPseudoData(seed=%d, PoI=%d, snapshot=%s)", 
	       seed, valPoI, snapshotName.Data()), false);
  
  // Load the parameters from profiling data to create UNCONDITIONAL ENSEMBLE:
  //TString snapshotName = Form("paramsProfilePoI%d", valPoI);
  if (m_workspace->getSnapshot(snapshotName)) {
    m_workspace->loadSnapshot(snapshotName);
    printer(Form("TestStat: Loaded snapshot %s",snapshotName.Data()), false);
  }
  else {
    // Load the original parameters from profiling:
    //m_workspace->loadSnapshot(m_nominalSnapshot);
    printer(Form("TestStat: ERROR! No snapshot %s", snapshotName.Data()),
	    true);
  }
  
  // First check that PDF is of simultaneous type:
  if (!((TString)(m_mc->GetPdf()->ClassName())).Contains("Simultaneous")) {
    printer("TestStat::createPseudoData() ERROR. PDF not RooSimultaneous",true);
  }
  
  RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  
  RooRandom::randomGenerator()->SetSeed(seed);
  statistics::constSet(nuisanceParameters, true);
  statistics::constSet(globalObservables, false);
  
  // Randomize the global observables and set them constant for now:
  statistics::randomizeSet(combPdf, globalObservables, seed); 
  statistics::constSet(globalObservables, true);
    
  // Set the values of parameters of interest as specified in the input.
  for (std::map<TString,double>::iterator iterPoI = namesAndValsPoI.begin();
       iterPoI != namesAndValsPoI.end(); iterPoI++) {
    TString currParName = iterPoI->first;
    double currParVal = iterPoI->second;
    if (m_workspace->var(currParName)) {
      m_workspace->var(currParName)->setVal(currParVal);
      m_workspace->var(currParName)->setConstant(true);
    }
    else {
      printer(Form("TestStat: ERROR! %s parameter missing",
		   (iterPoI->first).Data()), true);
    }
  }
  
  // Check if other parameter settings have been specified for toys:
  // WARNING! This overrides the randomization settings above!
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("TestStat: Error! Parameter %s not found in workspace! See printout above for clues...", (iterParam->first).Data()), true);
    }
  }
  
  // Store the toy dataset and number of events per dataset:
  map<string,RooDataSet*> toyDataMap; toyDataMap.clear();
  m_numEventsPerCate.clear();
  RooRealVar *wt = new RooRealVar("wt","wt",1.0);
  
  // Iterate over the categories:
  printer("TestStat: Iterate over categories to generate toy data.", false);
  //RooSimultaneous *simPdf = (RooSimultaneous*)m_workspace->pdf("combinedPdf");
  //TIterator *cateIter = simPdf->indexCat().typeIterator();
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  
  while ((cateType = (RooCatType*)cateIter->Next())) {
    
    RooAbsPdf *currPdf = combPdf->getPdf(cateType->GetName());
    //RooAbsPdf *currPdf = simPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPdf->getObservables(observables);
    
    // If you want to bin the pseudo-data (speeds up calculation):
    if (m_config->getBool("DoBinnedFit")) {
      currPdf->setAttribute("PleaseGenerateBinned");
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs, AutoBinned(true),
					 Extended(currPdf->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
      std::cout << "TOYNEVENTS = " 
		<< toyDataMap[(std::string)cateType->GetName()]->numEntries()
		<< std::endl;
    }
    /*
    // If you want to manually bin the pseudo-data (speeds up calculation):
    if (m_config->getBool("DoBinnedFit")) {
    RooAbsData *TestStat::binDataSet(RooAbsData *unbinnedDataSet, 
				       RooRealVar* observable, 
				       RooRealVar* weight);
    */ 
    
    // Construct unbinned pseudo-data by default:
    else {
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs,Extended(true));
      std::cout << "TOYNEVENTS " << cateType->GetName() << " = " 
		<< toyDataMap[(std::string)cateType->GetName()]->numEntries()
		<< std::endl;
    }
    
    double currEvt = toyDataMap[(std::string)cateType->GetName()]->sumEntries();
    m_numEventsPerCate.push_back(currEvt);

    // Add ghost events if requested:
    //if (m_config->isDefined("AddGhostEvents") &&
    //	m_config->getBool("AddGhostEvents")) {
    //addGhostEvents(toyDataMap[(std::string)cateType->GetName()],
    //		     (RooRealVar*)(currObs->first()), wt);
    //currObs->add(*wt);
    //}
    
    /*
    TCanvas *can = new TCanvas("can", "can", 800, 800);
    can->cd();
    RooPlot* frame = ((RooRealVar*)(currObs->first()))->frame(115);
    toyDataMap[(std::string)cateType->GetName()]->plotOn(frame);
    frame->SetXTitle(currObs->GetName());
    frame->SetYTitle(Form("Events / %2.1f GeV", m_geVPerBin));
    frame->Draw();
    gPad->SetLogy(); 
    frame->GetYaxis()->SetRangeUser(0.1, 10000);
    can->Print("data.eps");
    */    
  }
  
  // Create the combined toy RooDataSet:
  //RooCategory *categories = (RooCategory*)m_workspace->obj("categories");
  RooCategory *categories = NULL;
  TString nameRooCategory = m_config->getStr("WorkspaceRooCategory");
  if (m_workspace->obj(nameRooCategory)) {
    categories = (RooCategory*)m_workspace->obj(nameRooCategory);
  }
  else {
    printer(Form("TestStat: RooCategory object %s not found",
		 nameRooCategory.Data()), true);
  }
  RooDataSet* pseudoData = NULL;
  TString newDataName = (toyIndex < 0) ? "toyData" : Form("toyData%d",toyIndex);
  //if (m_config->isDefined("AddGhostEvents") &&
  //  m_config->getBool("AddGhostEvents")) {
  //pseudoData = new RooDataSet(newDataName, newDataName, *observables, 
  //				RooFit::Index(*categories),
  //				RooFit::Import(toyDataMap),
  //				RooFit::WeightVar(*wt));
  //}
  //else {
    pseudoData = new RooDataSet(newDataName, newDataName, *observables, 
				RooFit::Index(*categories),
				RooFit::Import(toyDataMap));
    //}
  
  // Save the parameters used to generate toys:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams(poi, m_mapPoI);

  // release nuisance parameters (but don't change the values!):
  //m_workspace->loadSnapshot(m_nominalSnapshot);
  statistics::constSet(nuisanceParameters, false);
  statistics::constSet(globalObservables, true);
  
  // Import into the workspace then return:
  m_workspace->import(*pseudoData);
  return pseudoData;
}

/*
   -----------------------------------------------------------------------------
   Get the expected limit for error band (CLs).
   @param muUpMed - The median 95% CLs mu value.
   @param qMuAsimov - The value of qMu from fitting mu=0 Asimov data
   @param N - The Gaussian quintile.
   @param alpha = The probability (alpha=0.05 for 95% limits).
   @return - The 95% CLs limit.
double TestStat::findMuUp(double muUpMed, double qMuAsimov, int N,
			  double alpha) {
  double sigma = sqrt( (muUpMed * muUpMed) / qMuAsimov);
  double muUpN = sigma *
    (TMath::NormQuantile(1.0 - (alpha * ROOT::Math::gaussian_cdf(N))) + N);
  return muUpN;
}
*/

/**
   -----------------------------------------------------------------------------
   Check if all of the fits done by this class have converged.
   @return - True iff. all of the fits have been successfully convergent.
*/
bool TestStat::fitsAllConverged() { 
  return m_allGoodFits;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of q0 (same as qMu).
   @param x - The value of the test statistic.
   @return - The value of the asymptotic q0 test statistic distribution.
*/
double TestStat::functionQ0(double x) {
  // This corresponds to the "special case" of mu'=0
  double result = TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMu.
   @param x - The value of the test statistic.
   @return - The value of the qMu asymptotic test statistic distribution.
*/
double TestStat::functionQMu(double x) {
  // This corresponds to the "special case" of mu'=mu
  double result = TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x));
  return result;
}

/**
   -----------------------------------------------------------------------------
   Implements the functional form of qMuTilde.
   @param x - The value of the test statistic.
   @param asimovTestStat - The test stat value on Asimov data with mu=0 but
   fitting under mu=1 hypothesis.
   @return - The value of the asymptotic test statistic distribution.
*/
double TestStat::functionQMuTilde(double x, double asimovTestStat) {
  // This corresponds to the "special case" of mu=mu'
  double result = 0.0;
  double cutoff = asimovTestStat; // asimov test stat...
  if (x == 0) result = 0.0;
  else if (x <= cutoff) {
    result = (TMath::Exp(-1*x/2.0) / (2.0*sqrt(2.0*TMath::Pi()*x)));
  }
  else {
    result = (TMath::Exp(-1*TMath::Power((x+cutoff),2) / (8*cutoff))
	      / (2*sqrt(2*TMath::Pi()*cutoff)));
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Get the CL value from CLs.
   @param CLs - The CLs value to convert to CL.
   @return - The corresponding CL value.
*/
double TestStat::getCLFromCLs(double CLs) {
  return (1.0 - CLs);
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value from CL.
   @param CL - The CL value to convert to CLs.
   @return - The corresponding CLs value.
*/
double TestStat::getCLsFromCL(double CL) {
  return (1.0 - CL);
}

/**
   -----------------------------------------------------------------------------
   Get the CL value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.

double TestStat::getCLFromQMu(double qMu, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, N));
  return CL;
}
*/

/**
   -----------------------------------------------------------------------------
   Get the CL value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param sigma - The sigma value...
   @param mu - The mu value... 
   @return - The CL value.
*/
double TestStat::getCLFromQMu(double qMu, double sigma, double mu) {
  return getCLFromCLs(getCLsFromQMu(qMu, sigma, mu));
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.

double TestStat::getCLsFromQMu(double qMu, double N) {
  // N = 0 for exp and obs
  double pMu = getPMuFromQMu(qMu);
  double pB = getPFromN(N);
  double CLs = pMu / (1.0 - pB);
  return CLs;
}
*/

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param sigma - The sigma value...
   @param mu - The mu value... 
   @return - The CLs value.
*/
double TestStat::getCLsFromQMu(double qMu, double sigma, double mu) {
  // N = 0 for exp and obs
  double pMu = getPMuFromQMu(qMu);
  double pB = getPbFromQMu(qMu, sigma, mu);
  double CLs = pMu / (1.0 - pB);
  return CLs;
}

/**
   -----------------------------------------------------------------------------
   Get the negative-log-likelihood for a fit of a specified type to a specified
   dataset.
   @param datasetName - The name of the dataset in the workspace.
   @param valPoI - The parameter of interest value to fix (0 or 1).
   @param fixPoI - True if PoI should be fixed to the specified value.
   @param namesAndValsPoI - Map of names and values of PoIs to set for fit.
   @param resetParams - True iff original parameter values are used at start. 
   @return - The NLL value.
*/
double TestStat::getFitNLL(TString datasetName, int valPoI, bool fixPoI,
			   std::map<TString,double> namesAndValsPoI,
			   bool resetParams) {
  printer(Form("TestStat: getFitNLL(%s, PoI=%d, fixPoI=%d, resetPars=%d)",
	       datasetName.Data(),valPoI,(int)fixPoI,(int)resetParams),false);
  
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  if (resetParams) m_workspace->loadSnapshot(m_nominalSnapshot);
  RooArgSet* origValNP
    = (RooArgSet*)m_workspace->getSnapshot(m_nominalSnapshot);
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  // Look for dataset:
  if (!m_workspace->data(datasetName)) {
    printer(Form("TestStat: Error! Requested data not available: %s",
		 datasetName.Data()), true);
  }
  
  // release nuisance parameters before fit and set to the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // the global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  // Set the values of parameters of interest as specified in the input.
  for (std::map<TString,double>::iterator iterPoI = namesAndValsPoI.begin();
       iterPoI != namesAndValsPoI.end(); iterPoI++) {
    if (m_workspace->var(iterPoI->first)) {
      m_workspace->var(iterPoI->first)->setVal(iterPoI->second);
      m_workspace->var(iterPoI->first)->setConstant(fixPoI);
    }
    else {
      printer(Form("TestStat: ERROR! %s parameter missing",
		   (iterPoI->first).Data()), true);
    }
  }
  
  // Check if other parameter settings have been specified for fit:
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("TestStat: Error! Parameter %s not found in workspace!",
		   ((TString)iterParam->first).Data()), true);
    }
  }
  
  // The actual fit command:
  RooNLLVar* varNLL
    = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				     Extended(combPdf->canBeExtended()));
  
  RooFitResult *fitResult
    = statistics::minimize(varNLL, m_fitOptions, NULL, true);
  if (!fitResult || fitResult->status() != 0) m_allGoodFits = false;
  else m_allGoodFits = true;

  // Save a snapshot if requested:
  if (m_doSaveSnapshot) {
    TString textValPoI = fixPoI ? (Form("%d",(int)valPoI)) : "Free";
    m_workspace->saveSnapshot(Form("paramsProfilePoI%s", textValPoI.Data()),
			      *poiAndNuis);
  }
  
  // Plot the fit result if the user has set an output directory for plots:
  if (m_doPlot) {
    if (fixPoI && ((int)valPoI) == 1) plotFits("Mu1", datasetName);
    else if (fixPoI && ((int)valPoI) == 0) plotFits("Mu0", datasetName);
    else plotFits("MuFree", datasetName);
  }
  
  double nllValue = varNLL->getVal();
  delete varNLL;
  delete fitResult;
  
  // Save names and values of nuisance parameters, globs, other parameters:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams(poi, m_mapPoI);
  
  // Release nuisance parameters after fit and recover the default values:
  if (m_doResetParamsAfterFit) {
    statistics::constSet(nuisanceParameters, false, origValNP);
  }
  else {
    statistics::constSet(nuisanceParameters, false);
  }
  
  // Finish up, return NLL value.
  printer("TestStat: Fit has completed. Returning NLL.", false);
  return nllValue;
}

/**
   -----------------------------------------------------------------------------
   Get a map of global observable names to values from the most recent fit.
   @return - A map of global observable names and most recent fit values.
*/
std::map<std::string,double> TestStat::getGlobalObservables() {
  return m_mapGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get the number of events in each category in the most recent pseudo-dataset.
   @return - A vector containing the weighted # events per category.
*/
std::vector<double> TestStat::getNEventsToys() {
  return m_numEventsPerCate;
}

/**
   -----------------------------------------------------------------------------
   Get a map of nuisance parameter names to values from the most recent fit.
   @return - A map of nuisance parameter names and most recent fit values.
*/
std::map<std::string,double> TestStat::getNuisanceParameters() {
  return m_mapNP;
}

/**
   -----------------------------------------------------------------------------
   Get a map of parameter of interest names to values from most recent fit.
   @return - A map of PoI names and most recent fit values.
*/
std::map<std::string,double> TestStat::getPoIs() {
  return m_mapPoI;
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of p0 based on the test statistic q0.
   @param q0 - The test statistic q0.
   @return - The value of p0.
*/
double TestStat::getP0FromQ0(double q0) {
  return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(q0))));
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of p0 based on the test statistic r0.
   @param r0 - The test statistic r0.
   @return - The value of p0.
*/
double TestStat::getP0FromR0(double r0) {
  if (r0 > 0) return (1 - ROOT::Math::gaussian_cdf(sqrt(r0)));
  else return (1 - ROOT::Math::gaussian_cdf(-1 * sqrt(fabs(r0))));
}

/**
   -----------------------------------------------------------------------------
   Calculate p-value based on the standard deviation.
   @param N - The standard deviation.
   @return - The p-value.
*/
double TestStat::getPFromN(double N) {
  return (1 - ROOT::Math::gaussian_cdf(N));
}

/**
   -----------------------------------------------------------------------------
   Calculate the pB value based on qMu.
   @param qMu - The test statistic qMu.
   @param sigma - The sigma value...
   @param mu - The mu value... 
   @return - The value of pB.
*/
double TestStat::getPbFromQMu(double qMu, double sigma, double mu) {
  return (1 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMu)));
}

/**
   -----------------------------------------------------------------------------
   Calculate the pB value based on qMuTilde.
   @param qMuTilde - The test statistic qMuTilde.
   @param sigma - The sigma value...
   @param mu - The mu value... 
   @return - The value of pB.
*/
double TestStat::getPbFromQMuTilde(double qMuTilde, double sigma, double mu) {
  if (qMuTilde > 0 && qMuTilde <= (mu*mu)/(sigma*sigma)) {
    return (1.0 - ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qMuTilde)));
  }
  else {
    return 1.0 - ROOT::Math::gaussian_cdf((mu*mu/(sigma*sigma) - qMuTilde) / 
					  (2.0*fabs(mu/sigma)));
  }
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of pMu.
   @param qMu - The test statistic qMu.
   @return - The value of pMu.
*/
double TestStat::getPMuFromQMu(double qMu) {
  return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu))));
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of pMu.
   @param qMuTilde - The test statistic qMuTilde.
   @param sigma - The sigma value...
   @param mu - The mu value...
   @return - The value of pMu.
*/
double TestStat::getPMuFromQMuTilde(double qMuTilde, double sigma, double mu) {
  if (qMuTilde <= (mu*mu)/(sigma*sigma)) {
    return (1.0 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMuTilde))));
  }
  else { //if (qMuTilde > (mu*mu)/(sigma*sigma))
    return (1.0 - ROOT::Math::gaussian_cdf((qMuTilde + ((mu*mu)/(sigma*sigma)))/
					   (2.0 * fabs(mu/sigma))));
  }
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of rMu.
   @param rMu - The test statistic rMu.
   @return - The value of pMu.
*/
double TestStat::getPMuFromRMu(double rMu) {
  if (rMu > 0) return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(rMu))));
  else return (1 - ROOT::Math::gaussian_cdf(-1 * sqrt(fabs(rMu))));
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic q0 based on the nll.
   @param nllMu0 - NLL of a fit with signal strength 0;
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @return - The value of q0.
*/
double TestStat::getQ0FromNLL(double nllMu0, double nllMuHat, double muHat) {
  double q0 = (muHat < 0) ? 0 : (2 * (nllMu0 - nllMuHat));
  if (q0 < 0.0) q0 = 0.0;
  return q0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMu based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of qMu.
*/
double TestStat::getQMuFromNLL(double nllMu, double nllMuHat, double muHat,
			       double muTest) {
  if (muHat <= muTest) return (2.0 * (nllMu - nllMuHat));
  else return 0.0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic qMuTilde based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMu0 - NLL of a fit with signal strength 0.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of qMuTilde.
*/
double TestStat::getQMuTildeFromNLL(double nllMu, double nllMu0,
				    double nllMuHat, double muHat,
				    double muTest) {
  if (muHat < 0) return (2.0 * (nllMu - nllMu0));
  else if (muHat >= 0 && muHat <= muTest) return (2.0 * (nllMu - nllMuHat));
  else return 0.0; //(muHat > muTest)
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic r0 (two-sided q0) based on the nll.
   @param nllMu0 - NLL of a fit with signal strength 0;
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @return - The value of r0.
*/
double TestStat::getR0FromNLL(double nllMu0, double nllMuHat, double muHat) {
  if (muHat <= 0.0) return (-2 * (nllMu0 - nllMuHat));
  else return (2.0 * (nllMu0 - nllMuHat));
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic rMu (two-sided qMu) based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of rMu.
*/
double TestStat::getRMuFromNLL(double nllMu, double nllMuHat, double muHat,
			       double muTest) {
  if (muHat <= muTest) return (2.0 * (nllMu - nllMuHat));
  else return (-2.0 * (nllMu - nllMuHat));
}

/**
   -----------------------------------------------------------------------------
   Calculate the test statistic rMuTilde (two-sided qMu tilde) based on the nll.
   @param nllMu - NLL of a fit with signal strength mu.
   @param nllMu0 - NLL of a fit with signal strength 0.
   @param nllMuHat - NLL of a fit with profiled signal strength.
   @param muHat - Profiled signal strength.
   @param muTest - Tested value of signal strength.
   @return - The value of rMuTilde.
*/
double TestStat::getRMuTildeFromNLL(double nllMu, double nllMu0,
				    double nllMuHat, double muHat,
				    double muTest) {
  if (muHat <= 0) return (2.0 * (nllMu - nllMu0));
  else if (muHat > 0 && muHat <= muTest) return (2.0 * (nllMu - nllMuHat));
  else return (-2.0 * (nllMu - nllMuHat)); //(muHat > muTest)
}

/**
   -----------------------------------------------------------------------------
   Calculate the sigma parameter.
   @param qMu - The test statistic qMu (or qMu tilde if doTilde=true).
   @param mu - Tested value of signal strength.
   @param muHat - Profiled signal strength.
   @param doTilde - True iff using qMuTilde instead of qMu.
   @param direction =  -1 if N < 0, else +1.
   @return - The value of sigma.

double TestStat::getSigma(double qMu, double mu, double muHat, bool doTilde,
			  int direction) {
  if (mu*direction < muHat) return (fabs(mu-muHat) / sqrt(qMu));
  else if (muHat < 0 && doTilde) {
    return sqrt(mu*mu-2*mu*muHat*direction) / sqrt(qMu);
  }
  else return (mu-muHat)*direction / sqrt(qMu);
}
*/

/**
   -----------------------------------------------------------------------------
   Calculate the sigma parameter. Note: it is always positive...
   @param qMu - The test statistic qMu (or qMu tilde if doTilde=true).
   @param mu - Tested value of signal strength.
   @param muHat - Profiled signal strength.
   @param doTilde - True iff using qMuTilde instead of qMu.
   @param direction =  -1 if N < 0, else +1.
   @return - The value of sigma.
*/
double TestStat::getSigma(double qMu, double mu, double muHat, bool doTilde) {
  if (qMu == 0) printer("getSigma: ERROR! Divide by qMu=0!",true); 
  
  if (mu > muHat) return sqrt((mu-muHat)*(mu-muHat) / qMu);
  else if (muHat < 0 && doTilde) return sqrt((mu*mu - 2*mu*muHat) / qMu);
  else return sqrt((mu-muHat)*(mu-muHat) / qMu);
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of z0 based on the test statistic q0.
   @param q0 - The test statistic q0.
   @return - The value of z0 (significance).
*/
double TestStat::getZ0FromQ0(double q0) {
  return sqrt(q0);
}

/**
   -----------------------------------------------------------------------------
   Calculate the value of zMu based on the test statistic qMu.
   @param qMu - The test statistic qMu.
   @return - The value of zMu (significance).
*/
double TestStat::getZMuFromQMu(double qMu) {
  return sqrt(qMu);
}

/**
   -----------------------------------------------------------------------------
   Calculate Z (significance) based on the p-value.
   @param p - The p-value.
   @return - The significance (# standard deviations).
*/
double TestStat::getZFromP(double p) {
  return TMath::NormQuantile(1.0 - p);
}

/**
   -----------------------------------------------------------------------------
   Get the intersection point for the graph.
   @param graph - The graph for which intercept points will be found.
   @param valueToIntercept - The y-value to intercept.
   @return - The x-value of the intercept.
*/
double TestStat::graphIntercept(TGraph *graph, double valueToIntercept) {
  
  // Loop over points in the graph to get search range:
  double rangeMin = 0.0; double rangeMax = 0.0;
  for (int i_p = 0; i_p < graph->GetN(); i_p++) {
    double xCurr; double yCurr;
    graph->GetPoint(i_p, xCurr, yCurr);
    if (i_p == 0) rangeMin = xCurr;
    if (i_p == (graph->GetN()-1)) rangeMax = xCurr;
    // Also exit prematurely if only one point:
    if (graph->GetN() < 2) return xCurr;
  }
  
  // Determine whether slope > 0 (increasing) or < 0 (decreasing):
  bool increasing = true;
  if (graph->Eval(rangeMin) > graph->Eval(rangeMax)) increasing = false;
  
  // Bisection method to search for intercept:
  double precision = 0.0001;
  int nIterations = 0;
  int maxIterations = 30;
  double stepSize = (rangeMax - rangeMin) / 2.0;
  double currXValue = (rangeMax + rangeMin) / 2.0;
  double currYValue = graph->Eval(currXValue);
  while ((fabs(currYValue - valueToIntercept)/valueToIntercept) > precision && 
	 nIterations <= maxIterations) {
    
    currYValue = graph->Eval(currXValue);
    
    nIterations++;
    stepSize = 0.5 * stepSize;

    if (increasing) {
      if (currYValue > valueToIntercept) currXValue -= stepSize;
      else currXValue += stepSize;
    }
    else {
      if (currYValue > valueToIntercept) currXValue += stepSize;
      else currXValue -= stepSize;
    }
  }
  
  // Print error message and return bad value if convergence not achieved:
  if (nIterations == maxIterations) {
    std::cout << "TestStat: ERROR! Intercept not found." << std::cout;
    return -999;
  }
  
  return currXValue;
}

/**
   -----------------------------------------------------------------------------
   Takes in a variable declaration such as "name[1,0,5]" and returns "name".
   @param varForm - The form of the variable declared.
   @return - The name of the variable without the rest of the expression.
*/
TString TestStat::nameOfVar(TString varForm) {
  TString name = varForm;
  name.Remove(name.First("["));
  return name;
}

/**
   -----------------------------------------------------------------------------
   Get the most recent NLL scan graph.
   @return - The most recent NLL scan TGraph.
*/
TGraph *TestStat::nllScanGraph() {
  return m_graphNLL;
}

/**
   -----------------------------------------------------------------------------
   Create a ratio plot:
   @param dataName - The name of the RooAbsData set in the workspace.
   @param pdfName - The name of the RooAbsPdf in the workspace.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* TestStat::plotDivision(TString dataName, TString pdfName, 
				       TString obsName, double xMin, 
				       double xMax, double xBins){
  printer(Form("TestStat::plotDivision(%s, %s, %s, %f, %f, %f)",
	       dataName.Data(),pdfName.Data(),obsName.Data(),xMin,xMax,xBins),
	  false);
  
  RooRealVar *observable = m_workspace->var(obsName);
  RooAbsData *data = m_workspace->data(dataName);
  RooAbsPdf *pdf = m_workspace->pdf(pdfName); 
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  double nEvents = data->sumEntries();
    
  observable->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    observable->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable), 
					 RooFit::NormSet(*observable), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = observable->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight / currPdfWeight;
    if (currDataWeight == 0) currWeight = 1.0;
    result->SetPoint(pointIndex, currMass, currWeight);
    
    double currError = originHist->GetBinError(pointIndex+1) / currPdfWeight;
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  return result;
}

/**
   -----------------------------------------------------------------------------
   Plot the fits produced by the specified model.
   @param fitType - The type of fit.
   @param datasetName - The name of the dataset to be plotted. 
*/
void TestStat::plotFits(TString fitType, TString datasetName) {
  printer(Form("TestStat:plotFits(%s, %s)",fitType.Data(),datasetName.Data()),
	  false);
  
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  TPad *pad1 = new TPad( "pad1", "pad1", 0.00, 0.33, 1.00, 1.00 );
  TPad *pad2 = new TPad( "pad2", "pad2", 0.00, 0.00, 1.00, 0.33 );
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->SetBorderMode(0);
  
  can->cd();
  pad1->Draw();
  pad2->Draw();

  RooSimultaneous* combPdf = (RooSimultaneous*)m_mc->GetPdf();
  RooArgSet* observables = (RooArgSet*)m_mc->GetObservables();
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  while ((cateType = (RooCatType*)cateIter->Next())) {
    pad1->cd();
    pad1->Clear();
    RooAbsData *currData 
      = m_workspace->data(datasetName)->getSimData(cateType->GetName());
    RooAbsPdf *currPdf = combPdf->getPdf(cateType->GetName());
    RooArgSet *currObsSet = currPdf->getObservables(observables);
    RooRealVar *currObs = (RooRealVar*)currObsSet->first();
    
    // Set the resonant analysis plot binning and axis scale to paper settings:
    int nBinsForPlot = currObs->getBins();
    //      = (int)((currObs->getMax() - currObs->getMin()) / m_geVPerBin);
    
    // Plot everything on RooPlot:
    RooPlot* frame = currObs->frame(nBinsForPlot);
    currData->plotOn(frame);
    currPdf->plotOn(frame, LineColor(2), LineStyle(1));
    frame->SetXTitle(currObs->GetName());
    frame->SetYTitle(Form("Events / %2.1f GeV", m_geVPerBin));
    frame->Draw();
    
    if (m_useLogScale) {
      gPad->SetLogy();
      frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    }
    else frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    
    // Print ATLAS text on the plot:    
    TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
    t.SetTextFont(72); t.SetTextSize(0.05);
    t.DrawLatex(0.20, 0.88, "ATLAS");
    t.SetTextFont(42); t.SetTextSize(0.05);
    t.DrawLatex(0.32, 0.88, m_config->getStr("ATLASLabel"));
    t.DrawLatex(0.20, 0.82, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
				 (m_config->getNum("AnalysisLuminosity")/
				  1000.0)));
    
    TString printCateName 
      = m_config->getStr(Form("PrintCateName_%s", cateType->GetName()));
    t.DrawLatex(0.20, 0.76, printCateName);
    
    // Second pad on the canvas (Ratio Plot):
    pad2->cd();
    pad2->Clear();
    
    double ratioMin = 0.0;//-0.2
    double ratioMax = 2.0;//2.2
    TH1F *medianHist = new TH1F("median", "median", nBinsForPlot, 
				currObs->getMin(), currObs->getMax());
    for (int i_b = 1; i_b <= nBinsForPlot; i_b++) {
      if (m_ratioOrSubtraction.EqualTo("Ratio")) {
	medianHist->SetBinContent(i_b, 1.0);
      }
      else {
	medianHist->SetBinContent(i_b, 0.0);
      }
    }
    medianHist->SetLineColor(kRed);
    medianHist->SetLineWidth(2);
    medianHist->GetXaxis()->SetTitle(currObs->GetName());
    if (m_ratioOrSubtraction.EqualTo("Ratio")) {
      medianHist->GetYaxis()->SetTitle("Data / Fit");
    }
    else {
      medianHist->GetYaxis()->SetTitle("Data - Fit");
    }
    medianHist->GetXaxis()->SetTitleOffset(0.95);
    medianHist->GetYaxis()->SetTitleOffset(0.7);
    medianHist->GetXaxis()->SetTitleSize(0.1);
    medianHist->GetYaxis()->SetTitleSize(0.1);
    medianHist->GetXaxis()->SetLabelSize(0.1);
    medianHist->GetYaxis()->SetLabelSize(0.1);
    if (m_ratioOrSubtraction.EqualTo("Ratio")) {
      medianHist->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    }
    medianHist->GetYaxis()->SetNdivisions(4);
    if (m_ratioOrSubtraction.EqualTo("Ratio")) medianHist->Draw();
    
    TGraphErrors* subData = NULL;
    if (m_ratioOrSubtraction.EqualTo("Ratio")) {
      subData = plotDivision(currData->GetName(), currPdf->GetName(),
			     currObs->GetName(), currObs->getMin(),
			     currObs->getMax(), nBinsForPlot);
      
      TLine *line = new TLine();
      line->SetLineStyle(1);
      line->SetLineWidth(2);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->SetLineStyle(2);
      line->DrawLine(currObs->getMin(), ((1.0+ratioMin)/2.0),
		     currObs->getMax(), ((1.0+ratioMin)/2.0));
      line->DrawLine(currObs->getMin(), ((1.0+ratioMax)/2.0),
		     currObs->getMax(), ((1.0+ratioMax)/2.0));
      
      subData->Draw("EPSAME");
    }
    else {
      subData = plotSubtract(currData->GetName(), currPdf->GetName(),
			     currObs->GetName(), currObs->getMin(),
			     currObs->getMax(), nBinsForPlot);
      medianHist->GetYaxis()->SetRangeUser(subData->GetYaxis()->GetXmin(),
					   subData->GetYaxis()->GetXmax());
      medianHist->Draw();
      subData->Draw("EPSAME");
    }
    
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_anaType.Data(), fitType.Data(), cateType->GetName()));
    can->Print(Form("%s/fitPlot_%s_%s_%s.C", m_plotDir.Data(),
		    m_anaType.Data(), fitType.Data(), cateType->GetName()));
    
    delete frame;
  }
  delete can;
}

/**
   -----------------------------------------------------------------------------
   Create a subtraction plot:
   @param dataName - The name of the RooAbsData set in the workspace.
   @param pdfName - The name of the RooAbsPdf in the workspace.
   @param xMin - The minimum value of the observable range.
   @param xMax - The maximum value of the observable range.
   @param xBins - The number of bins for the observable.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* TestStat::plotSubtract(TString dataName, TString pdfName, 
				     TString obsName, double xMin, double xMax, 
				     double xBins){
  printer(Form("TestStat::plotSubtract(%s, %s, %s, %f, %f, %f)",
	       dataName.Data(),pdfName.Data(),obsName.Data(),xMin,xMax,xBins),
	  false);
  
  RooRealVar *observable = m_workspace->var(obsName);
  RooAbsData *data = m_workspace->data(dataName);
  RooAbsPdf *pdf = m_workspace->pdf(pdfName); 
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  double nEvents = data->sumEntries();
    
  observable->setRange("fullRange", xMin, xMax);
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
  				   RooFit::Binning(xBins, xMin, xMax));
  TGraphErrors *result = new TGraphErrors();
  double increment = (xMax - xMin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = xMin; i_m < xMax; i_m += increment) {
    observable->setRange(Form("range%2.2f",i_m), i_m, (i_m+increment));
    RooAbsReal* intCurr
      = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable), 
					 RooFit::NormSet(*observable), 
					 RooFit::Range(Form("range%2.2f",i_m)));
    double valCurr = intCurr->getVal();
    
    double currMass = i_m + (0.5*increment);
    double currPdfWeight = nEvents * (valCurr / valTot);
    TString varName = observable->GetName();
    double currDataWeight = data->sumEntries(Form("%s>%f&&%s<%f",varName.Data(),
						  i_m,varName.Data(),
						  (i_m+increment)));
    double currWeight = currDataWeight - currPdfWeight;
    if (currDataWeight > 0.000001) {
      result->SetPoint(pointIndex, currMass, currWeight);
    }
    
    double currError = originHist->GetBinError(pointIndex+1);
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  return result;
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void TestStat::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Print the names and values of parameters in a RooArgSet:
   @param setName - The name of the RooArgSet.
   @param set - The RooArgSet.
*/
void TestStat::printSet(TString setName, RooArgSet* set) {
  printer(Form("TestStat::printSet(%s)",setName.Data()), false);
  TIterator *iterSet = set->createIterator();
  RooRealVar *curr = NULL;
  // Bin each of the observables:
  while ((curr = (RooRealVar*)iterSet->Next())) {
    TString currLine = Form("\t %s = %f", ((TString)curr->GetName()).Data(),
			    curr->getVal());
    printer(currLine, false);
  }
}

/**
   -----------------------------------------------------------------------------
   Choose whether or not to reset nuisance parameters after calling getFitNLL().
   @param doResetParamsAfterFit - True iff. parameters should be set to original
   values after fitting.
*/
void TestStat::resetParamsAfterFit(bool doResetParamsAfterFit) {
  m_doResetParamsAfterFit = doResetParamsAfterFit;
}

/**
   -----------------------------------------------------------------------------
   Choose whether or not to save snapshots from profiling data.
   @param doSaveSnapshot - True iff you want to save snapshots in future fits.
*/
void TestStat::saveSnapshots(bool doSaveSnapshot) {
  m_doSaveSnapshot = doSaveSnapshot;
}

/**
   -----------------------------------------------------------------------------
   Use if Asimov data should be scaled. Set the scale factor and the variable
   names for which this scale factor should be applied.
*/
void TestStat::scaleAsimovData(double scaleFactor,
			       std::vector<TString> varsToScale) {
  m_AsimovScaleFactor = scaleFactor;
  m_AsimovVarsToScale = varsToScale;
}

/**
   -----------------------------------------------------------------------------
   Scan the NLL of a particular variable.
   @param scanName - The name of the scan (for saving plots...).
   @param datasetName - The name of the dataset in the workspace.
   @param varToScan - The variable that will be scanned.
   @param fixPoI - True if PoI should be fixed to the specified value.
   @param namesAndValsPoI - Map of names and values of PoIs to set for fit.
   @param nScanPoints - The number of points to use in the scan (more=better).
   @return - The -1 sigma, median, and +1 sigma variable values.
*/
std::map<int,double> TestStat::scanNLL(TString scanName, TString datasetName,
				       TString varToScan, bool fixPoI,
				       std::map<TString,double> namesAndValsPoI,
				       int nScanPoints) {
  printer(Form("TestStat::scanNLL(datasetName = %s, varToScan = %s)",
	       datasetName.Data(), varToScan.Data()), false);
  
  // A map to store the scan results (median = 0, sigmas are +/-1, +/-2, etc...)
  std::map<int,double> result; result.clear();
    // Load objects for fitting:
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* origValNP
    = (RooArgSet*)m_workspace->getSnapshot(m_nominalSnapshot);
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstPoI = (RooRealVar*)poi->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  m_workspace->loadSnapshot(m_nominalSnapshot);
  
  // Look for dataset:
  if (!m_workspace->data(datasetName)) {
    printer(Form("TestStat: Error! Requested data not available: %s",
		 datasetName.Data()), true);
  }

  // Release nuisance parameters before fit and set to the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // The global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
  // Set the values of parameters of interest as specified in the input.
  for (std::map<TString,double>::iterator iterPoI = namesAndValsPoI.begin();
       iterPoI != namesAndValsPoI.end(); iterPoI++) {
    if (m_workspace->var(iterPoI->first)) {
      m_workspace->var(iterPoI->first)->setVal(iterPoI->second);
      m_workspace->var(iterPoI->first)->setConstant(fixPoI);
    }
    else {
      printer(Form("TestStat: ERROR! %s parameter missing",
		   (iterPoI->first).Data()), true);
    }
  }
  
  // Check if other parameter settings have been specified for fit:
  for (std::map<TString,double>::iterator iterParam = m_paramValToSet.begin();
       iterParam != m_paramValToSet.end(); iterParam++) {
    // Check that workspace contains parameter:
    if (m_workspace->var(iterParam->first)) {
      m_workspace->var(iterParam->first)->setVal(iterParam->second);
      m_workspace->var(iterParam->first)
	->setConstant(m_paramConstToSet[iterParam->first]);
    }
    else {
      m_workspace->Print("v");
      printer(Form("TestStat: Error! Parameter %s not found in workspace!",
		   ((TString)iterParam->first).Data()), true);
    }
  }
  
  // Make sure variable to scan is free in the fit:
  m_workspace->var(varToScan)->setConstant(false);

  // First do a fit to get the unconditional maximum likelihood fit result:
  RooNLLVar* varMLNLL
    = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				     Extended(combPdf->canBeExtended()));
  RooFitResult *fitResultML
    = statistics::minimize(varMLNLL, m_fitOptions, NULL, true);
  if (!fitResultML || fitResultML->status() != 0) m_allGoodFits = false;
  m_workspace->saveSnapshot(Form("paramsProfileML_%s",scanName.Data()),
			    *poiAndNuis);
  double minNLLVal = varMLNLL->getVal();
  delete fitResultML;
  delete varMLNLL;
  
  // Then retrieve the central value and errors as preliminary result:
  double errorLo = m_workspace->var(varToScan)->getError();
  double errorHi = m_workspace->var(varToScan)->getError();
  result[0] = m_workspace->var(varToScan)->getVal();
  result[-2] = result[0] - (2.0 * errorLo);
  result[-1] = result[0] - errorLo;
  result[1] = result[0] + errorHi;
  result[2] = result[0] + (2.0 * errorHi);
  
  double scanMin = result[0] - (3.0 * errorLo);
  double scanMax = result[0] + (3.0 * errorHi);
  /*
  if (scanMin < m_workspace->var(varToScan)->getMin()) {
    scanMin = m_workspace->var(varToScan)->getMin();
  }
  if (scanMax > m_workspace->var(varToScan)->getMax()) {
    scanMax = m_workspace->var(varToScan)->getMax();
  }
  */
  
  // Compile a list of points to scan:
  std::vector<double> scanValues; scanValues.clear();
  for (double currValue = scanMin; currValue <= scanMax; 
       currValue += ((scanMax-scanMin)/nScanPoints)){
    scanValues.push_back(currValue);
  }
  scanValues.push_back(result[0]);
  std::sort(scanValues.begin(), scanValues.end());
  
  //----------------------------------------//
  // Scan the variable of interest, and perform a fit at each point to get NLL:
  
  TGraph *gNLL = new TGraph(); 
  gNLL->SetNameTitle(Form("tmp%s",scanName.Data()),
		     Form("tmp%s",scanName.Data()));
  TGraph *gNLLLo = new TGraph();
  TGraph *gNLLHi = new TGraph();
  int indexLo = 0;
  int indexHi = 0;
  
  for (int scanIndex = 0; scanIndex < (int)scanValues.size(); scanIndex++) {
    
    // Set the value of the variable to scan:
    m_workspace->var(varToScan)->setVal(scanValues[scanIndex]);
    m_workspace->var(varToScan)->setConstant(true);
    
    // The actual fit command:
    RooNLLVar* currNLL
      = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				       Extended(combPdf->canBeExtended()));
    
    RooFitResult *currFitResult
      = statistics::minimize(currNLL, m_fitOptions, NULL, true);
    if (!currFitResult || currFitResult->status() != 0) m_allGoodFits = false;
    double currNLLValue = currNLL->getVal();
    delete currFitResult;
    delete currNLL;
    
    double deltaNLL = 2 * (currNLLValue - minNLLVal);
    // Add to the scan graph:
    gNLL->SetPoint(scanIndex, scanValues[scanIndex], deltaNLL);
    
    if (scanValues[scanIndex] <= result[0]) {
      gNLLLo->SetPoint(indexLo, scanValues[scanIndex], deltaNLL);
      indexLo++;
    }
    if (scanValues[scanIndex] >= result[0]) {
      gNLLHi->SetPoint(indexHi, scanValues[scanIndex], deltaNLL);
      indexHi++;
    }
    printer(Form("TestStat: NLL at %f is %f (%f - %f).",
		 scanValues[scanIndex], currNLLValue-minNLLVal,
		 currNLLValue, minNLLVal), false);
  }
  
  // Release nuisance parameters after fit and recover the default values:
  statistics::constSet(nuisanceParameters, false, origValNP);
  
  // Calculate the actual intercept values:
  result[-2] = graphIntercept(gNLLLo, 4.0);
  result[-1] = graphIntercept(gNLLLo, 1.0);
  result[1] = graphIntercept(gNLLHi, 1.0);
  result[2] = graphIntercept(gNLLHi, 4.0);

  if (m_doPlot) {
    // Draw the graph:
    TCanvas *can = new TCanvas("can", "can", 800, 600);
    can->cd();
    gNLL->SetLineColor(kRed+1);
    gNLL->SetLineWidth(2);
    gNLL->GetXaxis()->SetTitle(varToScan);
    gNLL->GetYaxis()->SetTitle("-2#DeltaNLL");
    gNLL->GetYaxis()->SetRangeUser(0, gNLL->GetMaximum());
    gNLL->Draw("AL");
    
    // 1 sigma and 2 sigma lines:
    TLine *line = new TLine();
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->DrawLine(gNLL->GetXaxis()->GetXmin(), 1.0,
		   gNLL->GetXaxis()->GetXmax(), 1.0);
    line->DrawLine(gNLL->GetXaxis()->GetXmin(), 4.0,
		   gNLL->GetXaxis()->GetXmax(), 4.0);
    
    // Print ATLAS text on the plot:    
    TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
    t.SetTextFont(72); t.SetTextSize(0.05);
    t.DrawLatex(0.30, 0.86, "ATLAS");
    t.SetTextFont(42); t.SetTextSize(0.05);
    t.DrawLatex(0.42, 0.86, m_config->getStr("ATLASLabel"));
    t.DrawLatex(0.30, 0.80, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
				 (m_config->getNum("AnalysisLuminosity")/
				  1000.0)));
    t.DrawLatex(0.30, 0.74, Form("%s = %2.2f_{-%2.2f}^{+%2.2f}", 
				 varToScan.Data(), result[0], result[-1],
				 result[1]));
    
    // Print the canvas:
    can->Print(Form("%s/scanNLL_%s_%s.eps", m_plotDir.Data(), scanName.Data(),
		    m_anaType.Data()));
    can->Print(Form("%s/scanNLL_%s_%s.C", m_plotDir.Data(), scanName.Data(),
		    m_anaType.Data()));
    delete line;
    delete can;
  }
  
  // Also return the NLL scan graph:
  m_graphNLL = gNLL;

  // Finish up, return NLL value.
  printer("TestStat: NLL Scan has completed.", false);
  delete gNLLLo;
  delete gNLLHi;
  
  // Return the median, +/-1, and +/-2 sigma values:
  return result;
}

/**
   -----------------------------------------------------------------------------
   Choose the options for the fit, to be passed to "statistics" class minimizer.
   @param fitOptions - The options for the fits done by this class.
*/
void TestStat::setFitOptions(TString fitOptions) {
  m_fitOptions = fitOptions;
}

/**
   -----------------------------------------------------------------------------
   Set the name of the snapshot to load by default in calls to getFitNLL()
   @param nominalSnapshot - The name of the nominal snapshot to load before
   running getFitNLL();
*/
void TestStat::setNominalSnapshot(TString nominalSnapshot) {
  m_nominalSnapshot = nominalSnapshot;
}

/**
   -----------------------------------------------------------------------------
   Choose to plot ratio or subtraction plot below plot of mass points.
   @param ratioOrSubtraction - "Ratio" or "Subtraction" plot below main plot.
*/
void TestStat::setSubPlot(TString ratioOrSubtraction) {
  m_ratioOrSubtraction = ratioOrSubtraction;
}

/**
   -----------------------------------------------------------------------------
   Set axis options for plots:
*/
void TestStat::setPlotAxis(bool useLogScale, double yMin, double yMax,
			     double GeVPerBin) {
  m_useLogScale = useLogScale;
  m_yMin = yMin;
  m_yMax = yMax;
  m_geVPerBin = GeVPerBin;
}

/**
   -----------------------------------------------------------------------------
   Set an output directory and enable plotting.
   @param directory - The output directory path.
*/
void TestStat::setPlotDirectory(TString directory) {
  system(Form("mkdir -vp %s", directory.Data()));
  m_plotDir = directory;
  m_doPlot = true;
}

/**
   -----------------------------------------------------------------------------
   Set the named parameter to a certain value and either fix or free it.
   @param paramName - The name of the fit parameter.
   @param paramVal - The new value of the fit parameter.
   @param doSetConstant - True iff the parameter should be set constant. 
*/
void TestStat::setParam(TString paramName, double paramVal,
			  bool doSetConstant) {
  m_paramValToSet[paramName] = paramVal;
  m_paramConstToSet[paramName] = doSetConstant;
  
  // Do this immediately?
  m_workspace->var(paramName)->setVal(paramVal);
  m_workspace->var(paramName)->setConstant(doSetConstant);
}

/**
   -----------------------------------------------------------------------------
   Store the parameter names and values.
   @param set - The RooArgSet to store.
   @param map - The map to hold the values.
*/
void TestStat::storeParams(RooArgSet *set, std::map<std::string,double>& map){
  map.clear();
  TIterator *iterSet = set->createIterator();
  RooRealVar *curr = NULL;
  while ((curr = (RooRealVar*)iterSet->Next())) {
    map[(std::string)curr->GetName()] = curr->getVal();
  }
}

/**
   -----------------------------------------------------------------------------
   Get the workspace.
   @return - A pointer to the class RooWorkspace object.
*/
RooWorkspace *TestStat::theWorkspace() {
  return m_workspace;
}

/**
   -----------------------------------------------------------------------------
   Get the ModelConfig.
   @return - A pointer to the class ModelConfig object.
*/

ModelConfig *TestStat::theModelConfig() {
  return m_mc;
}

/**
   -----------------------------------------------------------------------------
   Bisection method to search for 95% CL limit intercept. Basically adjust mu 
   for the numerator in qMu (lMu/lfree) until CLs=0.05. 
   @param mapPoI - Map of names and values of mu=1 PoIs to set for CL.
   @param datasetName - The name of the dataset to fit.
   @param poiForNorm - The name of the normalization poi.
   @param doTilde - True iff using qMuTilde instead of qMu.
   @return - The upper-limit on the signal strength. 
*/
double TestStat::upperLimitFinder(std::map<TString,double> mapPoI, 
				  TString datasetName, TString poiForNorm,
				  bool doTilde) {
  printer(Form("TestStat::upperLimitFinder(%s, %s)",
	       datasetName.Data(), poiForNorm.Data()), false);
  //
  // Also add tilde option to everything (qMu -> qMuTilde, etc.)
  //
  
  // Calculate nllMuHat only once:
  double nllMuHat = getFitNLL(datasetName, 0.0, false, mapPoI, false);
  double profiledNorm = (getPoIs())[(std::string)(poiForNorm)];
  
  double precision = 0.0001;
  int nIterations = 0;
  int maxIterations = 30;
  double stepSize = (m_workspace->var(poiForNorm)->getMax() - 
		     m_workspace->var(poiForNorm)->getMin()) / 2.0;
  double muLimit = (m_workspace->var(poiForNorm)->getMax() + 
		    m_workspace->var(poiForNorm)->getMin()) / 2.0;
  double CLs = 0.0;
  double intercept = 0.05;// For 95% exclusion
  while ((((CLs - intercept) / intercept) > precision) && 
	 (nIterations <= maxIterations)) {
    std::cout << "TestStat:asymptoticLimitBisector: Iteration " << nIterations 
	      << " starting, muLimit=" << muLimit << ", CLs=" << CLs
	      << std::endl;
    
    // Set the signal strength:
    mapPoI[poiForNorm] = muLimit;
    
    // For  now, use mu = 1 (qmu for mu=1), since we are adjusting the cross-
    // section parameter and not the mu parameter.
    double muForQMu = 1.0;
    
    // Calculate qMu:
    double nllMu1 = getFitNLL(datasetName, 1.0, true, mapPoI, false);
    double muHat = (profiledNorm / mapPoI[poiForNorm]);
    double qMu = getQMuFromNLL(nllMu1, nllMuHat, muHat, muForQMu);
    
    // Replace qMu with qMuTilde if requested (requires extra fit for nllMu0):
    if (doTilde) {
      mapPoI[poiForNorm] = 0.0;
      double nllMu0 = getFitNLL(datasetName, 0.0, true, mapPoI, false);
      mapPoI[poiForNorm] = muLimit;
      qMu = getQMuTildeFromNLL(nllMu1, nllMu0, nllMuHat, muHat, muForQMu);
    }
    
    // Get Sigma:
    //int direction = (N < 0) ? -1 : 1;
    double sigma = getSigma(qMu, muForQMu, muHat, doTilde);
    
    // Which Mu to use? Probably muHat, right?
    CLs = getCLsFromQMu(qMu, sigma, muForQMu);

    // Increase the iteration count and reduce the step size:
    nIterations++;
    stepSize = 0.5 * stepSize;
    
    // Update the value of the mu limit:
    if (CLs > intercept) muLimit -= stepSize;
    else muLimit += stepSize;
  } // At this point, muLimit should be the 95% CL limit on signal strength.
  
  if (nIterations == maxIterations) {
    printer("TestStat::asymptoticLimitBisector: Maximum Iterations!", true);
  }
  
  return muLimit;
}

/**
   -----------------------------------------------------------------------------
   Option to use single or two-sided test statistics for q0 and qMu.
   @param useTwoSided - True if test statistics should be two-sided (e.g. r0 
   instead of q0, rMu instead of qMu...).
*/
void TestStat::useTwoSidedTestStat(bool useTwoSided) {
  m_useTwoSided = useTwoSided;
}
