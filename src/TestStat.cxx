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
  
  // Map storing all calculations:
  m_calculatedValues.clear();
  
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
  
  // Make new or load old values:
  if (m_options.Contains("FromFile")) loadStatsFromFile();
  
  // Finished instantiating the test statistic class. 
  std::cout << "TestStat: Initialized Successfully!" << std::endl;
  return;
}

/**
   -----------------------------------------------------------------------------
   Get the value of one of the test statistics.
   @param testStat - The test stat. name (p0, CL, CLs).
   @param observed - True iff observed, false if expected. 
   @param N - The standard deviation (-2, -1, 0, +1, +2). 
   @return - The test statistic value.
*/
double TestStat::accessValue(TString testStat, bool observed, int N) {
  TString currMapKey = getKey(testStat, observed, N);
  // Check that corresponding entry exists:
  if (mapValueExists(currMapKey)) return m_calculatedValues[currMapKey];
  else return 0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the CL and CLs values using model fits.
   @param type - "Exp" "Obs" or "Both".

void TestStat::calculateNewCL() {
  printer("TestStat::calculateNewCL()", false);
  
  // Calculate observed qmu: 
  double muHatObs = 0.0;
  double nllMu1Obs = getFitNLL(m_dataForObsQMu, 1.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQMu, 1.0, false, muHatObs);
  double obsQMu = getQMuFromNLL(nllMu1Obs, nllMuHatObs, muHatObs, 1);
  
  // Calculate expected qmu:
  double muHatExp = 0.0;
  double nllMu1Exp = getFitNLL(m_dataForExpQMu, 1.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQMu, 0.0, false, muHatExp);
  double expQMu = getQMuFromNLL(nllMu1Exp, nllMuHatExp, muHatExp, 1);
  
  // Calculate CL:
  double expCLn2 = getCLFromQMu(expQMu, -2);
  double expCLn1 = getCLFromQMu(expQMu, -1);
  double expCLp1 = getCLFromQMu(expQMu, 1);
  double expCLp2 = getCLFromQMu(expQMu, 2);
  double expCL = getCLFromQMu(expQMu, 0);
  double obsCL = getCLFromQMu(obsQMu, 0);
  
  // Write CL values to file:
  ofstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt",
		   m_outputDir.Data(), m_anaType.Data()));
  textCL << m_anaType << " " << obsCL << " " << expCLn2 << " " << expCLn1
	 << " " << expCL << " " << expCLp1 << " " << expCLp2 << std::endl;
  textCL.close();
  
  // Print summary:
  std::cout << "  expected CL +2s = " << expCLp2 << std::endl;
  std::cout << "  expected CL +1s = " << expCLp1 << std::endl;
  std::cout << "  expected CL nom = " << expCL    << std::endl;
  std::cout << "  expected CL -1s = " << expCLn1 << std::endl;
  std::cout << "  expected CL -2s = " << expCLn2 << std::endl;
  std::cout << "  observed CL = " << obsCL << std::endl;
  std::cout << " " << std::endl;
  if (m_allGoodFits) std::cout << "All good fits? True" << std::endl;
  else std::cout << "All good fits? False" << std::endl;
  cout << " " << endl;
  if (obsQMu < 0) std::cout << "WARNING! obsQMu < 0 : " << obsQMu << std::endl;
  if (expQMu < 0) std::cout << "WARNING! expQMu < 0 : " << expQMu << std::endl;
  
  // Save CL and CLs for later access:
  m_calculatedValues[getKey("CL",0,-2)] = expCLn2;
  m_calculatedValues[getKey("CL",0,-1)] = expCLn1;
  m_calculatedValues[getKey("CL",0,0)] = expCL;
  m_calculatedValues[getKey("CL",0,1)] = expCLp1;
  m_calculatedValues[getKey("CL",0,2)] = expCLp2;
  m_calculatedValues[getKey("CL",1,0)] = obsCL;
  
  m_calculatedValues[getKey("CLs",0,-2)] = getCLsFromCL(expCLn2);
  m_calculatedValues[getKey("CLs",0,-1)] = getCLsFromCL(expCLn1);
  m_calculatedValues[getKey("CLs",0,0)] = getCLsFromCL(expCL);
  m_calculatedValues[getKey("CLs",0,1)] = getCLsFromCL(expCLp1);
  m_calculatedValues[getKey("CLs",0,2)] = getCLsFromCL(expCLp2);
  m_calculatedValues[getKey("CLs",1,0)] = getCLsFromCL(obsCL);
}
*/

/**
   -----------------------------------------------------------------------------
   Calculate the p0 value using model fits.

void TestStat::calculateNewP0() {
  printer("TestStat::calculateNewP0()", false);
  
  // Calculate observed q0: 
  double muHatObs = 0.0;
  double nllMu0Obs = getFitNLL(m_dataForObsQ0, 0.0, true, muHatObs);
  double nllMuHatObs = getFitNLL(m_dataForObsQ0, 0.0, false, muHatObs);
  double obsQ0 = getQ0FromNLL(nllMu0Obs, nllMuHatObs, muHatObs);
  
  // Calculate expected q0:
  double muHatExp = 0.0;
  double nllMu0Exp = getFitNLL(m_dataForExpQ0, 0.0, true, muHatExp);
  double nllMuHatExp = getFitNLL(m_dataForExpQ0, 0.0, false, muHatExp);
  double expQ0 = getQ0FromNLL(nllMu0Exp, nllMuHatExp, muHatExp);
  
  // Calculate p0 from q0:
  double expP0 = getP0FromQ0(expQ0);
  double obsP0 = getP0FromQ0(obsQ0);
  
  // Write p0 values to file:
  ofstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  textP0 << m_anaType << " " << expP0 << " " << obsP0 << std::endl;
  textP0.close();
  
  // Print summary:
  std::cout << "\n  Expected p0 = " << expP0 << std::endl;
  std::cout << "\n  Observed p0 = " << obsP0 << std::endl;
  if (fitsAllConverged()) std::cout << "All good fits? True\n" << std::endl;
  else std::cout << "All good fits? False\n" << std::endl;
  
  // Save p0 for later access:
  m_calculatedValues[getKey("p0", 1, 0)] = obsP0;
  m_calculatedValues[getKey("p0", 0, 0)] = expP0;
}
*/

/**
   -----------------------------------------------------------------------------
   Clears all data stored by the class, but does not modify the workspace.
*/
void TestStat::clearData() {
  m_allGoodFits = true;
  m_calculatedValues.clear();
  m_mapGlobs.clear();
  m_mapNP.clear();
  m_mapPoI.clear();
  m_doSaveSnapshot = false;
  m_doPlot = false;
  m_plotDir = "";
  clearFitParamSettings();
  m_graphNLL = NULL;
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
   Create a pseudo-dataset with a given value of DM and SM signal strength.
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
    //m_workspace->loadSnapshot("paramsOrigin");
    printer(Form("TestStat: ERROR! No snapshot %s", snapshotName.Data()),
	    true);
  }
  
  //RooAbsPdf* combPdf = (RooAbsPdf*)m_mc->GetPdf();
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
      TIterator *iterObs = currObs->createIterator();
      RooRealVar *currObs = NULL;
      
      // Bin each of the observables:
      while ((currObs = (RooRealVar*)iterObs->Next())) {
	currObs->setBins(m_config->getInt("WorkspaceObsBins"));
      }
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs, AutoBinned(true),
					 Extended(currPdf->canBeExtended()),
					 GenBinned("PleaseGenerateBinned"));
    }
    // Construct unbinned pseudo-data by default:
    else {
      toyDataMap[(std::string)cateType->GetName()]
	= (RooDataSet*)currPdf->generate(*currObs,Extended(true));
    }
    
    double currEvt = toyDataMap[(std::string)cateType->GetName()]->sumEntries();
    m_numEventsPerCate.push_back(currEvt);
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
  if (toyIndex < 0) {
    pseudoData = new RooDataSet("toyData", "toyData", *observables, 
				RooFit::Index(*categories),
				RooFit::Import(toyDataMap));
  }
  else {
    pseudoData = new RooDataSet(Form("toyData%d",toyIndex),
				Form("toyData%d",toyIndex), *observables, 
				RooFit::Index(*categories),
				RooFit::Import(toyDataMap));
  }
  
  // Save the parameters used to generate toys:
  storeParams(nuisanceParameters, m_mapNP);
  storeParams(globalObservables, m_mapGlobs);
  storeParams(poi, m_mapPoI);

  // release nuisance parameters (but don't change the values!):
  //m_workspace->loadSnapshot("paramsOrigin");
  statistics::constSet(nuisanceParameters, false);
  statistics::constSet(globalObservables, true);
  
  // Import into the workspace then return:
  m_workspace->import(*pseudoData);
  return pseudoData;
}

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
   Implements the functional form of qMu.
   @param x - The value of the test statistic.
   @return - The value of the asymptotic test statistic distribution.
*/
double TestStat::functionQMu(double x) {
  // This corresponds to the "special case" of mu=mu'
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
*/
double TestStat::getCLFromQMu(double qMu, double N) {
  double CL = getCLFromCLs(getCLsFromQMu(qMu, N));
  return CL;
}

/**
   -----------------------------------------------------------------------------
   Get the CLs value using qMu and the type.
   @param qMu - The value of the test statistic.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The CLs value.
*/
double TestStat::getCLsFromQMu(double qMu, double N) {
  // N = 0 for exp and obs
  double pMu = getPMuFromQMu(qMu);
  double pB = getPFromN(N);
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
  if (resetParams) m_workspace->loadSnapshot("paramsOrigin");
  RooArgSet* origValNP = (RooArgSet*)m_workspace->getSnapshot("paramsOrigin");
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
  
  RooFitResult *fitResult = statistics::minimize(varNLL, "", NULL, true);
  if (!fitResult || fitResult->status() != 0) m_allGoodFits = false;
    
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
   Get the key for the value map.
   @param testStat - The test statistic.
   @param observed - True iff. observed.
   @param N - The sigma value (-2,-1,0,1,2).
*/
TString TestStat::getKey(TString testStat, bool observed, int N) {
  TString currKey = testStat;
  
  if (observed) currKey += "_obs";
  else currKey += "_exp";
  
  if (N > 0) currKey += Form("_p%d",N);
  else if (N < 0) currKey += Form("_n%d",N);
  
  return currKey;
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
   Calculate the value of pMu.
   @param qMu - The test statistic qMu.
   @return - The value of pMu.
*/
double TestStat::getPMuFromQMu(double qMu) {
  return (1 - ROOT::Math::gaussian_cdf(sqrt(fabs(qMu))));
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
  double qMu = 0.0;
  if (muHat < muTest) qMu = 2 * (nllMu - nllMuHat);
  return qMu;
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
  double qMuTilde = 0;
  if (muHat <= 0) qMuTilde = 2 * (nllMu - nllMu0);
  else if (muHat > 0 && muHat <= muTest) qMuTilde = 2 * (nllMu - nllMuHat);
  else if (muHat > muTest) qMuTilde = 0;
  return qMuTilde;
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
  
  // 
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
    std::cout << "DHCLScan: ERROR! Intercept not found." << std::cout;
    return -999;
  }
  
  return currXValue;
}

/**
   -----------------------------------------------------------------------------
   Load the statistics files (p0 and CL) that were previously generated. If none
   are found, then create from scratch automatically.
*/
void TestStat::loadStatsFromFile() {
  
  // Load input p0 file:
  ifstream textP0;
  textP0.open(Form("%s/p0/p0_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  
  // Load input CL file:
  ifstream textCL;
  textCL.open(Form("%s/CL/CL_values_%s.txt", 
		   m_outputDir.Data(), m_anaType.Data()));
  
  // If the input files don't exist, create from scratch:
  if (!textCL || !textP0) {
    printer("TestStat: ERROR! Stats don't exist. Exiting.", true);
    //calculateNewCL();
    //calculateNewP0();
    return;
  }
  
  // Read p0 values:
  TString inName; double inExpP0; double inObsP0; 
  while (!textP0.eof()) {
    textP0 >> inName >> inExpP0 >> inObsP0;
    textP0.close();
  }
  m_calculatedValues[getKey("p0",1,0)] = inObsP0;
  m_calculatedValues[getKey("p0",0,0)] = inExpP0;
  
  // Read CL values:
  double inObsCL, inExpCLn2, inExpCLn1, inExpCL, inExpCLp1, inExpCLp2;
  while (!textCL.eof()) {
    textCL >> inName >> inObsCL >> inExpCLn2 >> inExpCLn1 >> inExpCL
	   >> inExpCLp1 >> inExpCLp2;
  }
  textCL.close();
  
  // save CL and CLs for later access:
  m_calculatedValues[getKey("CL", 0, -2)] = inExpCLn2;
  m_calculatedValues[getKey("CL", 0, -1)] = inExpCLn1;
  m_calculatedValues[getKey("CL", 0, 0)] = inExpCL;
  m_calculatedValues[getKey("CL", 0, 1)] = inExpCLp1;
  m_calculatedValues[getKey("CL", 0, 2)] = inExpCLp2;
  m_calculatedValues[getKey("CL", 1, 0)] = inObsCL;
  
  m_calculatedValues[getKey("CLs", 0, -2)] = getCLsFromCL(inExpCLn2);
  m_calculatedValues[getKey("CLs", 0, -1)] = getCLsFromCL(inExpCLn1);
  m_calculatedValues[getKey("CLs", 0, 0)] = getCLsFromCL(inExpCL);
  m_calculatedValues[getKey("CLs", 0, 1)] = getCLsFromCL(inExpCLp1);
  m_calculatedValues[getKey("CLs", 0, 2)] = getCLsFromCL(inExpCLp2);
  m_calculatedValues[getKey("CLs", 1, 0)] = getCLsFromCL(inObsCL);
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
    
  // Loop over categories:
  std::vector<TString> cateNames = m_config->getStrV("WorkspaceCateNames");
  for (int i_c = 0; i_c < (int)cateNames.size(); i_c++) {
    pad1->cd();
    pad1->Clear();

    std::vector<TString> obsNames = m_config->getStrV("WorkspaceObsNames");
    TString obsName = obsNames[i_c];
    double obsMin = (*m_workspace->var(obsName)).getMin();
    double obsMax = (*m_workspace->var(obsName)).getMax();
    
    // Set the resonant analysis plot binning and axis scale to paper settings:
    int nBinsForPlot = (int)((obsMax - obsMin) / m_geVPerBin);
    
    // Plot everything on RooPlot:
    RooPlot* frame = (*m_workspace->var(obsName)).frame(nBinsForPlot);
    if ((m_workspace->pdf("model_"+cateNames[i_c]))) {
      (*m_workspace->data(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()))).plotOn(frame);
      (*m_workspace->pdf("model_"+cateNames[i_c])).plotOn(frame, LineColor(2), LineStyle(1));
    }
    frame->SetXTitle(obsName);
    frame->SetYTitle(Form("Events / %2.1f GeV", m_geVPerBin));
    frame->Draw();
    
    if (m_useLogScale) {
      gPad->SetLogy();
      frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    }
    else frame->GetYaxis()->SetRangeUser(m_yMin, m_yMax);
    
    TH1F *histDH = new TH1F("histDH", "histDH", 1, 0, 1);
    TH1F *histSH = new TH1F("histSH", "histSH", 1, 0, 1);
    TH1F *histBkg = new TH1F("histBkg", "histBkg", 1, 0, 1);
    TH1F *histSig = new TH1F("histSig", "histSig", 1, 0, 1);
    histDH->SetLineColor(6);
    histSH->SetLineColor(3);
    histBkg->SetLineColor(4);
    histSig->SetLineColor(2);
    
    histDH->SetLineStyle(3);
    histSH->SetLineStyle(4);
    histBkg->SetLineStyle(2);
    histSig->SetLineStyle(1);
    
    TLegend leg(0.63, 0.69, 0.92, 0.92);
    leg.SetFillColor(0);
    leg.SetTextSize(0.05);
    leg.SetBorderSize(0);
    leg.SetTextFont(42);
    leg.AddEntry(histDH, "DiHiggs", "l");
    leg.AddEntry(histSH, "Single Higgs", "l");
    leg.AddEntry(histBkg, "Continuum Bkg.", "l");
    leg.AddEntry(histSig, "Sum", "l");
    leg.Draw("SAME");
    
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
      = m_config->getStr(Form("PrintCateName_%s", cateNames[i_c].Data()));
    t.DrawLatex(0.20, 0.76, printCateName);
    
    // Second pad on the canvas (Ratio Plot):
    pad2->cd();
    pad2->Clear();
    
    double ratioMin = 0.0;//-0.2
    double ratioMax = 2.0;//2.2
    TH1F *medianHist = new TH1F("median","median",nBinsForPlot,obsMin,obsMax);
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
    medianHist->GetXaxis()->SetTitle(obsName);
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
      subData 
	= plotDivision(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()),
		       Form("model_%s", cateNames[i_c].Data()), obsName, obsMin,
		       obsMax, nBinsForPlot);
      
      TLine *line = new TLine();
      line->SetLineStyle(1);
      line->SetLineWidth(2);
      line->SetLineColor(kBlack);
      line->SetLineWidth(1);
      line->SetLineStyle(2);
      line->DrawLine(obsMin,((1.0+ratioMin)/2.0), obsMax,((1.0+ratioMin)/2.0));
      line->DrawLine(obsMin,((1.0+ratioMax)/2.0), obsMax,((1.0+ratioMax)/2.0));
      
      subData->Draw("EPSAME");
    }
    else {
      subData
	= plotSubtract(Form("%s_%s",datasetName.Data(),cateNames[i_c].Data()),
		       Form("model_%s", cateNames[i_c].Data()), obsName,
		       obsMin, obsMax, nBinsForPlot);
      medianHist->GetYaxis()->SetRangeUser(subData->GetYaxis()->GetXmin(),
					   subData->GetYaxis()->GetXmax());
      medianHist->Draw();
      subData->Draw("EPSAME");
    }
    
    can->Print(Form("%s/fitPlot_%s_%s_%s.eps", m_plotDir.Data(),
		    m_anaType.Data(), fitType.Data(), cateNames[i_c].Data()));
    can->Print(Form("%s/fitPlot_%s_%s_%s.C", m_plotDir.Data(),
		    m_anaType.Data(), fitType.Data(), cateNames[i_c].Data()));
    delete histDH;
    delete histSH;
    delete histBkg;
    delete histSig;
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
   Check whether the specified map entry exists.
    @param mapKey - The key for which we are finding a value:
    @return - True iff the categorization has been defined. 
*/
bool TestStat::mapValueExists(TString mapKey) {

  // Checks if there is a key corresponding to mapKey in the map: 
  bool nonExistent = (m_calculatedValues.find(mapKey) ==
		      m_calculatedValues.end());
  if (nonExistent) {
    std::cout << "TestStat: key " << mapKey << " not defined!" << std::endl;
  }
  return !nonExistent;
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
   Scan the NLL of a particular variable.
   @param scanName - The name of the scan (for saving plots...).
   @param datasetName - The name of the dataset in the workspace.
   @param varToScan - The variable that will be scanned.
   @param varsToFix - The variables that will be fixed at max-likelihood values.
   @return - The -1 sigma, median, and +1 sigma variable values.
*/
std::map<int,double> TestStat::scanNLL(TString scanName, TString datasetName,
				       TString varToScan,
				       std::vector<TString> varsToFix) {
  printer(Form("TestStat::scanNLL(datasetName = %s, varToScan = %s)",
	       datasetName.Data(), varToScan.Data()), false);
  
  // A map to store the scan results (median = 0, sigmas are +/-1, +/-2, etc...)
  std::map<int,double> result; result.clear();
  
  // Load objects for fitting:
  RooAbsPdf* combPdf = m_mc->GetPdf();
  RooArgSet* nuisanceParameters = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooArgSet* globalObservables = (RooArgSet*)m_mc->GetGlobalObservables();
  RooArgSet* origValNP = (RooArgSet*)m_workspace->getSnapshot("paramsOrigin");
  RooArgSet* poi = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar* firstPoI = (RooRealVar*)poi->first();
  RooArgSet* poiAndNuis = new RooArgSet();
  poiAndNuis->add(*nuisanceParameters);
  poiAndNuis->add(*poi);
  
  m_workspace->loadSnapshot("paramsOrigin");
  
  // Look for dataset:
  if (!m_workspace->data(datasetName)) {
    printer(Form("TestStat: Error! Requested data not available: %s",
		 datasetName.Data()), true);
  }
  
  // Allow the PoI to float:
  firstPoI->setConstant(false);
  
  // Release nuisance parameters before fit and set to the default values
  statistics::constSet(nuisanceParameters, false, origValNP);
  // The global observables should be fixed to the nominal values...
  statistics::constSet(globalObservables, true);
  
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
    
  // First do a fit to get the unconditional maximum likelihood fit result:
  RooNLLVar* varMLNLL
    = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				     Extended(combPdf->canBeExtended()));
  RooFitResult *fitResultML = statistics::minimize(varMLNLL, "", NULL, true);
  if (!fitResultML || fitResultML->status() != 0) m_allGoodFits = false;
  m_workspace->saveSnapshot("paramsProfileML", *poiAndNuis);
  double minNLLVal = varMLNLL->getVal();
  delete fitResultML;
  delete varMLNLL;
  
  // Then retrieve the central value and error as "preliminary result":
  double errorLo = m_workspace->var(varToScan)->getError();
  double errorHi = m_workspace->var(varToScan)->getError();
  result[0] = m_workspace->var(varToScan)->getVal();
  result[-1] = result[0] - errorLo;
  result[1] = result[0] + errorHi;
  double scanMin = result[0] - (2.0 * errorLo);
  double scanMax = result[0] + (2.0 * errorHi);
  if (scanMin < m_workspace->var(varToScan)->getMin()) {
    scanMin = m_workspace->var(varToScan)->getMin();
  }
  if (scanMax > m_workspace->var(varToScan)->getMax()) {
    scanMax = m_workspace->var(varToScan)->getMax();
  }
  
  // Compile a list of points to scan:
  std::vector<double> scanValues; scanValues.clear();
  for (double currValue = scanMin; currValue <= scanMax; 
       currValue += ((scanMax-scanMin)/m_config->getNum("NumberOfScanPoints"))){
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
    // Load the snapshot from the ML fit:
    m_workspace->loadSnapshot("paramsProfileML");
    
    // Release nuisance parameters before fit and set to the default values
    statistics::constSet(nuisanceParameters, false);
    // The global observables should be fixed to the nominal values...
    statistics::constSet(globalObservables, true);
    
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
        
    // Also fix other variables to ML values:
    for (int i_v = 0; i_v < (int)varsToFix.size(); i_v++) {
      if (m_workspace->var(varsToFix[i_v])) {
	m_workspace->var(varsToFix[i_v])->setConstant(true);
      }
      else {
	printer(Form("TestStat: ERROR! Parameter %s not found in workspace!",
		     varsToFix[i_v].Data()), true);
      }
    }
    
    // Set the value of the variable to scan:
    m_workspace->var(varToScan)->setVal(scanValues[scanIndex]);
    m_workspace->var(varToScan)->setConstant(true);
    
    // The actual fit command:
    RooNLLVar* currNLL
      = (RooNLLVar*)combPdf->createNLL(*m_workspace->data(datasetName),
				       Extended(combPdf->canBeExtended()));
    
    RooFitResult *currFitResult = statistics::minimize(currNLL, "", NULL, true);
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
  result[-1] = result[0] - graphIntercept(gNLLLo, 1.0);
  result[1] = graphIntercept(gNLLHi, 1.0) - result[0];

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
  
  // Return the median and +/-1 sigma values:
  return result;
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
