////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  ToyAnalysis.cxx                                                           //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/02/2016                                                          //
//                                                                            //
//  This program compares test statistic values from pseudo-experiment        //
//  ensembles and asymptotic formulae.                                        //
//                                                                            //
//  Options: ForcePlot                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ToyAnalysis.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the ToyAnalysis class.
   @param newConfigFile - The name of the analysis config file.
   @param options - Options for the toy analysis: ForcePlot, CLScan
*/
ToyAnalysis::ToyAnalysis(TString newConfigFile, TString options) {
  
  // Load the config file:
  m_config = new Config(newConfigFile);
  TString jobName = m_config->getStr("JobName");
  TString anaType = m_config->getStr("AnalysisType");
  
  printer(Form("ToyAnalysis::ToyAnalysis(%s)",newConfigFile.Data()),false);

  // Set output directory:
  setOutputDir(Form("%s/%s/ToyAnalysis", 
		    (m_config->getStr("MasterOutput")).Data(),
		    jobName.Data()));
  
  // Set the internal (private) variable initial conditions:
  m_nBins = 500;
  m_binMin = 0;
  m_binMax = 20;
  
  m_valuesQMu_Mu0.clear();
  m_valuesQMu_Mu1.clear();
  m_valuesBestFit_AsymZ0_Mu0.clear();
  m_valuesBestFit_AsymCL_Mu0.clear();
  m_valuesBestFit_AsymZ0_Mu1.clear();
  m_valuesBestFit_AsymCL_Mu1.clear();
  
  // Set the default fit types:
  m_fitTypes.clear();
  m_fitTypes.push_back("0"); 
  m_fitTypes.push_back("1");
  m_fitTypes.push_back("Free");
    
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
  
  // Load the workspace and model config from file:
  //TFile workspaceFile(m_config->getStr("WorkspaceFile"), "read");
  m_workspaceFile = new TFile(m_config->getStr("WorkspaceFile"), "read");
  //m_workspace
  //= (RooWorkspace*)workspaceFile.Get(m_config->getStr("WorkspaceName"));
  m_workspace
    = (RooWorkspace*)m_workspaceFile->Get(m_config->getStr("WorkspaceName"));
  if (m_workspace->obj(m_config->getStr("WorkspaceModelConfig"))) {
    m_mc = (ModelConfig*)m_workspace
      ->obj(m_config->getStr("WorkspaceModelConfig"));
  }
  
  // Create the test statistic class:
  m_ts = new TestStat(newConfigFile, "new", m_workspace);
  if (!m_workspace) m_filesLoaded = false;
  else m_filesLoaded = true;
}

/**
   -----------------------------------------------------------------------------
   Check whether the input TTrees were successfully loaded.
   @return - True iff. all input files are OK.
*/
bool ToyAnalysis::areInputFilesOK() {
  return m_filesLoaded;
}

/**
   -----------------------------------------------------------------------------
   Calculate the qMu corresponding to the median or +/-1,2 sigma values in
   S+B fits to background only toys.
   @param N - The sigma value (-2,-1,0,1,2). Use 0 for median.
   @return - The value of qMu.
*/
double ToyAnalysis::calculateBkgQMuForN(double N) {
  // First get the probability corresponding to the Gaussian N:
  double probability = getPbFromN(N);
  
  // Then find first q_mu value above this probability:
  std::sort(m_valuesQMu_Mu0.begin(), m_valuesQMu_Mu0.end());
  int totalEvents = (int)m_valuesQMu_Mu0.size();
  for (int i_e = 0; i_e < totalEvents; i_e++) {
    double fraction = (((double)i_e) / ((double)totalEvents));
    if (fraction >= probability) {
      printer(Form("calculateBkgQMuForN(%f) = %f using i_e=%d / totalEvents=%d",
		   N, m_valuesQMu_Mu0[i_e], i_e, totalEvents), false);
      return m_valuesQMu_Mu0[i_e];
    }
  }
  printer("ToyAnalysis::ERROR in QMu bkg calculation", true);
  return 0.0;
}

/**
   -----------------------------------------------------------------------------
   Calculate the CLs limit from toys (bands with N=-2,-1,0,+1,+2)
   @param qMu - The value of the test statistic.
   @return - The CLs value.
*/
double ToyAnalysis::calculateCLsFromToy(double qMu) {
  double pMu = calculatePMuFromToy(qMu);
  //double pB = getPbFromN(N);
  double pB = calculatePBFromToy(qMu);
  double CLs = pMu / (1.0 - pB);
  printer(Form("ToyAnalysis::calculateCLsFromToy(%f) = (%f / (1 - %f)) = %f",
	       qMu, pMu, pB, CLs), false);
  
  // Check that it is finite:
  if (std::isnan(CLs) || !std::isfinite(CLs)) {
    CLs = 1.0;
    printer("ToyAnalysis::converting nan or inf to 1",false);
  }
  return CLs;
}

/**
   -----------------------------------------------------------------------------
   Calculate the CL from the toys.
   @param qMu - The value of the test statistic.
   @return - The CL value.
*/
double ToyAnalysis::calculateCLFromToy(double qMu) {
  return (1.0 - calculateCLsFromToy(qMu));
}

/**
   -----------------------------------------------------------------------------
   Calculate the error from toy MC statistics on calculated p-value, using 
   binomial errors.
   @param pValue - The p-value for which we are computing the error.
   @param nToys - The number of toy MC used to compute the p-value.
   @return - The binomial error on the p-value.
*/
double ToyAnalysis::calculateErrorPVal(double pValue, int nToys) {
  return sqrt(pValue * (1.0 - pValue) / ((double)nToys));
}

/**
   -----------------------------------------------------------------------------
   Calculate the error from toy MC statistics on calculated CL value, using the
   binomial errors.
   @param qMu - The value of the test statistic.
   @return - The binomial error on the p-value.
*/
double ToyAnalysis::calculateErrorCLVal(double qMu) {
  // First obtain the p-values:
  double pMu = calculatePMuFromToy(qMu);
  double pB = calculatePBFromToy(qMu);
  
  // Then calculate the error for each p-value:
  double pBErr = calculateErrorPVal(pB, ((int)m_valuesQMu_Mu0.size()));
  double pMuErr = calculateErrorPVal(pMu, ((int)m_valuesQMu_Mu1.size()));
  
  // Finally, sum the errors in quadrature:
  return sqrt((pBErr * pBErr) + (pMuErr * pMuErr));
}

/**
   -----------------------------------------------------------------------------
   Calculate the error from toy MC statistics on calculated CL value, using 
   counting errors on the p-value.
   @param pValue - The p-value for which we are computing the error.
   @param nToys - The number of toy MC used to compute the p-value.
   @return - The number counting error on p.
*/
double ToyAnalysis::calculateErrorFromCounting(double pValue, int nToys) {
  double nToysHi = pValue * ((double)nToys);
  double errNumerator = sqrt(nToysHi) / nToysHi;
  double errDenominator = sqrt((double)nToys) / ((double)nToys);
  double fracErr = sqrt((errNumerator*errNumerator) + 
			(errDenominator*errDenominator));
  return fracErr * pValue;
}

/**
   -----------------------------------------------------------------------------
   Calculate the fractional number of toys in the mu=0 toy dataset that give a
   value of QMu less than than that observed QMu in data.
   @param qMu - The value of the test statistic.
   @return - The value of pB.
*/
double ToyAnalysis::calculatePBFromToy(double qMu) {
  int totalEvents = (int)m_valuesQMu_Mu0.size();
  int passingEvents = 0;
  for (int i_e = 0; i_e < totalEvents; i_e++) {
    if (m_valuesQMu_Mu0[i_e] <= qMu) passingEvents++;
  }
  double probability = (totalEvents > 0) ?
    (((double)passingEvents) / ((double)totalEvents)) : 0.0;
  
  printer(Form("calculatePBFromToy(%f) using passingEvents=%d / totalEvents=%d",
	       qMu, passingEvents, totalEvents), false);
  
  return probability;
}

/**
   -----------------------------------------------------------------------------
   Calculate the fractional number of toys in the mu=1 toy dataset that give a
   value of QMu larger than that observed QMu in data.
   @param qMu - The value of the test statistic.
   @return - The value of pMu.
*/
double ToyAnalysis::calculatePMuFromToy(double qMu) {
  int totalEvents = (int)m_valuesQMu_Mu1.size();
  int passingEvents = 0;
  for (int i_e = 0; i_e < totalEvents; i_e++) {
    if (m_valuesQMu_Mu1[i_e] >= qMu) passingEvents++;
  }
  double probability = (totalEvents > 0) ?
    (((double)passingEvents) / ((double)totalEvents)) : 0.0;
  
  printer(Form("calculatePBFromToy(%f) using passingEvents=%d / totalEvents=%d",
	       qMu, passingEvents, totalEvents), false);
  
  return probability;
}

/**
   -----------------------------------------------------------------------------
   Calculate pB based on the standard deviation.
   @param N - The standard deviation.
   @return - The value of pB.
*/
double ToyAnalysis::getPbFromN(double N) {
  //return (1 - ROOT::Math::gaussian_cdf(N));
  return ROOT::Math::gaussian_cdf(N);
}

/**
   -----------------------------------------------------------------------------
   Fill the histograms containing toy Data.
   @param muValue - the mu hypothesis under which the toys were generated.
   @param toyTree - the TTree containing the pseudo data.
*/
void ToyAnalysis::fillToyHistograms(int muValue, ToyTree *toyTree) {
  printer(Form("ToyAnalysis::fillToyHistograms(%d)",muValue), false);
  
  // Instantiate the histograms:
  m_hMuProfiled[muValue] = new TH1F(Form("hMuProfiled%d",muValue),
				    Form("hMuProfiled%d",muValue), 
				    m_nBins, -2.0, 4.0);
  m_hQMu[muValue] = new TH1F(Form("hQMu%d",muValue),Form("hQMu%d",muValue),
			     m_nBins, m_binMin, m_binMax);
  m_hQ0[muValue] = new TH1F(Form("hQ0%d",muValue),Form("hQ0%d",muValue),
			    m_nBins, m_binMin, m_binMax);
  m_hZ0[muValue] = new TH1F(Form("hZ0%d",muValue),Form("hZ0%d",muValue),
			    m_nBins, 0, 5);
  m_hCL[muValue] = new TH1F(Form("hCL%d",muValue),Form("hCL%d",muValue),
			    m_nBins, 0, 1);
  
  // Store names and numbers of parameters:
  m_namesGlobs.clear();
  m_namesNuis.clear();
  m_namesPoIs.clear();
  
  // Clear the qMu storage for pMu calculation if looking at Mu1 toy:
  if (muValue > 0) {
    m_valuesQMu_Mu1.clear();
    m_valuesBestFit_AsymZ0_Mu1.clear();
    m_valuesBestFit_AsymCL_Mu1.clear();
  }
  else  {
    m_valuesQMu_Mu0.clear();
    m_valuesBestFit_AsymZ0_Mu0.clear();
    m_valuesBestFit_AsymCL_Mu0.clear();
  }
  
  // Get names of NP:
  RooArgSet *setNuis = (RooArgSet*)m_mc->GetNuisanceParameters();
  RooRealVar *nuis = NULL;
  TIterator *iterNuis = setNuis->createIterator();
  while ((nuis = (RooRealVar*)iterNuis->Next())) {
    TString nameNuis = nuis->GetName();
    for (int i_f = 0; i_f < (int)m_fitTypes.size(); i_f++) {
      TString nuisKey = Form("%s_Mu%sFit_Mu%dData",
			     nameNuis.Data(), m_fitTypes[i_f].Data(), muValue);
      m_histStorage[nuisKey] = new TH1F(nuisKey, nuisKey, 100, -5, 5);
    }
    m_namesNuis.push_back(nameNuis);
  }
  
  // Get names of Globs:
  RooArgSet *setGlobs = (RooArgSet*)m_mc->GetGlobalObservables();
  RooRealVar *globs = NULL;
  TIterator *iterGlobs = setGlobs->createIterator();
  while ((globs = (RooRealVar*)iterGlobs->Next())) {
    TString nameGlobs = globs->GetName();
    for (int i_f = 0; i_f < (int)m_fitTypes.size(); i_f++) {
      TString globsKey = Form("%s_Mu%sFit_Mu%dData",
			      nameGlobs.Data(),m_fitTypes[i_f].Data(), muValue);
      m_histStorage[globsKey] = new TH1F(globsKey, globsKey, 100, -5, 5);
    }
    m_namesGlobs.push_back(nameGlobs);
  }
  
  // Also get names of non-systematic parameters:
  RooArgSet *setPars = (RooArgSet*)m_mc->GetParametersOfInterest();
  RooRealVar *pars = NULL;
  TIterator *iterPars = setPars->createIterator();
  while ((pars = (RooRealVar*)iterPars->Next())) {
    TString namePars = pars->GetName();
    
    if (m_config->isDefined(Form("PoIRange_%s",namePars.Data()))) {
      std::vector<double> currRange
	= m_config->getNumV(Form("PoIRange_%s",namePars.Data()));
      pars->setRange(currRange[0], currRange[1]);
    }
        
    for (int i_f = 0; i_f < (int)m_fitTypes.size(); i_f++) {
      TString parsKey = Form("%s_Mu%sFit_Mu%dData",
			     namePars.Data(), m_fitTypes[i_f].Data(), muValue);
      m_histStorage[parsKey]
	= new TH1F(parsKey, parsKey, 100, pars->getMin(), pars->getMax());
    }
    m_namesPoIs.push_back(namePars);
  }
  
  //----------------------------------------//
  // Loop over events in the TTree:
  int nEvents = toyTree->fChain->GetEntries();
  printer(Form("ToyAnalysis: Looping over %d toy events.",nEvents), false);
  bool isFirstLoop = true;
  for (int i_e = 0; i_e < nEvents; i_e++) {
    toyTree->fChain->GetEntry(i_e);
    
    // Only plot successful fits:
    //if (!(toyTree->convergedMu0 && toyTree->convergedMu1 &&
    //toyTree->convergedMuFree)) continue;
    
    // Get the test statistic values:
    double valueQMu = m_ts->getQMuFromNLL(toyTree->nllMu1, toyTree->nllMuFree,
					    toyTree->profiledPOIVal, 1);
    double valueQ0 = m_ts->getQ0FromNLL(toyTree->nllMu0, toyTree->nllMuFree,
					  toyTree->profiledPOIVal);
    double valueZ0 = m_ts->getZ0FromQ0(valueQ0);
    double valueCL = m_ts->getCLFromQMu(valueQMu, 0);
        
    // Fill histograms for the test statistics and POI:
    m_hQMu[muValue]->Fill(valueQMu);
    m_hQ0[muValue]->Fill(valueQ0);
    m_hMuProfiled[muValue]->Fill(toyTree->profiledPOIVal);
    m_hZ0[muValue]->Fill(valueZ0);
    m_hCL[muValue]->Fill(valueCL);
  
    // Also fill the QMu vector for pMu calculation:
    if (muValue > 0) {
      m_valuesQMu_Mu1.push_back(valueQMu);
      m_valuesBestFit_AsymZ0_Mu1.push_back(valueZ0);
      m_valuesBestFit_AsymCL_Mu1.push_back(valueCL);
    }
    else {
      m_valuesQMu_Mu0.push_back(valueQMu);
      m_valuesBestFit_AsymZ0_Mu0.push_back(valueZ0);
      m_valuesBestFit_AsymCL_Mu0.push_back(valueCL);
    }
    
    // Loop over the nuis:
    for (int i_n = 0; i_n < (int)((*toyTree->namesNP).size()); i_n++) {
      TString currNPName = (*toyTree->namesNP)[i_n];
      m_histStorage[Form("%s_Mu0Fit_Mu%dData", currNPName.Data(), muValue)]
	->Fill((*toyTree->valuesNPMu0)[i_n]);
      //m_histStorage[Form("%s_Mu1Fit_Mu%dData", currNPName.Data(), muValue)]
      //->Fill((*toyTree->valuesNPMu1)[i_n]);
      m_histStorage[Form("%s_MuFreeFit_Mu%dData", currNPName.Data(), muValue)]
	->Fill((*toyTree->valuesNPMuFree)[i_n]);
    }
    
    // Loop over the globs:
    for (int i_g = 0; i_g < (int)((*toyTree->namesGlobs).size()); i_g++) {
      TString currGlobName = (*toyTree->namesGlobs)[i_g];
      m_histStorage[Form("%s_Mu0Fit_Mu%dData", currGlobName.Data(), muValue)]
	->Fill((*toyTree->valuesGlobsMu0)[i_g]);
      //m_histStorage[Form("%s_Mu1Fit_Mu%dData", currGlobName.Data(), muValue)]
      //->Fill((*toyTree->valuesGlobsMu1)[i_g]);
      m_histStorage[Form("%s_MuFreeFit_Mu%dData", currGlobName.Data(), muValue)]
	->Fill((*toyTree->valuesGlobsMuFree)[i_g]);
    }
    
    // Loop over the non-systematic parameters:
    for (int i_p = 0; i_p < (int)((*toyTree->namesPoIs).size()); i_p++) {
      TString currPoIName = (*toyTree->namesPoIs)[i_p];
      m_histStorage[Form("%s_Mu0Fit_Mu%dData", currPoIName.Data(), muValue)]
	->Fill((*toyTree->valuesPoIsMu0)[i_p]);
      //m_histStorage[Form("%s_Mu1Fit_Mu%dData", currPoIName.Data(), muValue)]
      //->Fill((*toyTree->valuesPoIsMu1)[i_p]);
      m_histStorage[Form("%s_MuFreeFit_Mu%dData", currPoIName.Data(), muValue)]
	->Fill((*toyTree->valuesPoIsMuFree)[i_p]);
    }
  }
  
  // Then scale the statistics histograms:
  m_hQMu[muValue]->Scale(1.0 / m_hQMu[muValue]->Integral(1, m_nBins));
  m_hQ0[muValue]->Scale(1.0 / m_hQ0[muValue]->Integral(1, m_nBins));
  m_hZ0[muValue]->Scale(1.0 / m_hZ0[muValue]->Integral(1, m_nBins));
  m_hCL[muValue]->Scale(1.0 / m_hCL[muValue]->Integral(1, m_nBins));
  
}

/**
   -----------------------------------------------------------------------------
   Get the asymptotic form of the test statistic. Stored in the m_hAsymptotic
   histogram.
   @param statistic - the test statistic for asymptotic formula.

void ToyAnalysis::getAsymptoticForm(TString statistic) {
  std::cout << "ToyAnalysis: Constructing Asymptotic form of "
	    << statistic << "." << std::endl;
  
  // The histogram to contain the asymptotic form:
  m_hAsymptotic 
    = new TH1F("hAsymptotic", "hAsymptotic", m_nBins, m_binMin, m_binMax);
  
  // First get the value from fitting Asimov data (asimovDataMu0):
  double muHat = 0.0;
  double nllMu1 = m_ts->getFitNLL("asimovDataMu0", 1.0, true, muHat);
  double nllMu0 = m_ts->getFitNLL("asimovDataMu0", 0.0, true, muHat);
  double nllMuHat = m_ts->getFitNLL("asimovDataMu0", 0.0, false, muHat);
  double qMu = m_ts->getQMuFromNLL(nllMu1, nllMuHat, muHat, 1);
  double qMuTilde 
    = m_ts->getQMuTildeFromNLL(nllMu1, nllMu0, nllMuHat, muHat,1);
  
  double asimovTestStat = 0.0;
  if (statistic.EqualTo("QMu")) asimovTestStat = qMu;
  else if (statistic.EqualTo("QMuTilde")) asimovTestStat = qMuTilde;
  
  // Construct the test statistic function:
  for (int i_b = 1; i_b <= 1000; i_b++) {
    double q = m_hAsymptotic->GetBinCenter(i_b);
    if (statistic.EqualTo("QMu")) {
      m_hAsymptotic->SetBinContent(i_b, m_ts->functionQMu(q));
    }
    else if (statistic.EqualTo("QMuTilde")) {
      m_hAsymptotic
	->SetBinContent(i_b, m_ts->functionQMuTilde(q,asimovTestStat));
    }
  }
  // Then scale:
  m_hAsymptotic->SetBinContent(0, 0);
  m_hAsymptotic->SetBinContent(m_nBins+1, 0);
  m_hAsymptotic->Scale(1.0 / m_hAsymptotic->Integral(1, m_nBins));
  m_hAsymptotic->SetBinContent(1,1+m_hAsymptotic->GetBinContent(1));// WHY?
  m_hAsymptotic->Scale(1.0 / m_hAsymptotic->Integral(1,m_nBins));
}
*/

/**
   -----------------------------------------------------------------------------
   Retrieve the functional form of the asymptotic test statistic approximation.
*/
TH1F* ToyAnalysis::getAsymptoticHist() {
  return m_hAsymptotic;
}

/**
   -----------------------------------------------------------------------------
   Get a list of fits performed for each toy:
*/
std::vector<TString> ToyAnalysis::getFitTypes() {
  return m_fitTypes;
}

/**
   -----------------------------------------------------------------------------
   Retrieve the histogram of a particular nuisance parameter, global observable,
   or general parameter.
   @param paramName - the name of the global observable.
   @param fitType - the type of fit generating the distribution.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* ToyAnalysis::getHist(TString paramName, TString fitType, int toyMu) {
  TString mapKey = Form("%s_Mu%sFit_Mu%dData", paramName.Data(),
			fitType.Data(), toyMu);
  if (m_histStorage.count(mapKey) > 0) return m_histStorage[mapKey];
  else printer(Form("ToyAnalysis:: ERROR no hist %s",mapKey.Data()), true);
  return NULL;
}

/**
   -----------------------------------------------------------------------------
   Get the signal strength histogram.
   @param toyMu - the mu value used to generate the toy data that was fitted.
*/
TH1F* ToyAnalysis::getMuHist(int toyMu) {
  return m_hMuProfiled[toyMu];
}

/**
   -----------------------------------------------------------------------------
   Get a list of global observable names.
   @return - A vector of global observable names.
*/
std::vector<TString> ToyAnalysis::getNamesGlobalObservables() {
  return m_namesGlobs;
}

/**
   -----------------------------------------------------------------------------
   Get a list of nuisance parameter names.
   @return - A vector of nuisance parameter names.
*/
std::vector<TString> ToyAnalysis::getNamesNuisanceParameters() {
  return m_namesNuis;
}

/**
   -----------------------------------------------------------------------------
   Get a list of parameter of interest names.
   @return - A vector of parameter of interest names.
*/
std::vector<TString> ToyAnalysis::getNamesPoI() {
  return m_namesPoIs;
}

/**
   -----------------------------------------------------------------------------
   Get the test statistic histogram.
   @param statistic - The name of the test statistic: Q0, QMu, QMuTilde, Z0, CL.
   @param toyMu - The mu value used to generate the toy data that was fitted.
   @return - A histogram with the statistical values.
*/
TH1F* ToyAnalysis::getStatHist(TString statistic, int toyMu) {
  if (statistic.EqualTo("Q0")) return m_hQ0[toyMu];
  else if (statistic.EqualTo("QMu")) return m_hQMu[toyMu];
  else if (statistic.EqualTo("Z0")) return m_hZ0[toyMu];
  else if (statistic.EqualTo("CL")) return m_hCL[toyMu];
  //else if (statistic.EqualTo("QMuTilde")) return m_hQMuTilde[toyMu];
  else printer(Form("ToyAnalysis: ERROR! no %s", statistic.Data()), true);
  return NULL;
}

/**
   -----------------------------------------------------------------------------
   Get the vector of statistics values.
   @param statistic - The name of the test statistic: Q0, QMu, QMuTilde, Z0, CL.
   @param toyMu - The mu value used to generate the toy data that was fitted.
   @return - A vector of the given statistics values for the specified ensemble.
*/
std::vector<double> ToyAnalysis::getStatValues(TString statistic, int toyMu) {

  if (statistic.EqualTo("QMu")) {
    if (toyMu == 0) return m_valuesQMu_Mu0;
    else if (toyMu == 1) return m_valuesQMu_Mu1;
  }
  else if (statistic.EqualTo("Z0")) {
    if (toyMu == 0) return m_valuesBestFit_AsymZ0_Mu0;
    else if (toyMu == 1) return m_valuesBestFit_AsymZ0_Mu1;
  }
  else if (statistic.EqualTo("CL")) {
    if (toyMu == 0) return m_valuesBestFit_AsymCL_Mu0;
    else if (toyMu == 1) return m_valuesBestFit_AsymCL_Mu1;
  }
  std::cout << "ToyAnalysis: could not return stat value for statistic=" 
	    << statistic << " and toyMu=" << toyMu << std::endl;
  exit(0);
}

/**
   -----------------------------------------------------------------------------
   Load a toy ensemble from file. 
   @param toyMu - The mu value used to generate this toy dataset.
   @param toyFileForm - The form of the toy file names (incuding wildcarding).
*/

void ToyAnalysis::loadToy(int toyMu, TString toyFileForm) {
  
  // Add all of the individual pseudoexperiment files together:
  TString toyFile = Form("toy_mu%d.root", toyMu);
  system(Form("hadd -f %s %s", toyFile.Data(), toyFileForm.Data()));
  
  // Create temporary file lists:
  TString list = Form("temp_list_mu%d.txt", toyMu);
  system(Form("echo %s | tee %s", toyFile.Data(), list.Data()));
  
  // Create TChain of ntuples:
  TChain* chain = CommonFunc::MakeChain("toy", list, "badfile");
  ToyTree *tree = new ToyTree(chain);
  
  if (!(tree->fChain->GetEntries() > 0)) m_filesLoaded = false;
  
  // Store tree info in histograms:
  fillToyHistograms(toyMu, tree);
  
  // Remove the temporary file list:
  system(Form("rm %s", list.Data()));
  
  std::cout << "ToyAnalysis::loadToy(" << toyMu << ", " << toyFileForm 
	    << ") Finished!" << std::endl;
  std::cout << "\t" << tree->fChain->GetEntries()
	    << " mu=" << toyMu << " toy experiments were analyzed" << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Plot the distributions of nuisance parameters and global observables
   @param paramName - The name of the global observable.
   @param toyMu - The mu value used to generate the toy data that was fitted.
*/
void ToyAnalysis::plotHist(TString paramName, int toyMu) {
  printer(Form("ToyAnalysis::plotHist(%s, %d)",paramName.Data(),toyMu),false);
  
  // Retrieve the histograms for the parameter after various fits:
  TH1F *histMu0 = getHist(paramName, "0", toyMu);
  //TH1F *histMu1 = getHist(paramName, "1", toyMu);
  TH1F *histMuFree = getHist(paramName, "Free", toyMu);
  
  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
  gPad->SetLogy();
  
  // Format histograms:
  double min = 0.0001;
  double max = 10;
  histMu0->Scale(1.0 / histMu0->Integral());
  //histMu1->Scale(1.0 / histMu1->Integral());
  histMuFree->Scale(1.0 / histMuFree->Integral());
  histMu0->SetLineColor(kRed);
  //histMu1->SetLineColor(kGreen+2);
  histMuFree->SetLineColor(kBlue);
  histMu0->SetLineWidth(3);
  //histMu1->SetLineWidth(3);
  histMuFree->SetLineWidth(3);
  histMu0->SetLineStyle(1);
  //histMu1->SetLineStyle(2);
  histMuFree->SetLineStyle(4);
  
  // Format axis titles:
  histMu0->GetYaxis()->SetTitle("Fraction of toys");
  histMu0->GetXaxis()->SetTitle(paramName);
  histMu0->GetYaxis()->SetRangeUser(min,max);
  
  // Draw histograms:
  histMu0->Draw("hist");
  //histMu1->Draw("histSAME");
  histMuFree->Draw("histSAME");
  
  // Create a legend:
  TLegend leg(0.20, 0.80, 0.60, 0.92);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.03);
  leg.SetFillColor(0);
  leg.AddEntry(histMu0,Form("#mu=0 fixed, mean=%f",histMu0->GetMean()),"l");
  //leg.AddEntry(histMu1,Form("#mu=1 fixed, mean=%f",histMu1->GetMean()),"l");
  leg.AddEntry(histMuFree,Form("#hat{#mu}, mean=%f",histMuFree->GetMean()),"l");
  leg.Draw("SAME");

  // Print the canvas:
  can->Print(Form("%s/plot_%s_toy%i.eps", m_outputDir.Data(), paramName.Data(), 
		  toyMu));
  can->Clear();
  gPad->SetLogy(0);
}

/**
   -----------------------------------------------------------------------------
   Plot the profiled values of the signal strength.  
*/
void ToyAnalysis::plotProfiledMu() {
 
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  
  TH1F *histMu0 = getMuHist(0);
  TH1F *histMu1 = getMuHist(1);
  histMu0->SetLineColor(kBlue);
  histMu1->SetLineColor(kRed);
  histMu0->SetLineWidth(2);
  histMu1->SetLineWidth(2);

  gPad->SetLogy();
  histMu0->GetXaxis()->SetTitle("#mu_{profiled}");
  histMu0->GetYaxis()->SetTitle("Fraction of toys");
  histMu0->Draw("");
  histMu1->Draw("SAME");
  
  TLegend leg(0.56,0.74,0.88,0.86);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.04);
  leg.SetFillColor(0);
  leg.AddEntry(histMu0,"#mu=0 toy MC","l");
  leg.AddEntry(histMu1,"#mu=1 toy MC","l");
  leg.Draw("SAME");
  
  // Also get the mean signal strengths from the toy profiling:
  double meanMu0 = histMu0->GetMean();
  TLine *lineMu0 = new TLine();
  lineMu0->SetLineStyle(2);
  lineMu0->SetLineWidth(3);
  lineMu0->SetLineColor(kBlue); 
  lineMu0->DrawLine(meanMu0, histMu0->GetYaxis()->GetXmin(),
		    meanMu0, histMu0->GetYaxis()->GetXmax());
  double meanMu1 = histMu1->GetMean();
  TLine *lineMu1 = new TLine();
  lineMu1->SetLineStyle(2);
  lineMu1->SetLineWidth(3);
  lineMu1->SetLineColor(kRed);
  lineMu1->DrawLine(meanMu1, histMu0->GetYaxis()->GetXmin(),
		    meanMu1, histMu0->GetYaxis()->GetXmax());
  TLatex textMu0;
  textMu0.SetNDC();
  textMu0.SetTextColor(kBlue);
  textMu0.SetTextFont(42);
  textMu0.SetTextSize(0.04); 
  textMu0.DrawLatex(0.56, 0.65,Form("#mu=0 toy mean = %2.2f",meanMu0));
  TLatex textMu1;
  textMu1.SetNDC();
  textMu1.SetTextColor(kRed); 
  textMu1.SetTextFont(42);
  textMu1.SetTextSize(0.04);
  textMu1.DrawLatex(0.56, 0.7,Form("#mu=1 toy mean = %2.2f",meanMu1));
  
  can->Print(Form("%s/plot_profiledMu.eps", m_outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);  
}

/**
   -----------------------------------------------------------------------------
   Plot the distributions of nuisance parameters and global observables
   @param statistic -the name of the test statistic to plot.
*/
void ToyAnalysis::plotTestStat(TString statistic) {
  TH1F *hStatMu0 = getStatHist(statistic, 0);
  TH1F *hStatMu1 = getStatHist(statistic, 1);

  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
  
  hStatMu0->GetXaxis()->SetTitle(printStatName(statistic));
  hStatMu1->GetXaxis()->SetTitle(printStatName(statistic));
  
  hStatMu0->GetYaxis()->SetTitle("Normalized entries");  
  hStatMu1->GetYaxis()->SetTitle("Normalized entries");
  
  hStatMu0->SetLineColor(kBlue);
  hStatMu1->SetLineColor(kRed);
  
  // Draw test statistic histograms:
  gPad->SetLogy();
  hStatMu0->GetYaxis()->SetRangeUser(0.0001,1.0);
  hStatMu0->Draw("");
  hStatMu1->Draw("SAME");
  
  m_hAsymptotic->SetLineColor(kBlack);
  m_hAsymptotic->Draw("SAME");
  
  TLegend leg(0.49, 0.76, 0.84, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(hStatMu1, "#mu=1 toy MC","l");
  leg.AddEntry(hStatMu0, "#mu=0 toy MC","l");
  leg.AddEntry(m_hAsymptotic, "Asyptotic distribution","l");
  leg.Draw("SAME");
  
  can->Print(Form("%s/plot_%s.eps", m_outputDir.Data(), statistic.Data()));
  can->Clear();
  gPad->SetLogy(0);
}

/**
   -----------------------------------------------------------------------------
   Draw the asymptotic function along with toy test statistic, and comparisons.
   @param statistic - the name of the test statistic.
*/
void ToyAnalysis::plotTestStatComparison(TString statistic) {
  
  // Construct canvas and sub-pads:
  TCanvas *can = new TCanvas("can", "can", 1000, 1200);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.6, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.4, 1.0, 0.6);
  TPad *pad3 = new TPad("pad2", "pad2", 0.0, 0.2, 1.0, 0.4);
  TPad *pad4 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.2);
  pad1->SetBottomMargin(0.00001);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.00001);
  pad3->SetTopMargin(0.00001);
  pad3->SetBottomMargin(0.00001);
  pad4->SetTopMargin(0.00001);
  pad4->SetBottomMargin(0.4);
  pad1->SetBorderMode(0);
  pad2->SetBorderMode(0);
  pad3->SetBorderMode(0);
  pad4->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  
  // Get the toy and asymptotic distributions:
  TH1F *hStatMu1 = getStatHist(statistic, 1);
  //m_hAsymptotic;
  
  // Calculate the histograms to plot:
  TH1F *hIntegralToy = new TH1F("hIntegralToy", "hIntegralToy", 
				m_nBins, m_binMin, m_binMax);
  TH1F *hIntegralAsym = new TH1F("hIntegralAsym", "hIntegralAsym",
				 m_nBins, m_binMin, m_binMax);
  TH1F *hRatio = new TH1F("hRatio", "hRatio", m_nBins, m_binMin, m_binMax);
  TH1F *hSignificanceToy = new TH1F("hSignificanceToy", "hSignificanceToy",
				    m_nBins, m_binMin, m_binMax);
  TH1F *hSignificanceAsym = new TH1F("hSignificanceAsym", "hSignificanceAsym",
				     m_nBins, m_binMin, m_binMax);
  for (int i_b = 1; i_b <= m_nBins; i_b++) {
    double valueAsym = m_hAsymptotic->Integral(i_b, m_nBins);
    double valueToy = hStatMu1->Integral(i_b, m_nBins);
    hIntegralToy->SetBinContent(i_b, valueToy);
    hIntegralAsym->SetBinContent(i_b, valueAsym);
    hSignificanceToy->SetBinContent(i_b, -1.0*TMath::NormQuantile(valueToy));
    hSignificanceAsym->SetBinContent(i_b, -1.0*TMath::NormQuantile(valueAsym));
    
    double ratio = (TMath::NormQuantile(valueToy) == 0) ? 1.0 :
      TMath::NormQuantile(valueAsym) / TMath::NormQuantile(valueToy);
    
    hRatio->SetBinContent(i_b, ratio);
  }
   
  // Pad 1: Draw test statistic under mu=1 hypothesis and asymptotic function.
  pad1->cd();
  m_hAsymptotic->SetLineColor(kBlue);
  hStatMu1->Draw("");
  m_hAsymptotic->Draw("SAME");
  gPad->SetLogy();
  gPad->SetLogx();
  
  TLegend leg(0.65, 0.65, 0.9, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.06);
  TString printName = printStatName(statistic);
  leg.AddEntry(hStatMu1,Form("Toy MC %s",printName.Data()),"l");
  leg.AddEntry(m_hAsymptotic,Form("Asymptotic %s",printName.Data()),"l");
  leg.Draw("SAME");
  TLatex textToy;
  textToy.SetNDC();
  textToy.SetTextColor(kRed); 
  textToy.SetTextFont(42);
  textToy.SetTextSize(0.05);
  textToy.DrawLatex(0.2, 0.2, Form("Toy Mean = %2.2f", hStatMu1->GetMean()));
  TLatex textAsym;
  textAsym.SetNDC();
  textAsym.SetTextColor(kBlue); 
  textAsym.SetTextFont(42);
  textAsym.SetTextSize(0.05);
  textAsym.DrawLatex(0.2, 0.15, Form("Asymptotics Mean = %2.2f", 
				     m_hAsymptotic->GetMean()));
  
  // Pad 2: Draw the CDFs for asymptotics and pseudo-experiments
  pad2->cd();
  hIntegralToy->SetLineColor(kRed);
  hIntegralAsym->SetLineColor(kBlue);
  hIntegralToy->GetXaxis()->SetTitle(printName);
  hIntegralAsym->GetXaxis()->SetTitle(printName);
  hIntegralToy->GetYaxis()->SetTitle("CDF");
  hIntegralAsym->GetYaxis()->SetTitle("CDF");
  hIntegralToy->GetYaxis()->SetTitleSize(0.1);
  hIntegralToy->GetYaxis()->SetLabelSize(0.1);
  hIntegralToy->GetYaxis()->SetTitleOffset(0.65);
  hIntegralToy->Draw("");
  hIntegralAsym->Draw("SAME");
  gPad->SetLogy();
  gPad->SetLogx();
  
  // Pad 3: Plot the calculated significance corresponding to the CDF on Pad 2.
  pad3->cd();
  hSignificanceToy->SetLineWidth(2);
  hSignificanceAsym->SetLineWidth(2);
  hSignificanceToy->SetLineColor(kRed);
  hSignificanceAsym->SetLineColor(kBlue);
  hSignificanceToy->GetXaxis()->SetTitle(printName);
  hSignificanceAsym->GetXaxis()->SetTitle(printName);
  hSignificanceToy->GetYaxis()->SetTitle("Z [#sigma]");
  hSignificanceAsym->GetYaxis()->SetTitle("Z [#sigma]");
  hSignificanceToy->SetBinContent(1, 0);
  hSignificanceAsym->SetBinContent(1, 0);
  hSignificanceToy->GetYaxis()->SetTitleSize(0.1);
  hSignificanceToy->GetYaxis()->SetLabelSize(0.1);
  hSignificanceToy->GetYaxis()->SetTitleOffset(0.65);
  hSignificanceToy->GetYaxis()->SetNdivisions(6);
  hSignificanceToy->Draw("");
  hSignificanceAsym->Draw("SAME");
  gPad->SetLogx();
  
  // Pad 4: Calculate the ratio of Z values from Pad 3.
  pad4->cd();
  hRatio->SetLineColor(kBlack);
  hRatio->SetLineWidth(2);
  hRatio->GetYaxis()->SetTitle("Z_{Asym.} / Z_{Toy}");
  hRatio->GetXaxis()->SetTitle(printName);
  hRatio->GetYaxis()->SetTitleSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetTitleOffset(0.65);
  hRatio->GetYaxis()->SetNdivisions(6);
  hRatio->GetXaxis()->SetTitleSize(0.1);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->Draw();
  gPad->SetLogx();
  TLine *line2 = new TLine();
  line2->SetLineStyle(1);
  line2->SetLineWidth(1);
  line2->SetLineColor(kBlack);
  line2->DrawLine(m_binMin, 1.0, m_binMax, 1.0);
  line2->SetLineStyle(2);
  line2->DrawLine(m_binMin, 1.2, m_binMax, 1.2);
  line2->DrawLine(m_binMin, 0.8, m_binMax, 0.8);
  hRatio->Draw("SAME");
  
  can->Print(Form("%s/plot_comp_%s.eps", m_outputDir.Data(), statistic.Data()));
  can->Print(Form("%s/plot_comp_%s.png", m_outputDir.Data(), statistic.Data()));
  can->Clear();
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void ToyAnalysis::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Convert the statistic name into a LaTex formatted name.
   @param statistic - the name of the test statistic.
*/
TString ToyAnalysis::printStatName(TString statistic) {
  if (statistic.EqualTo("Q0")) return TString("q_{0}");
  else if (statistic.EqualTo("QMu")) return TString("q_{#mu}");
  else if (statistic.EqualTo("QMuTilde")) return TString("#tilde{q}_{#mu}");
  else return TString("q");
}

/**
   -----------------------------------------------------------------------------
   Set the fit types that are used in the analysis.
   @param fitTypes - A vector of fit types to be analyzed for each toy.
*/
void ToyAnalysis::setFitTypes(std::vector<TString> fitTypes) {
  m_fitTypes.clear();
  m_fitTypes = fitTypes;
}

/**
   -----------------------------------------------------------------------------
   Change the output directory for this class.
   @param outputDirectory - The new output directory.
*/
void ToyAnalysis::setOutputDir(TString outputDirectory) {
  // Set output directory:
  m_outputDir = outputDirectory;
  // Create output directory:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
}
