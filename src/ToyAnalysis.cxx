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
//  Options: StudyRetries                                                     //
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
  m_options = options;
  m_config = new Config(newConfigFile);
  TString jobName = m_config->getStr("JobName");
  TString anaType = m_config->getStr("AnalysisType");
  
  printer(Form("ToyAnalysis::ToyAnalysis(%s, %s)", 
	       newConfigFile.Data(), options.Data()), false);
  
  // Set output directory:
  setOutputDir(Form("%s/%s/ToyAnalysis", 
		    (m_config->getStr("MasterOutput")).Data(),
		    jobName.Data()));
  
  // Set the internal (private) variable initial conditions:
  m_nBins = 500;
  m_binMin = 0;
  m_binMax = 20;
  
  m_hAsymptoticQ0 = NULL;
  m_hAsymptoticQMu = NULL;

  m_weightsIS_Mu0.clear();
  m_weightsIS_Mu1.clear();
  m_valuesQMu_Mu0.clear();
  m_valuesQMu_Mu1.clear();
  m_valuesMuHat_Mu0.clear();
  m_valuesMuHat_Mu1.clear();
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
  m_workspaceFile = new TFile(m_config->getStr("WorkspaceFile"), "read");
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
  
  // Then sort the qMu vector (and associated weight vector):
  //std::sort(m_valuesQMu_Mu0.begin(), m_valuesQMu_Mu0.end());
  sortPairedVectors(m_valuesQMu_Mu0, m_weightsIS_Mu0);
  
  // Then find the first q_mu value above this probability:
  double denominator = 0.0;
  for (int i_e = 0; i_e < (int)m_valuesQMu_Mu0.size(); i_e++) {
    denominator += m_weightsIS_Mu0[i_e];
  }
  double numerator = 0.0;
  for (int i_e = 0; i_e < (int)m_valuesQMu_Mu0.size(); i_e++) {
    numerator += m_weightsIS_Mu0[i_e];
    double fraction = numerator / denominator;
    if (fraction >= probability) {
      printer(Form("calculateBkgQMuForN(%f) = %f using (%f pass / %f total)",
		   N, m_valuesQMu_Mu0[i_e], numerator, denominator), false);
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
  double totalEvents = 0.0;
  double passingEvents = 0.0;
  for (int i_e = 0; i_e < (int)m_valuesQMu_Mu0.size(); i_e++) {
    if (m_valuesQMu_Mu0[i_e] <= qMu) passingEvents += m_weightsIS_Mu0[i_e];
    totalEvents += m_weightsIS_Mu0[i_e];
  }
  double probability = (totalEvents > 0.0) ? (passingEvents/totalEvents) : 0.0;
  
  printer(Form("calculatePBFromToy(%f) using (%f passing / %f total)",
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
  double totalEvents = 0.0;
  double passingEvents = 0.0;
  for (int i_e = 0; i_e < (int)m_valuesQMu_Mu1.size(); i_e++) {
    if (m_valuesQMu_Mu1[i_e] >= qMu) passingEvents += m_weightsIS_Mu1[i_e];
    totalEvents += m_weightsIS_Mu1[i_e];
  }
  double probability = (totalEvents > 0.0) ? (passingEvents/totalEvents) : 0.0;
  
  printer(Form("calculatePBFromToy(%f) using (%f passing / %f total)",
	       qMu, passingEvents, totalEvents), false);
  return probability;
}

/**
   -----------------------------------------------------------------------------
   Check whether a particular type of fit is done in these toys.
   @param fitType - The type of fit to do (or not do).
   @return - True iff. the fit is done.
*/
bool ToyAnalysis::doThisFit(TString fitType) {
  for (int i_f = 0; i_f < (int)m_fitTypes.size(); i_f++) {
    if ((m_fitTypes[i_f]).EqualTo(fitType)) return true;
  }
  return false;
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
				    m_nBins, 0.0, 4.0);
  m_hQMu[muValue] = new TH1F(Form("hQMu%d",muValue),Form("hQMu%d",muValue),
			     m_nBins, m_binMin, m_binMax);
  m_hQ0[muValue] = new TH1F(Form("hQ0%d",muValue),Form("hQ0%d",muValue),
			    m_nBins, m_binMin, m_binMax);
  m_hZ0[muValue] = new TH1F(Form("hZ0%d",muValue),Form("hZ0%d",muValue),
			    m_nBins, 0, 5);
  m_hCL[muValue] = new TH1F(Form("hCL%d",muValue),Form("hCL%d",muValue),
			    m_nBins, 0, 1);
  
  // For investigations into fit retries:
  if (m_options.Contains("StudyRetries")) {
    m_hRetries[muValue] = new TH1F(Form("hRetries%d",muValue),
				   Form("hRetries%d",muValue), 50, 0, 50);
    m_hImprovement[muValue] = new TH1F(Form("hImprove%d",muValue),
				       Form("hImprove%d",muValue), 50, 0, 50);
    m_hMedImprovement[muValue] = new TH1F(Form("hMedImprove%d",muValue),
					  Form("hMedImprove%d",muValue),
					  50, 0, 50);
    m_hCounter[muValue] = new TH1F(Form("hCounter%d",muValue),
				   Form("hCounter%d",muValue), 50, 0, 50);
    
    m_h2RetriesZ[muValue] = new TH2F(Form("h2RetriesZ%d",muValue),
				     Form("h2RetriesZ%d",muValue),
				     50, 0, 50, 50, 0, 5);
    m_h2ZImprovement[muValue] = new TH2F(Form("h2ZImprove%d",muValue),
					 Form("h2ZImprove%d",muValue),
					 50, 0, 5, 50, 0.001, 5);
    m_hImpAtThisStep[muValue] = new TH1F(Form("hImpAtThisStep%d",muValue),
					 Form("hImpAtThisStep%d",muValue),
					 50, 0, 50);
    
    // Plot the Gaussian Z0 as a function of # retries:
    for (int i_r = 0; i_r < 50; i_r++) {
      m_hZ0Retries[muValue][i_r] = new TH1F(Form("hZ0%d_%d", muValue, i_r), 
					    Form("hZ0%d_%d", muValue, i_r), 
					    m_nBins, 0, 5);
    }
  }
  
  // An NLL map to keep track of the NLL per ML fit-retry:
  std::map<int,std::vector<double> > nllMap; nllMap.clear();
  for (int i_b = 1; i_b <= 50; i_b++) (nllMap[i_b]).clear();
  
  // Store names and numbers of parameters:
  m_namesGlobs.clear();
  m_namesNuis.clear();
  m_namesPoIs.clear();
  
  // Clear the qMu storage for pMu calculation if looking at Mu1 toy:
  if (muValue > 0) {
    m_weightsIS_Mu1.clear();
    m_valuesMuHat_Mu1.clear();
    m_valuesQMu_Mu1.clear();
    m_valuesBestFit_AsymZ0_Mu1.clear();
    m_valuesBestFit_AsymCL_Mu1.clear();
  }
  else  {
    m_weightsIS_Mu0.clear();
    m_valuesMuHat_Mu0.clear();
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
  
  // Also get names of parameters of interest:
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
    if (!(toyTree->convergedMu0 && toyTree->convergedMuFree)) continue;
        
    // Get the test statistic values:
    double valueQMu = m_ts->getQMuFromNLL(toyTree->nllMu1, toyTree->nllMuFree,
					    toyTree->profiledPOIVal, 1);
    double valueQ0 = m_ts->getQ0FromNLL(toyTree->nllMu0, toyTree->nllMuFree,
					  toyTree->profiledPOIVal);
    double valueZ0 = m_ts->getZ0FromQ0(valueQ0);
    double valueCL = m_ts->getCLFromQMu(valueQMu, 0);
    
    //double toyWeight = 1.0;
    double toyWeight = toyTree->weight;
    
    // Fill histograms for the test statistics and POI:
    m_hQMu[muValue]->Fill(valueQMu, toyWeight);
    m_hQ0[muValue]->Fill(valueQ0, toyWeight);
    m_hMuProfiled[muValue]->Fill(toyTree->profiledPOIVal, toyWeight);
    m_hZ0[muValue]->Fill(valueZ0, toyWeight);
    m_hCL[muValue]->Fill(valueCL, toyWeight);
    
    // Plots to make only in the case that retries are being studied:
    if (m_options.Contains("StudyRetries")) {
      m_hRetries[muValue]->Fill(toyTree->bestFitUpdate);
      m_h2RetriesZ[muValue]->Fill(toyTree->bestFitUpdate, valueZ0);
      
      // See how big the change is as a function of:
      double currMinimum=0.0; double prevMinimum=0.0; double improvement=0.0;
      for (int i_r = 0; i_r < (int)((*toyTree->nllPerRetry).size()); i_r++) { 
	
	if (i_r == 0) {
	  currMinimum = (*toyTree->nllPerRetry)[i_r];
	}
	else if ((*toyTree->nllPerRetry)[i_r] < currMinimum) {
	  prevMinimum = currMinimum;
	  currMinimum = (*toyTree->nllPerRetry)[i_r];
	  improvement = fabs(prevMinimum-currMinimum);
	  m_hImpAtThisStep[muValue]->Fill(i_r, improvement);
	  
	  if (i_r == toyTree->bestFitUpdate) {
	    m_hImprovement[muValue]->Fill(i_r, improvement);
	    (nllMap[i_r]).push_back(improvement);
	  }
	}
      }
      
      m_h2ZImprovement[muValue]->Fill(valueZ0, improvement);
      
      // For each event, find the best significance below a certain point:
      for (int i_r = 0; i_r < 50; i_r++) {
	
	// Loop over the stored NLLs and find the 
	double currMinNll = 0.0;
	for (int i_b = 0; i_b <= i_r; i_b++) {
	  if (currMinNll > (*toyTree->nllPerRetry)[i_b]) {
	    currMinNll = (*toyTree->nllPerRetry)[i_b];
	  }
	}
	
	double currValueQ0 = m_ts->getQ0FromNLL(toyTree->nllMu0, currMinNll,
						toyTree->profiledPOIVal);
	double currValueZ0 = m_ts->getZ0FromQ0(currValueQ0);
	m_hZ0Retries[muValue][i_r]->Fill(currValueZ0);
      }
    }// end of plots for retry study
    
    // Also fill the QMu vector for pMu calculation:
    if (muValue > 0) {
      m_weightsIS_Mu1.push_back(toyWeight);
      m_valuesMuHat_Mu1.push_back(toyTree->profiledPOIVal);
      m_valuesQMu_Mu1.push_back(valueQMu);
      m_valuesBestFit_AsymZ0_Mu1.push_back(valueZ0);
      m_valuesBestFit_AsymCL_Mu1.push_back(valueCL);
    }
    else {
      m_weightsIS_Mu0.push_back(toyWeight);
      m_valuesMuHat_Mu0.push_back(toyTree->profiledPOIVal);
      m_valuesQMu_Mu0.push_back(valueQMu);
      m_valuesBestFit_AsymZ0_Mu0.push_back(valueZ0);
      m_valuesBestFit_AsymCL_Mu0.push_back(valueCL);
    }
    
    // Loop over the nuis:
    for (int i_n = 0; i_n < (int)((*toyTree->namesNP).size()); i_n++) {
      TString currNPName = (*toyTree->namesNP)[i_n];
      if (doThisFit("0")) {
	m_histStorage[Form("%s_Mu0Fit_Mu%dData", currNPName.Data(), muValue)]
	  ->Fill((*toyTree->valuesNPMu0)[i_n], toyWeight);
      }
      if (doThisFit("1")) {
	m_histStorage[Form("%s_Mu1Fit_Mu%dData", currNPName.Data(), muValue)]
	  ->Fill((*toyTree->valuesNPMu1)[i_n], toyWeight);
      }
      if (doThisFit("Free")) {
	m_histStorage[Form("%s_MuFreeFit_Mu%dData", currNPName.Data(), muValue)]
	  ->Fill((*toyTree->valuesNPMuFree)[i_n], toyWeight);
      }
    }
    
    // Loop over the globs:
    for (int i_g = 0; i_g < (int)((*toyTree->namesGlobs).size()); i_g++) {
      TString currGlobName = (*toyTree->namesGlobs)[i_g];
      if (doThisFit("0")) {
	m_histStorage[Form("%s_Mu0Fit_Mu%dData", currGlobName.Data(), muValue)]
	  ->Fill((*toyTree->valuesGlobsMu0)[i_g], toyWeight);
      }
      if (doThisFit("1")) {
	m_histStorage[Form("%s_Mu1Fit_Mu%dData", currGlobName.Data(), muValue)]
	  ->Fill((*toyTree->valuesGlobsMu1)[i_g], toyWeight);
      }
      if (doThisFit("Free")) {
	m_histStorage[Form("%s_MuFreeFit_Mu%dData",currGlobName.Data(),muValue)]
	  ->Fill((*toyTree->valuesGlobsMuFree)[i_g], toyWeight);
      }
    }
    
    // Loop over the non-systematic parameters:
    for (int i_p = 0; i_p < (int)((*toyTree->namesPoIs).size()); i_p++) {
      TString currPoIName = (*toyTree->namesPoIs)[i_p];
      if (doThisFit("0")) {
	m_histStorage[Form("%s_Mu0Fit_Mu%dData", currPoIName.Data(), muValue)]
	  ->Fill((*toyTree->valuesPoIsMu0)[i_p], toyWeight);
      }
      if (doThisFit("1")) {
	m_histStorage[Form("%s_Mu1Fit_Mu%dData", currPoIName.Data(), muValue)]
	  ->Fill((*toyTree->valuesPoIsMu1)[i_p], toyWeight);
      }
      if (doThisFit("Free")) {
	m_histStorage[Form("%s_MuFreeFit_Mu%dData",currPoIName.Data(),muValue)]
	  ->Fill((*toyTree->valuesPoIsMuFree)[i_p], toyWeight);
      }
    }
  }
  
  // Then scale the statistics histograms:
  m_hMuProfiled[muValue]
    ->Scale(1.0 / m_hMuProfiled[muValue]->Integral(1, m_nBins));
  m_hQMu[muValue]->Scale(1.0 / m_hQMu[muValue]->Integral(1, m_nBins));
  m_hQ0[muValue]->Scale(1.0 / m_hQ0[muValue]->Integral(1, m_nBins));
  m_hZ0[muValue]->Scale(1.0 / m_hZ0[muValue]->Integral(1, m_nBins));
  m_hCL[muValue]->Scale(1.0 / m_hCL[muValue]->Integral(1, m_nBins));
  
  // Finish some calculations for the retry study:
  if (m_options.Contains("StudyRetries")) {

    for (int i_r = 0; i_r < 50; i_r++) {
      m_hZ0Retries[muValue][i_r]
	->Scale(1.0 / m_hZ0Retries[muValue][i_r]->Integral());
    }
    
    // Also calculate median improvement:
    for (int i_b = 1; i_b <= 50; i_b++) {
      int size = (int)((nllMap[i_b]).size());
      if (size > 0) {
	double medianVal = (nllMap[i_b][(int)((double)size/2.0)]);
	m_hMedImprovement[muValue]->SetBinContent(i_b, medianVal);
	m_hMedImprovement[muValue]
	  ->SetBinError(i_b, (medianVal * sqrt((double)size)/((double)size)));
      }
    }
    
    // And calculate the averages:
    for (int i_b = 1; i_b <= 50; i_b++) {
      m_hImpAtThisStep[muValue]
	->SetBinContent(i_b, (m_hImpAtThisStep[muValue]->GetBinContent(i_b) /
			      m_hRetries[muValue]->GetBinContent(i_b)));
      
      m_hImprovement[muValue]
	->SetBinContent(i_b, (m_hImprovement[muValue]->GetBinContent(i_b) /
			      m_hRetries[muValue]->GetBinContent(i_b)));
    }
  }
  
}

/**
   -----------------------------------------------------------------------------
   Get the asymptotic test statistic distribution. Stored in the 
   m_hAsymptoticQ0 and m_hAsymptoticQMu histograms.
   @param statistic - the test statistic for asymptotic formula.
*/
TH1F* ToyAnalysis::getAsymptoticHist(TString statistic) {
  printer(Form("ToyAnalysis: Constructing Asymptotic form of %s", 
	       statistic.Data()), false);
  
  // If asymptotic histogram has already been made, return it:
  if (statistic.EqualTo("QMu") && m_hAsymptoticQMu) return m_hAsymptoticQMu;
  else if (statistic.EqualTo("Q0") && m_hAsymptoticQ0) return m_hAsymptoticQ0;
  
  // Otherwise, create a new asymptotic histogram:
  TH1F *hAsymptotic = new TH1F(Form("hAsymptotic%s",statistic.Data()),
			       Form("hAsymptotic%s",statistic.Data()),
			       m_nBins, m_binMin, m_binMax);
  
  // First get the value from fitting Asimov data (asimovDataMu0):
  //double muHat = 0.0;
  //double nllMu1 = m_ts->getFitNLL("asimovDataMu0", 1.0, true, muHat);
  //double nllMu0 = m_ts->getFitNLL("asimovDataMu0", 0.0, true, muHat);
  //double nllMuHat = m_ts->getFitNLL("asimovDataMu0", 0.0, false, muHat);
  //double qMu = m_ts->getQMuFromNLL(nllMu1, nllMuHat, muHat, 1);
  //double qMuTilde 
  //= m_ts->getQMuTildeFromNLL(nllMu1, nllMu0, nllMuHat, muHat,1);
  
  //double asimovTestStat = 0.0;
  //if (statistic.EqualTo("QMu")) asimovTestStat = qMu;
  //else if (statistic.EqualTo("QMuTilde")) asimovTestStat = qMuTilde;
  
  // Construct the test statistic function:
  for (int i_b = 1; i_b <= m_nBins; i_b++) {
    double q = hAsymptotic->GetBinCenter(i_b);
    if (statistic.EqualTo("QMu")){
      hAsymptotic->SetBinContent(i_b, m_ts->functionQMu(q));
    }
    else if (statistic.EqualTo("Q0")) {
      hAsymptotic->SetBinContent(i_b, m_ts->functionQ0(q));
    }
    //else if (statistic.EqualTo("QMuTilde")) {
    //hAsymptotic
    //	->SetBinContent(i_b, m_ts->functionQMuTilde(q,asimovTestStat));
    //}
  }
  
  // Scale so that shape is 1/2 chi^2 and 1/2 delta at 0:
  hAsymptotic->SetBinContent(0, 0);
  hAsymptotic->SetBinContent(m_nBins+1, 0);
  hAsymptotic->Scale(1.0 / hAsymptotic->Integral(1, m_nBins));
  hAsymptotic->SetBinContent(1, 1 + hAsymptotic->GetBinContent(1));
  hAsymptotic->Scale(1.0 / hAsymptotic->Integral(1, m_nBins));
  
  // Then return a pointer to the desired asymptotic histogram:
  if (statistic.EqualTo("QMu")) {
    m_hAsymptoticQMu = hAsymptotic;
    return m_hAsymptoticQMu;
  }
  else if (statistic.EqualTo("Q0")) {
    m_hAsymptoticQ0 = hAsymptotic;
    return m_hAsymptoticQ0;
  }
  else {
    printer(Form("ToyAnalysis: Error returning asymptotic form of %s",
		 statistic.Data()), true);
  }
  return NULL;
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
  else if (statistic.EqualTo("ImportanceWeight")) {
    if (toyMu == 0) return m_weightsIS_Mu0;
    else if (toyMu == 1) return m_weightsIS_Mu1;
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
  TH1F *histMu0 = doThisFit("0") ? getHist(paramName,"0",toyMu) : NULL;
  TH1F *histMu1 = doThisFit("1") ? getHist(paramName,"1",toyMu) : NULL;
  TH1F *histMuFree = doThisFit("Free") ? getHist(paramName,"Free",toyMu) : NULL;
  
  TCanvas *can = new TCanvas("can", "can",800, 800);
  can->cd();
  gPad->SetLogy();

  // Format histograms:
  double min = 0.0001;
  double max = 10;
  if (doThisFit("0")) {
    histMu0->Scale(1.0 / histMu0->Integral());
    histMu0->SetLineColor(kRed);
    histMu0->SetLineWidth(3);
    histMu0->SetLineStyle(1);
    histMu0->GetYaxis()->SetTitle("Fraction of toys");
    histMu0->GetXaxis()->SetTitle(paramName);
    histMu0->GetYaxis()->SetRangeUser(min,max);
    histMu0->Draw("hist");
  }
  if (doThisFit("1")) {  
    histMu1->Scale(1.0 / histMu1->Integral());
    histMu1->SetLineColor(kGreen+2);
    histMu1->SetLineWidth(3);
    histMu1->SetLineStyle(4);
    histMu1->Draw("histSAME");
  }
  if (doThisFit("Free")) {
    histMuFree->Scale(1.0 / histMuFree->Integral());
    histMuFree->SetLineColor(kBlue);
    histMuFree->SetLineWidth(3);
    histMuFree->SetLineStyle(2);
    histMuFree->Draw("histSAME");
  }
    
  // Create a legend:
  TLegend leg(0.20, 0.80, 0.60, 0.92);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.03);
  leg.SetFillColor(0);
  leg.AddEntry(histMu0, Form("#mu=0 fixed, mean=%f", histMu0->GetMean()), "l");
  if (doThisFit("1")) {
    leg.AddEntry(histMu1, 
		 Form("#mu=1 fixed, mean=%f", histMu1->GetMean()), "l");
  }
  if (doThisFit("Free")) {
    leg.AddEntry(histMuFree,
		 Form("#hat{#mu}, mean=%f",histMuFree->GetMean()),"l");
  }
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
 
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  TH1F *histMu0 = getMuHist(0);
  TH1F *histMu1 = getMuHist(1);
  histMu0->SetLineColor(kBlue+1);
  histMu0->SetFillStyle(3354);
  histMu0->SetFillColor(kBlue);
  histMu0->SetLineWidth(2);
  
  histMu1->SetFillColor(kRed);
  histMu1->SetFillStyle(3345);
  histMu1->SetLineColor(kRed+1);
  histMu1->SetLineWidth(2);
  
  double yMinimum = 0.51 / ((double)histMu0->GetEntries());
  double yMaximum = 0.99;
  
  double xWidth = ((histMu1->GetXaxis()->GetXmax() - 
		    histMu1->GetXaxis()->GetXmin()) /
		   ((double)histMu1->GetNbinsX()));
  gPad->SetLogy();
  histMu0->GetXaxis()->SetTitle("#mu_{profiled}");
  histMu0->GetYaxis()->SetTitle(Form("Fraction per %2.3f", xWidth));
  
  histMu0->GetYaxis()->SetRangeUser(yMinimum, yMaximum);
  histMu0->Draw("");
  histMu1->Draw("SAME");
    
  // Also get the mean signal strengths from the toy profiling:
  double medMu0 =m_valuesMuHat_Mu0[(int)(0.5*(double)m_valuesMuHat_Mu0.size())];
  TLine *lineMu0 = new TLine();
  lineMu0->SetLineStyle(2);
  lineMu0->SetLineWidth(3);
  lineMu0->SetLineColor(kBlue+1); 
  lineMu0->DrawLine(medMu0, yMinimum, medMu0, yMaximum);
  double medMu1 =m_valuesMuHat_Mu1[(int)(0.5*(double)m_valuesMuHat_Mu1.size())];
  TLine *lineMu1 = new TLine();
  lineMu1->SetLineStyle(2);
  lineMu1->SetLineWidth(3);
  lineMu1->SetLineColor(kRed+1);
  lineMu1->DrawLine(medMu1, yMinimum, medMu1, yMaximum);
  
  // Create a legend of histograms and median lines:
  TLegend leg(0.56, 0.50, 0.88, 0.80);
  leg.SetBorderSize(0);
  leg.SetTextSize(0.05);
  leg.SetFillColor(0);
  leg.SetTextFont(42);
  leg.AddEntry(histMu0, "#mu=0 toy MC", "F");
  leg.AddEntry(lineMu0, Form("#mu=0 median = %2.2f",medMu0), "l");
  leg.AddEntry(histMu1, "#mu=1 toy MC", "F");
  leg.AddEntry(lineMu1, Form("#mu=1 median = %2.2f",medMu1), "l");
  leg.Draw("SAME");
  
  // Draw the ATLAS text:
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.58, 0.88, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.71, 0.88, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.58, 0.82, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (m_config->getNum("AnalysisLuminosity")/
				1000.0)));

  can->Print(Form("%s/plot_profiledMu.eps", m_outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);  
}

/**
   -----------------------------------------------------------------------------
   Plot the number of retries before the minimum NLL was found.
   @param muValue - The mu-value of the histogram.
*/
void ToyAnalysis::plotRetries(int muValue) {
  
  TCanvas *can = new TCanvas("can", "can", 800, 600);
  can->cd();
  
  m_hRetries[muValue]->SetLineColor(kRed+2);
  m_hRetries[muValue]->SetFillColor(kRed-10);
  m_hRetries[muValue]->SetLineWidth(2);
  m_hRetries[muValue]
    ->GetXaxis()->SetTitle("Retry with minimum -log(#it{L})");
  m_hRetries[muValue]->GetYaxis()->SetTitle("Number of toy MC");
  //gPad->SetLogy();
  m_hRetries[muValue]->Draw("hist");
  can->Print(Form("%s/plot_Retries.eps", m_outputDir.Data()));
  can->Clear();
  //gPad->SetLogy(0);  
  
  m_hImprovement[muValue]->SetLineColor(kRed+2);
  m_hImprovement[muValue]->SetFillColor(kRed-10);
  m_hImprovement[muValue]->SetLineWidth(2);
  m_hImprovement[muValue]
    ->GetXaxis()->SetTitle("Retry with minimum -log(#it{L})");
  m_hImprovement[muValue]->GetYaxis()->SetTitle("Mean -#Delta log(#it{L})");
  gPad->SetLogy();
  m_hImprovement[muValue]->Draw("hist");
  can->Print(Form("%s/plot_Improvement.eps", m_outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);
  
  m_hMedImprovement[muValue]->SetLineColor(kRed+2);
  m_hMedImprovement[muValue]->SetFillColor(kRed-10);
  m_hMedImprovement[muValue]->SetLineWidth(2);
  m_hMedImprovement[muValue]
    ->GetXaxis()->SetTitle("Retry with minimum -log(#it{L})");
  m_hMedImprovement[muValue]
    ->GetYaxis()->SetTitle("Median -#Delta log(#it{L})");
  gPad->SetLogy();
  m_hMedImprovement[muValue]->Draw("hist");
  //m_hMedImprovement[muValue]->Draw("E1");
  can->Print(Form("%s/plot_MedImprovement.eps", m_outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);

  
  m_hImpAtThisStep[muValue]->SetLineColor(kRed+2);
  m_hImpAtThisStep[muValue]->SetFillColor(kRed-10);
  m_hImpAtThisStep[muValue]->SetLineWidth(2);
  m_hImpAtThisStep[muValue]
    ->GetXaxis()->SetTitle("Retry giving lower -log(#it{L})");
  m_hImpAtThisStep[muValue]->GetYaxis()->SetTitle("Mean -#Delta log(#it{L})");
  gPad->SetLogy();
  m_hImpAtThisStep[muValue]->Draw("hist");
  can->Print(Form("%s/plot_AvgImprovementPerRetry.eps", m_outputDir.Data()));
  can->Clear();
  gPad->SetLogy(0);
  

  m_h2RetriesZ[muValue];
  m_h2RetriesZ[muValue]
    ->GetXaxis()->SetTitle("Retry with minimum -log(#it{L})");
  m_h2RetriesZ[muValue]->GetYaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  m_h2RetriesZ[muValue]->Draw("colz");
  can->Print(Form("%s/plot_2D_RetriesVsZ0.eps", m_outputDir.Data()));
  can->Clear();
  

  m_h2ZImprovement[muValue];
  m_h2ZImprovement[muValue];
  m_h2ZImprovement[muValue]->GetYaxis()->SetTitle("Improvement on minimum");
  m_h2ZImprovement[muValue]->GetXaxis()->SetTitle("Z_{0}^{Local} [#sigma]");
  m_h2ZImprovement[muValue]->Draw("colz");
  can->Print(Form("%s/plot_2D_Z0VsImprovement.eps", m_outputDir.Data()));
  can->Clear();
 
  
  ///////////////////////
  TGraph *gGaussMean = new TGraph();
  TGraph *gGaussSigma = new TGraph();
  TGraph *gGaussZ0 = new TGraph();
  // m_hZ0Retries[muValue][i_r]->Fill(currValueZ0);
  for (int i_r = 0; i_r < 50; i_r++) {
    TF1 *fGauss = new TF1("fGauss", "gaus", 
			  m_hZ0Retries[muValue][i_r]->GetXaxis()->GetXmin(), 
			  m_hZ0Retries[muValue][i_r]->GetXaxis()->GetXmax());
    m_hZ0Retries[muValue][i_r]->Fit(fGauss, "0");
    gGaussMean->SetPoint(i_r, i_r, fGauss->GetParameter(1));
    gGaussSigma->SetPoint(i_r, i_r, fGauss->GetParameter(2));
    
    
    double totalIntegral 
      = fGauss->Integral(m_hZ0Retries[muValue][i_r]->GetXaxis()->GetXmin(), 
			 m_hZ0Retries[muValue][i_r]->GetXaxis()->GetXmax());
    double gaussIntegral 
      = fGauss->Integral(m_config->getNum("GlobalP0AnalysisSigma"),
			 m_hZ0Retries[muValue][i_r]->GetXaxis()->GetXmax());
    gGaussZ0->SetPoint(i_r, i_r,
		       TMath::NormQuantile(1.0-(gaussIntegral/totalIntegral)));
    
    delete fGauss;
  }
  
  gGaussMean->GetXaxis()->SetTitle("Retry");
  gGaussMean->GetYaxis()->SetTitle("Z_{0} Mean");
  gGaussMean->Draw("AEP");
  can->Print(Form("%s/graph_Z0MeanVsRetry.eps", m_outputDir.Data()));
  can->Clear();

  gGaussSigma->GetXaxis()->SetTitle("Retry");
  gGaussSigma->GetYaxis()->SetTitle("Z_{0} Standard deviation");
  gGaussSigma->Draw("AEP");
  can->Print(Form("%s/graph_Z0SigmaVsRetry.eps", m_outputDir.Data()));
  can->Clear();
  
  gGaussZ0->GetXaxis()->SetTitle("Retry");
  gGaussZ0->GetYaxis()->SetTitle("Z_{0}^{Global}");
  gGaussZ0->Draw("AEP");
  can->Print(Form("%s/graph_GlobalZ0VsRetry.eps", m_outputDir.Data()));
  can->Clear();
  
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
  //hStatMu0->GetYaxis()->SetRangeUser(0.0001,1.0);
  hStatMu0->Draw("");
  hStatMu1->Draw("SAME");
  
  TH1F *hAsymptotic = getAsymptoticHist(statistic);
  hAsymptotic->SetLineColor(kBlack);
  hAsymptotic->Draw("SAME");
  
  TLegend leg(0.49, 0.76, 0.84, 0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  leg.AddEntry(hStatMu1, "#mu=1 toy MC","l");
  leg.AddEntry(hStatMu0, "#mu=0 toy MC","l");
  leg.AddEntry(hAsymptotic, "Asyptotic distribution","l");
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
  
  // Plot settings:
  double yTitleOffset = 0.4;
  double yTitleSize = 0.15;
  double yLabelOffset = 0.005;
  double yLabelSize = 0.12;
  
  // Colors:
  Color_t fillColorAsym = kBlue-9;
  Color_t lineColorAsym = kBlue+1;
  Color_t fillColorToy = kRed-9;
  Color_t lineColorToy = kRed+1;
  
  // Get the toy and asymptotic distributions:
  TH1F *hStatNominal = NULL;
  if (statistic.EqualTo("QMu")) hStatNominal = getStatHist(statistic, 1);
  else if (statistic.EqualTo("Q0")) hStatNominal = getStatHist(statistic, 0);
  else printer("ToyAnalysis: ERROR don't know which hist to get", true);
  TH1F *hAsymptotic = getAsymptoticHist(statistic);
  
  // Create histograms for sub-plot axes:
  TH1F *hCDF_Toy = new TH1F("hCDF_Toy","hCDF_Toy", m_nBins, m_binMin, m_binMax);
  TH1F *hRatio = new TH1F("hRatio", "hRatio", m_nBins, m_binMin, m_binMax);
  TH1F *hZ0_Toy = new TH1F("hZ0_Toy", "hZ0_Toy", m_nBins, m_binMin, m_binMax);
  
  // Test statistics arrays:
  double value_Q[m_nBins] = {0.0};
  double error_Q[m_nBins] = {(0.5*(m_binMax-m_binMin)/((double)m_nBins))};
  
  // CDF arrays:
  double value_CDF_Asym[m_nBins] = {0.0};
  double value_CDF_Toy[m_nBins] = {0.0};
  double errorLo_CDF_Toy[m_nBins] = {0.0};
  double errorHi_CDF_Toy[m_nBins] = {0.0};
  
  // Z0 arrays:
  double value_Z0_Asym[m_nBins] = {0.0};
  double value_Z0_Toy[m_nBins] = {0.0};
  double errorLo_Z0_Toy[m_nBins] = {0.0};
  double errorHi_Z0_Toy[m_nBins] = {0.0};
  
  // Ratio arrays:
  double value_Ratio[m_nBins] = {0.0};
  double errorLo_Ratio[m_nBins] = {0.0};
  double errorHi_Ratio[m_nBins] = {0.0};
  
  // Asymptotic Graphs:
  TGraph *gCDF_Asym = new TGraphAsymmErrors();
  TGraph *gZ0_Asym = new TGraphAsymmErrors();
  
  int nToys = hStatNominal->GetEntries();
  double yMinimum = 0.51 / ((double)nToys);
  double yMaximum = 0.99;

  // Loop over histogram bins to set values:
  int toyPoints = 0; 
  for (int i_b = 1; i_b <= m_nBins; i_b++) {
    
    // p-value calculations for CDF graphs:
    double p_Asym = hAsymptotic->Integral(i_b, m_nBins);
    double p_Toy = hStatNominal->Integral(i_b, m_nBins);
    double pError_Toy = calculateErrorPVal(p_Toy, nToys);
    
    // Z0 graphs for significance graphs:
    double z0_Asym = TMath::NormQuantile(1-p_Asym);
    double z0_Toy = TMath::NormQuantile(1-p_Toy);
    double z0Ratio = (z0_Toy == 0) ? 1.0 : (z0_Asym / z0_Toy);
    
    // Histograms for plot axes:
    hCDF_Toy->SetBinContent(i_b, p_Toy);
    hZ0_Toy->SetBinContent(i_b, z0_Toy);
    hRatio->SetBinContent(i_b, z0Ratio);
        
    // Then fill graph values:
    gCDF_Asym->SetPoint(i_b-1, hStatNominal->GetBinCenter(i_b), p_Asym);
    gZ0_Asym->SetPoint(i_b-1, hStatNominal->GetBinCenter(i_b), z0_Asym);
    
    // Fill arrays for graphs:
    if (p_Toy > 0 && i_b > 0) {
      
      value_Q[toyPoints] = hStatNominal->GetBinCenter(i_b);
      
      // P0 (CDF) values:
      value_CDF_Asym[toyPoints] = p_Asym;
      value_CDF_Toy[toyPoints] = p_Toy;
      errorLo_CDF_Toy[toyPoints] = (pError_Toy < p_Toy) ? pError_Toy : p_Toy;
      errorHi_CDF_Toy[toyPoints] = pError_Toy;
      
      // Z and p errors are flipped:
      double z0ErrorLo_Toy
	= TMath::NormQuantile(1.0 - (p_Toy+errorHi_CDF_Toy[toyPoints]));
      double z0ErrorHi_Toy
	= TMath::NormQuantile(1.0 - (p_Toy-errorLo_CDF_Toy[toyPoints]));
      
      if (!std::isfinite(z0ErrorLo_Toy) || std::isnan(z0ErrorLo_Toy)) {
	z0ErrorLo_Toy = 0.0;
      }
      if (!std::isfinite(z0ErrorHi_Toy) || std::isnan(z0ErrorHi_Toy)) {
	z0ErrorHi_Toy = 0.0;
      }

      double z0RatioLo = z0_Asym / z0ErrorLo_Toy;
      double z0RatioHi = z0_Asym / z0ErrorHi_Toy;
      
      // Z0 graph values:
      value_Z0_Asym[toyPoints] = z0_Asym;
      value_Z0_Toy[toyPoints] = z0_Toy;
      errorLo_Z0_Toy[toyPoints] = fabs(z0_Toy - z0ErrorLo_Toy);
      errorHi_Z0_Toy[toyPoints] = fabs(z0_Toy - z0ErrorHi_Toy);
      
      // Ratio graph values:
      value_Ratio[toyPoints] = z0Ratio;
      errorLo_Ratio[toyPoints] = fabs(z0Ratio - z0RatioHi);
      errorHi_Ratio[toyPoints] = fabs(z0Ratio - z0RatioLo);
      
      if (!std::isfinite(errorLo_Ratio[toyPoints]) || 
	  std::isnan(errorLo_Ratio[toyPoints])) {
	errorLo_Ratio[toyPoints] = 0.0;
      }
      if (!std::isfinite(errorHi_Ratio[toyPoints]) || 
	  std::isnan(errorHi_Ratio[toyPoints])) {
	errorHi_Ratio[toyPoints] = 0.0;
      }
      
      toyPoints++;
    }
  }
    
  //---------- Pad 1: Draw test statistic ----------//
  pad1->cd();
  hStatNominal->SetLineWidth(2);
  hStatNominal->SetLineColor(lineColorToy);
  hStatNominal->SetFillColor(lineColorToy);
  hStatNominal->SetFillStyle(3354);
  hStatNominal->GetYaxis()->SetRangeUser(yMinimum, yMaximum);
  hStatNominal->GetYaxis()->SetTitleSize(yTitleSize/2.0);
  hStatNominal->GetYaxis()->SetTitleOffset(yTitleOffset*2.0);
  hStatNominal->GetYaxis()->SetLabelSize(yLabelSize/2.0);
  hStatNominal->GetYaxis()->SetLabelOffset(yLabelOffset*2.0);
  double xWidth = ((m_binMax- m_binMin) / ((double)m_nBins));
  hStatNominal->GetYaxis()->SetTitle(Form("Fraction / %2.2f",xWidth));
  hStatNominal->Draw("");
  hAsymptotic->SetLineWidth(2);
  hAsymptotic->SetLineColor(lineColorAsym);
  hAsymptotic->SetFillColor(lineColorAsym);
  hAsymptotic->SetFillStyle(3345);
  
  hAsymptotic->Draw("SAME");
  gPad->SetLogy();
  //gPad->SetLogx();
 
  // Draw the ATLAS text:
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.07);
  t.DrawLatex(0.63, 0.86, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.07);
  t.DrawLatex(0.75, 0.86, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.63, 0.78, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (m_config->getNum("AnalysisLuminosity")/
				1000.0)));
  
  // Print the legend:
  TLegend leg(0.63, 0.58, 0.92, 0.75);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.07);
  leg.SetTextFont(42);
  TString printName = printStatName(statistic);
  if (statistic.EqualTo("QMu")) {
    leg.AddEntry(hStatNominal, Form("#mu=1 Toy MC %s",printName.Data()), "F");
  }
  else if (statistic.EqualTo("Q0")) {
    leg.AddEntry(hStatNominal, Form("#mu=0 Toy MC %s",printName.Data()), "F");
  }
  
  leg.AddEntry(hAsymptotic, Form("Asymptotic %s",printName.Data()), "L");
  leg.Draw("SAME");
  
  //---------- Pad 2: Draw the CDFs ----------//
  pad2->cd();
  
  // Axis formatting:
  hCDF_Toy->SetLineColor(0);    
  hCDF_Toy->GetXaxis()->SetTitle(printName);
  hCDF_Toy->GetYaxis()->SetTitle("CDF");
  hCDF_Toy->GetYaxis()->SetTitleSize(yTitleSize);
  hCDF_Toy->GetYaxis()->SetTitleOffset(yTitleOffset);
  hCDF_Toy->GetYaxis()->SetLabelSize(yLabelSize);
  hCDF_Toy->GetYaxis()->SetLabelOffset(yLabelOffset);
  hCDF_Toy->GetYaxis()->SetRangeUser(yMinimum, yMaximum);
  hCDF_Toy->Draw("");
  
  TGraphAsymmErrors *gCDF_Toy
    = new TGraphAsymmErrors(toyPoints, value_Q, value_CDF_Toy, error_Q,
			    error_Q, errorLo_CDF_Toy, errorHi_CDF_Toy);
  gCDF_Toy->SetLineWidth(2);
  gCDF_Toy->SetLineColor(lineColorToy);
  gCDF_Toy->SetFillColor(fillColorToy);
  gCDF_Toy->Draw("3SAME");
  gCDF_Asym->SetLineColor(lineColorAsym);
  gCDF_Asym->SetLineWidth(2);
  gCDF_Asym->Draw("LSAME");
  
  gPad->SetLogy();
  //gPad->SetLogx();
  
  //---------- Pad 3: Draw the Z0 ----------//
  pad3->cd();
  
  double yMin3 = 0.01; double yMax3 = 4.49;
  hZ0_Toy->SetLineColor(0);
  hZ0_Toy->GetXaxis()->SetTitle(printName);
  hZ0_Toy->GetYaxis()->SetTitle("Z [#sigma]");
  hZ0_Toy->GetYaxis()->SetTitleSize(yTitleSize);
  hZ0_Toy->GetYaxis()->SetTitleOffset(yTitleOffset);
  hZ0_Toy->GetYaxis()->SetLabelSize(yLabelSize);
  hZ0_Toy->GetYaxis()->SetLabelOffset(yLabelOffset);
  hZ0_Toy->GetYaxis()->SetRangeUser(yMin3, yMax3);
  hZ0_Toy->GetYaxis()->SetNdivisions(6);
  hZ0_Toy->Draw("");
  
  TGraphAsymmErrors *gZ0_Toy
    = new TGraphAsymmErrors(toyPoints, value_Q, value_Z0_Toy, error_Q, error_Q, 
			    errorLo_Z0_Toy, errorHi_Z0_Toy);
  gZ0_Toy->SetLineWidth(2);
  gZ0_Toy->SetLineColor(lineColorToy);
  gZ0_Toy->SetFillColor(fillColorToy);
  gZ0_Toy->Draw("3SAME");
  gZ0_Asym->SetLineWidth(2);
  gZ0_Asym->SetLineColor(lineColorAsym);
  gZ0_Asym->GetYaxis()->SetTitle("Z [#sigma]");
  gZ0_Asym->Draw("LSAME");
  
  //gPad->SetLogx();
  
  //---------- Pad 4: Plot ratio of Z0 ----------//
  pad4->cd();
  
  double yMin4 = 0.61; double yMax4 = 1.39;
  hRatio->SetLineWidth(2);
  hRatio->SetLineColor(lineColorAsym);
  hRatio->GetYaxis()->SetTitle("Z_{Asym.} / Z_{Toy}");
  hRatio->GetXaxis()->SetTitle(printName);
  hRatio->GetYaxis()->SetTitleSize(yTitleSize);
  hRatio->GetYaxis()->SetTitleOffset(yTitleOffset);
  hRatio->GetYaxis()->SetLabelSize(yLabelSize);
  hRatio->GetYaxis()->SetLabelOffset(yLabelOffset);
  hRatio->GetYaxis()->SetNdivisions(6);
  hRatio->GetXaxis()->SetTitleSize(yTitleSize);
  hRatio->GetXaxis()->SetLabelSize(yLabelSize);
  hRatio->GetXaxis()->SetTitleOffset(1.0);
  hRatio->GetYaxis()->SetRangeUser(yMin4, yMax4);
  hRatio->Draw("");
  
  TGraphAsymmErrors *gRatio
    = new TGraphAsymmErrors(toyPoints, value_Q, value_Ratio, error_Q, error_Q, 
			    errorLo_Ratio, errorHi_Ratio);
  gRatio->SetLineWidth(2);
  gRatio->SetLineColor(lineColorAsym);
  gRatio->SetFillColor(fillColorAsym);
  gRatio->Draw("3SAME");
  
  //gPad->SetLogx();
  
  // Draw a line at zero significance
  TLine *line2 = new TLine();
  line2->SetLineStyle(1);
  line2->SetLineWidth(1);
  line2->SetLineColor(kBlack);
  line2->DrawLine(m_binMin, 1.0, m_binMax, 1.0);
  line2->SetLineColor(kGray+2);
  line2->SetLineStyle(3);
  line2->DrawLine(m_binMin, 1.3, m_binMax, 1.3);
  line2->DrawLine(m_binMin, 1.2, m_binMax, 1.2);
  line2->DrawLine(m_binMin, 1.1, m_binMax, 1.1);
  line2->DrawLine(m_binMin, 0.9, m_binMax, 0.9);
  line2->DrawLine(m_binMin, 0.8, m_binMax, 0.8);
  line2->DrawLine(m_binMin, 0.7, m_binMax, 0.7);
  //hRatio->Draw("SAME");
  
  // Draw a line at each sigma value:
  TLatex t2; t2.SetTextColor(kRed+1); t.SetTextFont(42); t.SetTextSize(0.1);
  
  TLine *line3 = new TLine();
  line3->SetLineWidth(1);
  line3->SetLineStyle(2);
  line3->SetLineColor(kRed+1);
  double sigmaValues[3] = {1.0, 4.0, 9.0};
  for (int i_s = 0; i_s < 3; i_s++) {
    pad1->cd();
    line3->DrawLine(sigmaValues[i_s], yMinimum, sigmaValues[i_s], yMaximum);
    t2.DrawLatex(sigmaValues[i_s]+0.25, 0.01, Form("%d#sigma",i_s+1));
    pad2->cd();
    line3->DrawLine(sigmaValues[i_s], yMinimum, sigmaValues[i_s], yMaximum);
    pad3->cd();
    line3->DrawLine(sigmaValues[i_s], yMin3, sigmaValues[i_s], yMax3);
    pad4->cd();
    line3->DrawLine(sigmaValues[i_s], yMin4, sigmaValues[i_s], yMax4);
  }
  
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

/**
   -----------------------------------------------------------------------------
   Change the binning and x-axis range for the statistical histogram.
   @param nBins - The number of histogram bins.
   @param binMin - The minimum of the x-axis.
   @param binMax - The maximum of the x-axis.
*/
void ToyAnalysis::setStatHistRanges(int nBins, int binMin, int binMax) {
  m_nBins = nBins;
  m_binMin = binMin;
  m_binMax = binMax;
}

/**
   -----------------------------------------------------------------------------
   Simultaneously sort vector 1 and vector 2. Vector1's values will be used for
   sorting comparison. Resulting vector will be sorted in ascending order. The
   input vectors are passed by reference and their ordering is modified. 
   @param vec1 - The first vector to be sorted. It contains values for ordering.
   @param vec2 - The second vector to be sorted. Values don't affect ordering.
*/
void ToyAnalysis::sortPairedVectors(std::vector<double> &vec1, 
				    std::vector<double> &vec2) {
  
  std::vector<double> vecSorted1; vecSorted1.clear();
  std::vector<double> vecSorted2; vecSorted2.clear();
  
  // Loop over the original vector:
  for (int i_n = 0; i_n < (int)vec1.size(); i_n++) {
    // If sorted vector is empty, add the current unsorted elements:
    if (i_n == 0) {
      vecSorted1.push_back(vec1[i_n]);
      vecSorted2.push_back(vec2[i_n]);
    }
    // If sorted vector not empty, loop over sorted vector:
    else {
      std::vector<double>::iterator iterVec1 = vecSorted1.begin();
      std::vector<double>::iterator iterVec2 = vecSorted2.begin();
      
      bool inserted = false;
      while (iterVec2 != vecSorted2.end()) {
	// Insert current unsorted value if it is less than current sorted value
	if (vec1[i_n] < *iterVec1) {
	  vecSorted1.insert(iterVec1, vec1[i_n]);
	  vecSorted2.insert(iterVec2, vec2[i_n]);
	  inserted = true;
	  break;
	}
	else {
	  iterVec1++;
	  iterVec2++;
	}
      }
      // If current unsorted values weren't inserted yet, insert them now. 
      if (!inserted) {
	vecSorted1.push_back(vec1[i_n]);
	vecSorted2.push_back(vec2[i_n]);
      }
    }
  }
  
  // Then equate the original vector with the sorted vector:
  vec1 = vecSorted1;
  vec2 = vecSorted2;
}
