////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: GenericToys.cxx                                                     //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 22/03/2016                                                          //
//                                                                            //
//  This program tosses toy Monte Carlo for various studies.                  //
//                                                                            //
//  Options:                                                                  //
//    "ToyImportSamp"                                                         //
//    "ToyForScan"                                                            //
//    "ToyMLPoint"                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Package includes:
#include "GenericToys.h"

/**
   -----------------------------------------------------------------------------
   Convert key and value in a map to two separate vectors that are passed by 
   reference. Stupid, yeah...
   @param map - The input map with keys and values to be split.
   @names - A vector of names.
   @values - A vector of values.
*/
void mapToVectors(std::map<std::string,double> map, 
		  std::vector<std::string>& names, std::vector<double>& values){
  names.clear();
  values.clear();
  for (std::map<std::string,double>::iterator mapIter = map.begin(); 
       mapIter != map.end(); mapIter++) {
    TString currName = TString(mapIter->first);
    if (!currName.Contains("gamma_stat_channel_bin")) {
      names.push_back(mapIter->first);
      values.push_back(mapIter->second);
    }
  }
}

/**
   -----------------------------------------------------------------------------
   Create a single pseudo-experiment and fit it. If option "ToyMLPoint" is used,
   the program will attempt to load the ML snapshot to get the values of the
   parameters of interest. On the other hand, if "ToyForScan" is used in the
   option, the parameters of interest are taken directly from the arguments
   passed to the main() method.
   @param doImportanceSampling - True iff you want to use importance sampling
   to inflate the statistics in the tail of the LLR distributions.
*/
void runToysSinglePoint(bool doImportanceSampling) {
  std::cout << "GenericToys::runToysSinglePoint(" << (int)doImportanceSampling
	    << ")" << std::endl;
  
  // Default snapshot to be used in the toy generation and prior to fitting:
  TString snapshotName = (m_inputPoIVal == 0) ? 
    m_config->getStr("WorkspaceSnapshotMu0") :
    m_config->getStr("WorkspaceSnapshotMu1");
  
  // Load model, data, etc. from workspace:
  TString workspaceFileName = m_config->getStr("WorkspaceFile");
  TObjArray *array = workspaceFileName.Tokenize("/");
  workspaceFileName
    = ((TObjString*)array->At(array->GetEntries()-1))->GetString();
  TFile *inputFile = new TFile(workspaceFileName);
  if (inputFile->IsZombie()) {
    workspaceFileName = m_config->getStr("WorkspaceFile");
    inputFile = new TFile(workspaceFileName, "read");
  }
  RooWorkspace *workspace
    = (RooWorkspace*)inputFile->Get(m_config->getStr("WorkspaceName"));
    
  //----------------------------------------//
  // Prepare model parameters for tossing and fitting toy MC:

  // The statistics class, for calculating qMu etc. 
  TestStat *testStat = new TestStat(m_configFile, "new", workspace);
    
  // Set the PoI ranges for this study:
  std::vector<TString> listPoI = m_config->getStrV("WorkspacePoIs");
  for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
    std::vector<double> currRange
      = m_config->getNumV(Form("PoIRange_%s", (listPoI[i_p]).Data()));
    if (testStat->theWorkspace()->var(listPoI[i_p])) {
      testStat->theWorkspace()->var(listPoI[i_p])
	->setRange(currRange[0], currRange[1]);
    }
    else {
      std::cout << "GlobalP0Toys: Workspace has no variable " << listPoI[i_p]
		<< std::endl;
      exit(0);
    }
  }
    
  // Turn off MC stat errors if requested for Graviton jobs:
  if (m_config->isDefined("TurnOffTemplateStat") && 
      m_config->getBool("TurnOffTemplateStat")) {
    testStat->theWorkspace()->loadSnapshot(snapshotName);
    const RooArgSet *nuisanceParameters
      = testStat->theModelConfig()->GetNuisanceParameters();
    TIterator *nuisIter = nuisanceParameters->createIterator();
    RooRealVar *nuisCurr = NULL;
    while ((nuisCurr = (RooRealVar*)nuisIter->Next())) {
      TString currName = nuisCurr->GetName();
      if (currName.Contains("gamma_stat_channel_bin")) {
	testStat->setParam(currName, nuisCurr->getVal(), true);
      }
    }
  }
  
  // Get the name of the normalization parameter:
  TString poiForNorm = m_config->getStr("PoIForNormalization");
  TString poiForMass = m_config->getStr("PoIForMass");
  TString poiForWidth = m_config->getStr("PoIForWidth");
  
  // For setting the PoI values in toy creation and fitting:
  std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
  std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
  
  // For a maximum-likelihood fit, set the parameters of interest using mu = 1 
  // snapshot. The mu=1 fit parameters of interest are used even for the mu=0 
  // parameters of interest because the background-only fits and toys should 
  // still have the spurious signal randomized and fitted at the proper mass.
  if (m_options.Contains("ToyMLPoint")) {
    std::cout << "GenericToys: Creating toy with ToyMLPoint option."
	      << std::endl;

    // Load the background-only snapshot (signal + background snapshot)
    // for background-only toys (signal + background toys).
    snapshotName = (m_inputPoIVal == 0) ? 
      m_config->getStr("WorkspaceSnapshotMu0") :
      m_config->getStr("WorkspaceSnapshotMu1");
    testStat->setNominalSnapshot(snapshotName);
    
    const RooArgSet *snapshotMu1
      = workspace->getSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
    TIterator *snapMu1Iter = snapshotMu1->createIterator();
    RooRealVar *snapParMu1 = NULL;
    while ((snapParMu1 = (RooRealVar*)snapMu1Iter->Next())) {
      TString currName = snapParMu1->GetName();
      for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
	if (currName.EqualTo(listPoI[i_p])) {
	  mapPoIMu1[listPoI[i_p]] = snapParMu1->getVal();
	  if (currName.EqualTo(poiForNorm)) mapPoIMu0[listPoI[i_p]] = 0.0;
	  else mapPoIMu0[listPoI[i_p]] = snapParMu1->getVal();
	}
      }
      // Set the parameters for storage in the TTree:
      if (currName.EqualTo(poiForNorm)) m_toyXSection = snapParMu1->getVal();
      if (currName.EqualTo(poiForMass)) m_toyMass = snapParMu1->getVal();
      if (currName.EqualTo(poiForWidth)) m_toyWidth = snapParMu1->getVal();
    }
  }
  // The snaphshots for fits at the scan points probably don't exist, so we must
  // perform these fits before doing anything else. 
  else if (m_options.Contains("ToyForScan")) {
    std::cout << "GenericToys: Creating toy with ToyForScan option and mass="
	      << m_toyMass << ", width=" << m_toyWidth << ", xs= " 
	      << m_toyXSection << std::endl;
    
    // Map for Mu = 0 parameters of interest:
    mapPoIMu0[m_config->getStr("PoIForNormalization")] = 0.0;
    mapPoIMu0[m_config->getStr("PoIForMass")] = m_toyMass;
    mapPoIMu0[m_config->getStr("PoIForWidth")] = m_toyWidth;
    
    // Map for Mu = 1 parameters of interest:
    mapPoIMu1[m_config->getStr("PoIForNormalization")] = m_toyXSection;
    mapPoIMu1[m_config->getStr("PoIForMass")] = m_toyMass;
    mapPoIMu1[m_config->getStr("PoIForWidth")] = m_toyWidth;
    
    // Be sure that snapshot is saved in workspace during fit:
    testStat->setNominalSnapshot(snapshotName);
    testStat->saveSnapshots(true);
    
    // Perform the fit:
    std::cout << "GenericToys: Pre-toy fit (for scan) starting" << std::endl;
    if (m_inputPoIVal == 0) {
      testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
			  0, true, mapPoIMu0, false);
    }
    else {
      testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
			  1, true, mapPoIMu1, false);
    }
    
    // Set the new snapshot as default:
    snapshotName = Form("paramsProfilePoI%d", m_inputPoIVal);
    testStat->setNominalSnapshot(snapshotName);
  }
  
  //----------------------------------------//
  // Preparations for importance sampling:
  int nImportanceDensities = 1;
  int importanceDensityToUse = 0;
  std::map<TString,double> mapPoIIS; mapPoIIS.clear();
  std::map<TString,double> mapPoIIS_Lo; mapPoIIS_Lo.clear();
  std::map<TString,double> mapPoIIS_Hi; mapPoIIS_Hi.clear();
  if (doImportanceSampling) {
    
    // Use ML snapshot (most extreme p0 case) to derive # importance densities:
    const RooArgSet *snapshotMu1
      = workspace->getSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
    TIterator *snapMu1Iter = snapshotMu1->createIterator();
    RooRealVar *snapParMu1 = NULL;
    while ((snapParMu1 = (RooRealVar*)snapMu1Iter->Next())) {
      TString currName = snapParMu1->GetName();
      if (currName.EqualTo(poiForNorm)) {
	// Use muhat to get number of importance densities (muhat / error):
	nImportanceDensities
	  = (int)std::round(snapParMu1->getVal() / snapParMu1->getError());
	break;
      }
    }
    
    // Generate a random importance density to use in this toy:
    TRandom3 randomIDGen = TRandom3(19*m_seed+7);
    importanceDensityToUse = randomIDGen.Integer(nImportanceDensities+1);
    
    // Then set the mapPoI normalization parameters for pseudo-data:
    if (m_inputPoIVal == 0) {
      mapPoIIS = mapPoIMu0;
      mapPoIIS_Lo = mapPoIMu0;
      mapPoIIS_Hi = mapPoIMu0;
      mapPoIIS[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)importanceDensityToUse / (double)nImportanceDensities);
      mapPoIIS_Lo[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)(importanceDensityToUse-1) / (double)nImportanceDensities);
      mapPoIIS_Hi[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)(importanceDensityToUse+1) / (double)nImportanceDensities);
    }
    else {
      mapPoIIS = mapPoIMu1;
      mapPoIIS_Lo = mapPoIMu1;
      mapPoIIS_Hi = mapPoIMu1;
      mapPoIIS[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)(nImportanceDensities-importanceDensityToUse) / 
	 (double)nImportanceDensities); 
      mapPoIIS_Lo[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)(nImportanceDensities-(importanceDensityToUse-1)) / 
	 (double)nImportanceDensities); 
      mapPoIIS_Hi[m_config->getStr("PoIForNormalization")] = m_toyXSection * 
	((double)(nImportanceDensities-(importanceDensityToUse+1)) / 
	 (double)nImportanceDensities); 
    }
      
    // Importance sampling is unnecessary in absence of extreme excess:
    if (nImportanceDensities < 1) doImportanceSampling = false;
  }
  
  //----------------------------------------//
  // Create the pseudo data:
  if (doImportanceSampling) {
    testStat->createPseudoData(m_seed, m_inputPoIVal, snapshotName, mapPoIIS);
  }
  else if (m_inputPoIVal == 0) {
    testStat->createPseudoData(m_seed, m_inputPoIVal, snapshotName, mapPoIMu0);
  }
  else {
    testStat->createPseudoData(m_seed, m_inputPoIVal, snapshotName, mapPoIMu1);
  }
  
  m_numEvents = workspace->data("toyData")->sumEntries();
  m_numEventsPerCate = testStat->getNEventsToys();
  
  // Globs are only randomized once for each dataset:
  mapToVectors(testStat->getGlobalObservables(), m_namesGlobs,
	       m_valuesGlobsMu0);
  mapToVectors(testStat->getGlobalObservables(), m_namesGlobs,
	       m_valuesGlobsMu1);
  mapToVectors(testStat->getGlobalObservables(), m_namesGlobs, 
	       m_valuesGlobsMuFree);
  
  // Strategy settings (if defined):
  if (m_config->isDefined("FitOptions")) {
    testStat->setFitOptions(m_config->getStr("FitOptions"));
  }
  
  //----------------------------------------//
  // Perform the fits to toys:
  
  // Mu = 0 fits:
  std::cout << "GenericToys: Mu=0 fit starting" << std::endl;
  m_nllMu0 = testStat->getFitNLL("toyData", 0, true, mapPoIMu0, false);
  m_convergedMu0 = testStat->fitsAllConverged();
  mapToVectors(testStat->getNuisanceParameters(), m_namesNP, m_valuesNPMu0);
  mapToVectors(testStat->getPoIs(), m_namesPoIs, m_valuesPoIsMu0);
  
  // Mu = 1 fits:
  std::cout << "GenericToys: Mu=1 fit starting" << std::endl;
  m_nllMu1 = testStat->getFitNLL("toyData", 1, true, mapPoIMu1, false);
  m_convergedMu1 = testStat->fitsAllConverged();
  mapToVectors(testStat->getNuisanceParameters(), m_namesNP, m_valuesNPMu1);
  mapToVectors(testStat->getPoIs(), m_namesPoIs, m_valuesPoIsMu1);
  double nominalNormalization
    = (testStat->getPoIs())[(std::string)(poiForNorm)];
  
  // Fix the mass and width before the mu-free fit:
  // Loop over the parameters of interest, and set all of those that are
  // not the cross-section constant and to the proper values
  for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
    // Set every parameter constant to mu=1 value EXCEPT cross-section:
    if (!(listPoI[i_p]).EqualTo(m_config->getStr("PoIForNormalization"))) {
      // get the value of the PoI to set from the previous fit:
      for (int i_n = 0; i_n < (int)m_namesPoIs.size(); i_n++) {
	if (TString(m_namesPoIs[i_n]).EqualTo(listPoI[i_p])) {
	  testStat->setParam(listPoI[i_p], m_valuesPoIsMu1[i_n], true);
	  break;
	}
      }
    }
  }
  
  // Mu free fits:
  std::cout << "GenericToys: Mu-free fit starting" << std::endl;
  m_nllMuFree = testStat->getFitNLL("toyData", 1, false, mapPoIMu1, false);
  m_convergedMuFree = testStat->fitsAllConverged();
  mapToVectors(testStat->getNuisanceParameters(), m_namesNP, m_valuesNPMuFree);
  mapToVectors(testStat->getPoIs(), m_namesPoIs, m_valuesPoIsMuFree);
  double profiledNormalization
    = (testStat->getPoIs())[(std::string)(poiForNorm)];

  //----------------------------------------//
  // Calculate weight for importance sampling:
  if (doImportanceSampling) {
    double nllIS = 1.0; double nllIS_Lo = 1.0; double nllIS_Hi = 1.0;
    double nllNom = (m_inputPoIVal == 0) ? m_nllMu0 : m_nllMu1;
    
    // Do fit with IS
    std::cout << "GenericToys: IS fit starting" << std::endl;
    nllIS = testStat->getFitNLL("toyData", 1, true, mapPoIIS, false);
    // Do fit with IS_Lo if necessary:
    if (importanceDensityToUse > 0) {
      std::cout << "GenericToys: IS-low fit starting" << std::endl;
      nllIS_Lo = testStat->getFitNLL("toyData", 1, true, mapPoIIS_Lo, false);
    }
    // Do fit with IS_Hi if necessary:
    if (importanceDensityToUse < nImportanceDensities) {
      std::cout << "GenericToys: IS_hi fit starting" << std::endl;
      nllIS_Hi = testStat->getFitNLL("toyData", 1, true, mapPoIIS_Hi, false);
    }
    
    // Upper constraint only:
    if (importanceDensityToUse == 0) {
      if (nllIS <= nllIS_Hi) m_weight = 1.0;
      else m_weight = 0.0;
    }
    // Lower constraint only:
    else if (importanceDensityToUse == nImportanceDensities) {
      if (nllIS < nllIS_Lo) m_weight = TMath::Exp(nllIS - nllNom);
      else m_weight = 0.0;
    }
    // One of the sandwiched importance densities (upper and lower constraint):
    else {
      if (nllIS < nllIS_Lo && nllIS <= nllIS_Hi) {
	m_weight = TMath::Exp(nllIS - nllNom);
      }
      else m_weight = 0.0;
    }  
    
    std::cout << "GenericToys: Importance sample " << importanceDensityToUse
	      << " with weight of " << m_weight << std::endl;
    std::cout << " \tnllNom = " << nllNom << std::endl;
    std::cout << " \tnllIS = " << nllIS << std::endl;
    
  }
  
  //----------------------------------------//
  // Post-fit calculations:
  
  // Get the profiled PoI value:
  m_profiledPOIVal = (profiledNormalization / nominalNormalization);
  
  // Calculate profile likelihood ratios:
  m_llrL1L0 = m_nllMu1 - m_nllMu0;
  m_llrL1Lfree = m_profiledPOIVal > 1.0 ? 0.0 : (m_nllMu1 - m_nllMuFree);
  m_llrL0Lfree = m_profiledPOIVal < 0.0 ? 0.0 : (m_nllMu0 - m_nllMuFree);
  
  // Fill the tree:
  m_outputTree->Fill();
  
  // Count the toys:
  m_seed++;
  m_outputTree->AutoSave("SaveSelf");
  
  // Close the input file before the loop repeats:
  inputFile->Close();
}

/**
   -----------------------------------------------------------------------------
   The main method tosses toys and saves data in a TTree.
   @param configFile - The name of the analysis config file.
   @param options - The options (see header note).
   @param seed - The random seed for pseudoexperiment creation.
   @param toysPerJob - The number of pseudoexperiments to create per job.
   @param poiVal - The value of the parameter of interest to use.
*/
int main(int argc, char **argv) {
  if (argc < 6) {
    std::cout << "Usage: " << argv[0] 
	      << " <configFile> <options> <seed> <toysPerJob> <poiVal> " 
	      << std::endl;
    exit(0);
  }
  
  // Clock the toys:
  clock_t time;
  time = clock();
  
  // Assign input parameters:
  m_configFile = argv[1];
  m_options = argv[2];
  m_seed = atoi(argv[3]);
  m_nToysPerJob = atoi(argv[4]);
  m_inputPoIVal = atoi(argv[5]);
  
  // Use the options to see whether mass, width, and cross-section are set
  if (m_options.Contains("ToyForScan")) {
    if (argc < 9) {
      std::cout << "GenericToys: ERROR! ToyForScan needs 8 args." << std::endl;
      exit(0);
    }
    else {
      // Get the mass, width, xsection 
      m_toyMassInt = atoi(argv[6]);
      m_toyWidthInt = atoi(argv[7]);
      m_toyXSectionInt = atoi(argv[8]);
    }
  }
  else {
    m_toyMassInt = 0;
    m_toyWidthInt = 0;
    m_toyXSectionInt = 0;
  }
  
  // Load the analysis configurations from file:
  m_config = new Config(m_configFile);
    
  // Construct the output directories:
  m_outputDir = Form("%s/%s/GenericToys", 
		     (m_config->getStr("MasterOutput")).Data(),
		     (m_config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s/err", m_outputDir.Data()));
  system(Form("mkdir -vp %s/log", m_outputDir.Data()));
  system(Form("mkdir -vp %s/single_files", m_outputDir.Data()));
  
  // Name and create output TFile and TTree:
  TString tempFName = Form("%s/single_files/toy_mu%i_%i",
			   m_outputDir.Data(), m_inputPoIVal, m_seed);
  if (m_options.Contains("ToyForScan")) {
    tempFName.Append(Form("_ForScan_mass%d_width%d_xs%d.root",
			  m_toyMassInt, m_toyWidthInt, m_toyXSectionInt));
  }
  else if (m_options.Contains("ToyMLPoint")) tempFName.Append("_MLPoint.root");
  else tempFName.Append(".root");
  m_outputFile = new TFile(tempFName, "recreate");
  m_outputTree = new TTree("toy", "toy");
  
  // Initialize variables to store in the TTree:
  m_toyMass = (double)m_toyMassInt;
  m_toyWidth = ((double)m_toyWidthInt) / 100.0;
  m_toyXSection = ((double)m_toyXSectionInt) / 1000.0;
  
  m_bestFitUpdate = 0;
  m_weight = 1.0;
  m_numEvents = 0.0;
  m_convergedMu1 = true;
  m_convergedMu0 = true;
  m_convergedMuFree = true;
  m_profiledPOIVal = 0.0;

  m_nllMu0 = 0.0;
  m_nllMu1 = 0.0;
  m_nllMuFree = 0.0;
  m_llrL1L0 = 0.0;
  m_llrL0Lfree = 0.0;
  m_llrL1Lfree = 0.0;

  m_namesNP.clear();
  m_valuesNPMu0.clear();
  m_valuesNPMu1.clear();
  m_valuesNPMuFree.clear();
  m_namesGlobs.clear();
  m_valuesGlobsMu0.clear();
  m_valuesGlobsMu1.clear();
  m_valuesGlobsMuFree.clear();
  m_namesPoIs.clear();
  m_valuesPoIsMu0.clear();
  m_valuesPoIsMu1.clear();
  m_valuesPoIsMuFree.clear();

  m_numEventsPerCate.clear();
  m_nllPerRetry.clear();
  
  // Set the output TTree branch addresses:
  m_outputTree->Branch("seed", &m_seed, "seed/I");
  m_outputTree->Branch("weight", &m_weight, "weight/D");
  m_outputTree->Branch("toyMass", &m_toyMass, "toyMass/D");
  m_outputTree->Branch("toyWidth", &m_toyWidth, "toyWidth/D");
  m_outputTree->Branch("toyXSection", &m_toyXSection, "toyXSection/D");
  m_outputTree->Branch("bestFitUpdate", &m_bestFitUpdate, "bestFitUpdate/I");
  m_outputTree->Branch("numEvents", &m_numEvents, "numEvents/D");
  m_outputTree->Branch("numEventsPerCate", &m_numEventsPerCate);
  m_outputTree->Branch("nllPerRetry", &m_nllPerRetry);
  m_outputTree->Branch("profiledPOIVal", &m_profiledPOIVal, "profiledPOIVal/D");
  m_outputTree->Branch("convergedMu0", &m_convergedMu0, "convergedMu0/O");
  m_outputTree->Branch("convergedMu1", &m_convergedMu1, "convergedMu1/O");
  m_outputTree->Branch("convergedMuFree", &m_convergedMuFree, 
		       "convergedMuFree/O");
  m_outputTree->Branch("nllMu0", &m_nllMu0, "nllMu0/D");
  m_outputTree->Branch("nllMu1", &m_nllMu1, "nllMu1/D");
  m_outputTree->Branch("nllMuFree", &m_nllMuFree, "nllMuFree/D");
  m_outputTree->Branch("llrL1L0", &m_llrL1L0, "llrL1L0/D");
  m_outputTree->Branch("llrL0Lfree", &m_llrL0Lfree, "llrL0Lfree/D");
  m_outputTree->Branch("llrL1Lfree", &m_llrL1Lfree, "llrL1Lfree/D");
  m_outputTree->Branch("namesNP", &m_namesNP);
  m_outputTree->Branch("valuesNPMu0", &m_valuesNPMu0);
  m_outputTree->Branch("valuesNPMu1", &m_valuesNPMu1);
  m_outputTree->Branch("valuesNPMuFree", &m_valuesNPMuFree);
  m_outputTree->Branch("namesGlobs", &m_namesGlobs);
  m_outputTree->Branch("valuesGlobsMu1", &m_valuesGlobsMu1);
  m_outputTree->Branch("valuesGlobsMu0", &m_valuesGlobsMu0);
  m_outputTree->Branch("valuesGlobsMuFree", &m_valuesGlobsMuFree);
  m_outputTree->Branch("namesPoIs", &m_namesPoIs);
  m_outputTree->Branch("valuesPoIsMu1", &m_valuesPoIsMu1);
  m_outputTree->Branch("valuesPoIsMu0", &m_valuesPoIsMu0);
  m_outputTree->Branch("valuesPoIsMuFree", &m_valuesPoIsMuFree);
    
  //----------------------------------------//
  // Loop over seeds to generate pseudo experiments:
  std::cout << "GenericToys: Generating " << m_nToysPerJob << " toys." << endl;
  for (int i_t = 0; i_t < m_nToysPerJob; i_t++) {
    std::cout << "GenericToys: Starting toy " << i_t << " of " << m_nToysPerJob
	      << std::endl;
    
    if (m_options.Contains("ToyForScan") || m_options.Contains("ToyMLPoint")) {
      runToysSinglePoint(m_options.Contains("ToyImportSamp"));
    }
    else {
      std::cout << "GenericToys: Option " << m_options << " not recognized."
		<< std::endl;
    }
  }
  
  // Write the output file, delete local file copies:
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();
  
  // Clock the toys:
  time = clock() - time;
  
  if (m_config->getBool("Verbose")) {
    std::cout << "\nGenericToys: Toy procedure concluded." << std::endl;
    printf("\t%d toys required %d clock cycles (%f seconds).\n\n",
	   m_nToysPerJob, (int)time, ((float)time/CLOCKS_PER_SEC));
  }
  return 0;
}
