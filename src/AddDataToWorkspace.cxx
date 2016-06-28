////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  AddDataToWorkspace.cxx                                                    //
//                                                                            //
//  Author: Andrew Hard                                                       //
//  Date: 24/05/2016                                                          //
//  Email: ahard@cern.ch                                                      //
//                                                                            //
//  This macro processes data MxAODs and creates a new dataset for the        //
//  workspaces.                                                               //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "HGammaMxAOD.h"
#include "TestStat.h"

// Globally scoped variables:
HGammaMxAOD *m_treeMxAOD;
TString m_categoryName;

/**
   -----------------------------------------------------------------------------
   Select the proper category for the current event in the HGammaMxAOD. The 
   categorization to use is based on the global variable m_categoryName.
   @return - The category index.
*/
int chooseCategory() {

  // Eta categorization:
  if (m_categoryName.EqualTo("EtaCate")) {
    double barrelEnd = 1.37;
    double eta1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[0];
    double eta2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_eta)[1];
    if (fabs(eta1) < barrelEnd && fabs(eta2) < barrelEnd) return 0;
    else if (fabs(eta1) < barrelEnd) return 1;
    else if (fabs(eta2) < barrelEnd) return 2;
    else return 3;
  }

  // Conversion categorization:
  else if (m_categoryName.EqualTo("ConversionCate")) {
    int conv1 = (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[0];
    int conv2 = (*m_treeMxAOD->HGamPhotonsAuxDyn_conversionType)[1];
    if (conv1 == 0 && conv2 == 0) return 0;
    else if (conv1 == 0) return 1;
    else if (conv2 == 0) return 2;
    else return 3;
  }
 
  // Exit because inputs are unknown:
  else {
    std::cout << "AddDataToWorkspace: ERROR! Categorization not found" 
	      << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Get the name for the current category. The categorization to use is based on 
   the global variable m_categoryName.
   @param categoryIndex - The index of the current category.
   @return - The category name.
*/
TString nameCategory(int categoryIndex) {
  if (m_categoryName.EqualTo("EtaCate")) {
    if (categoryIndex == 0) return "BarrelBarrel";
    else if (categoryIndex == 1) return "BarrelEndcap";
    else if (categoryIndex == 2) return "EndcapBarrel";
    else return "EndcapEndcap";
  }
  else if (m_categoryName.EqualTo("ConversionCate")) {
    if (categoryIndex == 0) return "UnconvUnconv";
    else if (categoryIndex == 1) return "UnconvConv";
    else if (categoryIndex == 2) return "ConvUnconv";
    else return "ConvConv";
  }
  
  // Exit because inputs are unknown:
  else {
    std::cout << "AddDataToWorkspace: ERROR! Categorization not found" << std::endl;
    exit(0);
  }
}

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @return - An updated list of file names.
*/
std::vector<TString> makeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "AddDataToWorkspace: Making local copies of inputs."
	    << std::endl;
  std::vector<TString> result; result.clear();
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    TString newName = Form("tempFile%d.root", i_f);
    if (fileNames[i_f].Contains("root://eosatlas/")) {
      system(Form("xrdcp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else if (fileNames[i_f].Contains("/eos/atlas/")) {
      system(Form("eos cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    else {
      system(Form("cp %s %s", fileNames[i_f].Data(), newName.Data()));
    }
    result.push_back(newName);
  }
  return result;
}

/**
   -----------------------------------------------------------------------------
   Remove any files that were copied over for speed.
   @param fileNames - The original file names.
*/
void removeLocalFileCopies(std::vector<TString> fileNames) {
  std::cout << "AddDataToWorkspace: Removing local copies of inputs."
	    << std::endl;
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    system(Form("rm %s", fileNames[i_f].Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void printProgressBar(int index, int total) {
  if (index%100000 == 0) {
    TString print_bar = " [";
    for (int bar = 0; bar < 20; bar++) {
      double current_fraction = double(bar) / 20.0;
      if (double(index)/double(total) > current_fraction) print_bar.Append("/");
      else print_bar.Append(".");
    }
    print_bar.Append("] ");
    double percent = 100.0 * (double(index) / double(total));
    TString text = Form("%s %2.2f ", print_bar.Data(), percent);
    std::cout << text << "%\r" << std::flush; 
  }
}

/**
   -----------------------------------------------------------------------------
   The main method for this utility. Provide 1 argument - the location of the 
   config (.cfg) file, which should be stored in the data/ directory. The main()
   method runs over the samples provided, performs the fits requests, and gives
   comparisons of parameterized and non-parameterized fits. 
*/
int main(int argc, char *argv[])
{
  // Check that arguments are provided.
  if (argc < 2) {
    std::cout << "\nUsage: " << argv[0] << " <configFile> <options>"
	      << std::endl;
    exit(0);
  }
  
  Config *config = new Config(TString(argv[1]));
  TString options = argv[2];
  
  // Check that output directory exists:
  TString outputDir = Form("%s/%s/AddDataToWorkspace",
			   (config->getStr("MasterOutput")).Data(),
			   (config->getStr("JobName")).Data());
  system(Form("mkdir -vp %s", outputDir.Data()));
  
  // Set the plot Style to ATLAS defaults:
  CommonFunc::SetAtlasStyle();
  
  
  // See if a local copy of the workspace has been made:
  TString workspaceFileName = config->getStr("WorkspaceFile");
  TObjArray *array = workspaceFileName.Tokenize("/");
  workspaceFileName
    = ((TObjString*)array->At(array->GetEntries()-1))->GetString();
  TFile *workspaceFile = new TFile(workspaceFileName);
  if (workspaceFile->IsZombie()) {
    workspaceFileName = config->getStr("WorkspaceFile");
    workspaceFile = new TFile(workspaceFileName, "read");
  }

  // Load the workspace and all likelihood model objects:
  RooWorkspace *workspace = (RooWorkspace*)workspaceFile
    ->Get(config->getStr("WorkspaceName"));
  ModelConfig *model = (ModelConfig*)workspace
    ->obj(config->getStr("WorkspaceModelConfig"));
  RooSimultaneous *combPdf = (RooSimultaneous*)model->GetPdf();
  RooArgSet *observables = (RooArgSet*)model->GetObservables();
  TString nameRooCategory = config->getStr("WorkspaceRooCategory");
  RooCategory* categories = NULL;
  if (workspace->obj(nameRooCategory)) {
    categories = (RooCategory*)workspace->obj(nameRooCategory);
  }
  
  // Instantiate the RooDataSets:
  std::map<std::string,RooDataSet*> dataMap; dataMap.clear();
  std::string cate1 = "";
  RooRealVar *obs1 = NULL;
  
  // Iterate over categories:
  TIterator *cateIter = combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  while ((cateType = (RooCatType*)cateIter->Next())) {
    RooAbsPdf *currPdf = combPdf->getPdf(cateType->GetName());
    RooArgSet *currObs = currPdf->getObservables(observables);
    // Define the (inclusive) category name and observable variable:
    if (cate1 =="") cate1 = (std::string)(cateType->GetName());
    if (obs1 == NULL) obs1 = (RooRealVar*)currObs->first();
    // Define a new dataset:
    dataMap[(std::string)cateType->GetName()]
      = new RooDataSet(Form("data_%s", TString(cateType->GetName()).Data()),
		       Form("data_%s", TString(cateType->GetName()).Data()),
		       RooArgSet(*currObs));
  }
  
  // Define per-event histograms:
  int nCategories = config->getInt("MxAODNCategories");
  m_categoryName = config->getStr("MxAODCategorization");
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = config->getStrV("MxAODsForData");
  // Make local copies of files if requested, to improve speed:
  if (config->getBool("MakeLocalMxAODCopies")) {
    fileNames = makeLocalFileCopies(fileNames);
  }
  // Create TChain of input files:
  TChain *chain = new TChain(config->getStr("MxAODTreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  m_treeMxAOD = new HGammaMxAOD(chain, config->getStr("MxAODTag"));
  
  // Simple counter for rudimentary checks:
  int countPass = 0;
  int countCate[50] = {0};
    
  //--------------------------------------//
  // Loop over events to build dataset:
  int nEvents = m_treeMxAOD->fChain->GetEntries();
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {
    
    // Load event from MxAOD:
    m_treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);
    
    //---------- Select the event ----------//
    if ((config->getStr("AnalysisType")).Contains("Graviton") &&
	!(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedExotic && 
	  m_treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(m_treeMxAOD->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       m_treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
    
    // Choose the category:
    int category = chooseCategory();
    
    // Add to event counts:
    countPass++;
    countCate[category]++;
    
    // Set the mass observable:
    obs1->setVal(m_treeMxAOD->HGamEventInfoAuxDyn_m_yy/1000.0);
    
    // Fill the RooDataSet:
    dataMap[cate1]->add(RooArgSet(*obs1));
  
  }// End of loop over events

  // Add Ghost events:
  RooRealVar weight("weight", "weight", 1.0);
  if (config->isDefined("AddGhostEventsToData") && 
      config->isDefined("AddGhostEventsToData")) {
    for (double mass = obs1->getMin(); mass <= obs1->getMax(); mass += 1.0) {
      obs1->setVal(mass);
      weight.setVal(0.000001);
      dataMap[cate1]->add(RooArgSet(*obs1, weight), weight.getVal());
    }
  }
  
  // Remove files that were copied:
  if (config->getBool("MakeLocalMxAODCopies")) removeLocalFileCopies(fileNames);
  
  // Category and event yield summary:
  std::cout << "\nThere were " << countPass << " in the final sample."
	    << std::endl;
  for (int i_c = 0; i_c < nCategories; i_c++) {
    std::cout << "\t" << nameCategory(i_c)<< "\t" << countCate[i_c] 
	      << std::endl;
  }
  
  // Loop over the frames, creating combined datasets for each:
  TString dataName = "Data2016";
  
  RooDataSet *newData = NULL;
  if (config->isDefined("AddGhostEventsToData") && 
      config->isDefined("AddGhostEventsToData")) {
    observables->add(weight);
    newData = new RooDataSet(dataName, dataName, *observables, 
			     RooFit::Index(*categories),
			     RooFit::Import(dataMap), WeightVar(weight));
  }
  else {
    newData = new RooDataSet(dataName, dataName, *observables, 
			     RooFit::Index(*categories),
			     RooFit::Import(dataMap));    
  }
  
  workspace->import(*newData);
  
  // Save the workspace with the new dataset:
  TString outputFileName = Form("%s/workspace_%s.root", outputDir.Data(),
				(config->getStr("AnalysisType")).Data());
  workspace->writeToFile(outputFileName);
  std::cout << "AddDataToWorkspace: Workspace saved to " << outputFileName
	    << "\n\t New dataset named " << dataName << std::endl;

 
  delete workspace;
  workspaceFile->Close();
  delete workspaceFile;
  delete newData;
  delete m_treeMxAOD;
  delete chain;
  delete config;  
  return 0;
}
