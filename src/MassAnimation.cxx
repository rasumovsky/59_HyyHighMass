////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: MassAnimation.cxx                                                   //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 06/04/2016                                                          //
//                                                                            //
//  Creates an animated diphoton mass specturm .gif.                          //
//                                                                            //
//  Procedure:                                                                //
//    animation->setOutputDirectory                                           //
//    animation->setNFrames(200);                                             //
//    animation->getDataForFrames();                                          //
//    in loop, do animation->makeSingleFrame(i_f);                            //
//    animation->makeGIF();                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "MassAnimation.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the MassAnimation class.
   @param configFileName - The name of the analysis config file.
   @param options - Options for the CL scan.
*/
MassAnimation::MassAnimation(TString configFileName, TString options) {
    
  // Load the config file:
  m_configFileName = configFileName;
  m_config = new Config(m_configFileName);
  m_options = options;

  printer(Form("MassAnimation::MassAnimation(%s, %s)", configFileName.Data(),
	       options.Data()), false);

  // Set output directory:
  setOutputDirectory(Form("%s/%s/MassAnimation",
			  (m_config->getStr("MasterOutput")).Data(),
			  (m_config->getStr("JobName")).Data()));
  
  // Load the workspace:
  // See if a local copy of the workspace has been made:
  TString workspaceFileName = m_config->getStr("WorkspaceFile");
  TObjArray *array = workspaceFileName.Tokenize("/");

  workspaceFileName
    = ((TObjString*)array->At(array->GetEntries()-1))->GetString();
  m_workspaceFile = new TFile(workspaceFileName);
  //if (m_workspaceFile->IsOpen()) {
  if (!m_workspaceFile->IsZombie()) {
    printer(Form("MassAnimation: Loaded ws from %s", workspaceFileName.Data()),
	    false);
  }
  else {
    workspaceFileName = m_config->getStr("WorkspaceFile");
    printer(Form("MassAnimation: Load ws from %s", workspaceFileName.Data()),
	    false);
    m_workspaceFile = new TFile(workspaceFileName, "read");
  }
  
  //m_workspaceFile = new TFile(m_config->getStr("WorkspaceFile"), "read");
  m_workspace = (RooWorkspace*)m_workspaceFile
    ->Get(m_config->getStr("WorkspaceName"));
  m_model = (ModelConfig*)m_workspace
    ->obj(m_config->getStr("WorkspaceModelConfig"));
  m_combPdf = (RooSimultaneous*)m_model->GetPdf();
  m_observables = (RooArgSet*)m_model->GetObservables();
  TString nameRooCategory = m_config->getStr("WorkspaceRooCategory");
  if (m_workspace->obj(nameRooCategory)) {
    m_categories = (RooCategory*)m_workspace->obj(nameRooCategory);
  }
  else {
    printer(Form("MassAnimation: RooCategory object %s not found",
		 nameRooCategory.Data()), true);
  }
  
  // Get observable and PDF names:  
  m_obsName = "";
  m_pdfName = "";
  TIterator *cateIter = m_combPdf->indexCat().typeIterator();
  RooCatType *cateType = NULL;
  while ((cateType = (RooCatType*)cateIter->Next())) {
    RooAbsPdf *currPdf = m_combPdf->getPdf(cateType->GetName());
    m_obsName = currPdf->getObservables(m_observables)->first()->GetName();
    m_pdfName = currPdf->GetName();
    break;
  }

  // Set the range for the observable:
  std::vector<double> massRangeVector
    = m_config->getNumV(Form("PoIRange_%s",
			     (m_config->getStr("PoIForMass")).Data()));
  m_workspace->var(m_obsName)->setMin(massRangeVector[0]);
  m_workspace->var(m_obsName)->setMax(massRangeVector[1]);
  m_workspace->var(m_obsName)->setRange("fullRange", massRangeVector[0], 
					massRangeVector[1]);
  

  // Settings for this class:
  m_geVPerBin = m_config->getInt("AnimationGeVPerBin");
  m_nFrames = 0;
  m_times.clear();
  
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
}

/**
   -----------------------------------------------------------------------------
   Load the data from an MxAOD for all frames that have been set using the
   MassAnimation::setNFrames() method.
*/
void MassAnimation::getDataForFrames() {
  
  // Instantiate the RooDataSets:
  std::map<int,std::map<std::string,RooDataSet*> > mapOfDataMaps;
  mapOfDataMaps.clear();
  std::string cate1 = "";
  RooRealVar *obs1 = NULL;
  for (int i_f = 0; i_f < m_nFrames; i_f++) {
    mapOfDataMaps[i_f].clear();
    
    // Iterate over categories:
    TIterator *cateIter = m_combPdf->indexCat().typeIterator();
    RooCatType *cateType = NULL;
    while ((cateType = (RooCatType*)cateIter->Next())) {
      RooAbsPdf *currPdf = m_combPdf->getPdf(cateType->GetName());
      RooArgSet *currObs = currPdf->getObservables(m_observables);
      // Define the (inclusive) category name and observable variable:
      if (cate1 =="") cate1 = (std::string)(cateType->GetName());
      if (obs1 == NULL) obs1 = (RooRealVar*)currObs->first();
      // Define a new dataset:
      mapOfDataMaps[i_f][(std::string)cateType->GetName()]
	= new RooDataSet(Form("data_%d_%s", i_f,
			      TString(cateType->GetName()).Data()),
			 Form("data_%d_%s", i_f,
			      TString(cateType->GetName()).Data()),
			 RooArgSet(*currObs));
    }
  }
  
  // Prepare for loop over input MxAOD/TTree:
  std::vector<TString> fileNames = m_config->getStrV("MxAODsForData");
  
  // Make local copies of files if requested, to improve speed:
  if (m_config->getBool("MakeLocalMxAODCopies")) {
    fileNames = makeLocalMxAODCopies(fileNames);
  }
  // Create TChain of input files:
  TChain *chain = new TChain(m_config->getStr("MxAODTreeName"));
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    chain->AddFile(fileNames[i_f]);
  }
  HGammaMxAOD *treeMxAOD = new HGammaMxAOD(chain, m_config->getStr("MxAODTag"));
  
  //--------------------------------------//
  // Loop over events to build dataset for signal parameterization:
  //int nEvents = (int)(0.1 * treeMxAOD->fChain->GetEntries());
  int nEvents = (int)(treeMxAOD->fChain->GetEntries());
  std::cout << "There are " << nEvents << " events to process." << std::endl;
  for (int index = 0; index < nEvents; index++) {

    // Load event from MxAOD:
    treeMxAOD->fChain->GetEntry(index);
    printProgressBar(index, nEvents);

    // Record the time stamp:
    for (int i_f = 0; i_f < m_nFrames; i_f++) {
      if (((((double)index+1)/((double)nEvents)) >= 
	   (((double)i_f+1)/((double)m_nFrames))) &&
	  ((((double)index)/((double)nEvents)) <
	   (((double)i_f+1)/((double)m_nFrames)))) {
	m_times[i_f]=timeToString((time_t)(treeMxAOD->EventInfoAux_timeStamp));
	printer(Form("MassAnimation::getDataForFrames: Setting timestamp %d", 
		     treeMxAOD->EventInfoAux_timeStamp), false);
      }
    }
    
    //---------- Event selection ----------//
    if ((m_config->getStr("AnalysisType")).Contains("Graviton") &&
	!(treeMxAOD->HGamEventInfoAuxDyn_isPassedExotic && 
	  treeMxAOD->HGamEventInfoAuxDyn_isPassedIsolationLowHighMyy &&
	  treeMxAOD->HGamEventInfoAuxDyn_m_yy > 200000)) {
      continue;
    }
    else if ((m_config->getStr("AnalysisType")).EqualTo("Scalar") &&
	     !(treeMxAOD->HGamEventInfoAuxDyn_isPassedLowHighMyy && 
	       treeMxAOD->HGamEventInfoAuxDyn_m_yy > 150000)) {
      continue;
    }
    
    // Set the mass observable:
    obs1->setVal(treeMxAOD->HGamEventInfoAuxDyn_m_yy/1000.0);
    
    // Fill the RooDataSets for each frame:
    for (int i_f = 0; i_f < m_nFrames; i_f++) {
      if ((((double)index+1)/((double)nEvents)) < 
	  (((double)i_f+1)/((double)m_nFrames))) {
	mapOfDataMaps[i_f][cate1]->add(RooArgSet(*obs1));
      }
    }
  }// End of event loop.
  
  // Remove files that were copied:
  if (m_config->getBool("MakeLocalMxAODCopies")) {
    removeLocalMxAODCopies(fileNames);
  }
  delete treeMxAOD;
  delete chain;
    
  // Loop over the frames, creating combined datasets for each:
  for (int i_f = 0; i_f < m_nFrames; i_f++) {
    m_data[i_f] = new RooDataSet(Form("dataForFrame%d",i_f),
				 Form("dataForFrame%d",i_f), *m_observables, 
				 RooFit::Index(*m_categories),
				 RooFit::Import(mapOfDataMaps[i_f]));
  }
}

/**
   -----------------------------------------------------------------------------
   Get the number of frames that are used in the GIF.
   @return - number of GIF frames.
*/
int MassAnimation::getNFrames() {
  return m_nFrames;
}

/**
   -----------------------------------------------------------------------------
   Make a GIF of the data by loading a canvas saved for each frame.
*/
void MassAnimation::makeGIF() {
  printer("MassAnimation::makeGIF()", false);
  

  // Remove prior animations:
  system(Form("rm %s/mass_animation.gif", m_outputDir.Data()));
  
  // Loop over the frames:
  for (int i_f = 0; i_f < m_nFrames; i_f++) {
    TFile currFile(Form("%s/file_frame%d.root", m_outputDir.Data(), i_f));
    if (currFile.IsOpen()) {
      TCanvas *currCan = (TCanvas*)currFile.Get("can");
      if (i_f == m_nFrames - 1) {
	currCan->Print(Form("%s/mass_animation.gif+500", m_outputDir.Data()));
	currCan->Print(Form("%s/mass_animation.gif++", m_outputDir.Data()));
      }
      else {
	currCan->Print(Form("%s/mass_animation.gif+10", m_outputDir.Data()));
      }
      currFile.Close();
      delete currCan;
    }
  }
  std::cout << "\tGIF printed to " << m_outputDir << "/mass_animation.gif"
	    << std::endl;
}

/**
   -----------------------------------------------------------------------------
   Loops over the total number of frames to create them all.
*/
void MassAnimation::makeAllFrames() {
  printer("MassAnimation::makeAllFrames()", false);
  
  for (int i_f = 0; i_f < m_nFrames; i_f++) makeSingleFrame(i_f);
}

/**
   -----------------------------------------------------------------------------
   Make a single frame of the animation (saves a canvas).
   @param frame - The index of the current frame.
*/
void MassAnimation::makeSingleFrame(int frame) {
  printer(Form("MassAnimation::makeSingleFrame(%d)",frame), false);

  // Create a canvas with two pads (one main plot, one subtraction plot)
  TCanvas *can = new TCanvas("can", "can", 800, 800);
  can->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.00, 0.33, 1.00, 1.00);
  TPad *pad2 = new TPad("pad2", "pad2", 0.00, 0.00, 1.00, 0.33);
  pad1->SetBottomMargin(0.00001);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.4);
  pad2->SetBorderMode(0);
  can->cd();
  pad1->Draw();
  pad2->Draw();
  
  //---------- Pad 1: The mass points ----------//
  pad1->cd();
  
  // Import the current dataset:
  TString dataName = m_data[frame]->GetName();
  m_workspace->import(*m_data[frame]);
  
  // Load background-only fit snapshot if available, otherwise do fit:
  TString snapshotNameMu0 = m_config->getStr("WorkspaceSnapshotMu0");
  if (m_workspace->getSnapshot(snapshotNameMu0)) {
    m_workspace->loadSnapshot(snapshotNameMu0);
  }
  else {
        
    // Background-only fit to the current dataset:
    std::vector<TString> namesPoI = m_config->getStrV("WorkspacePoIs");
    for (int i_p = 0; i_p < (int)namesPoI.size(); i_p++) {
      if ((namesPoI[i_p]).EqualTo(m_config->getStr("PoIForNormalization"))) {
	m_workspace->var(namesPoI[i_p])->setVal(0.0);
      }
      m_workspace->var(namesPoI[i_p])->setConstant(true);
    }
    
    // Turn off statistical errors of template:
    if (m_config->isDefined("TurnOffTemplateStat") && 
	m_config->getBool("TurnOffTemplateStat")) {
      RooArgSet *set_temp=new RooArgSet();
      set_temp->add(*(RooArgSet*)m_model->GetNuisanceParameters()
		    ->selectByName("*gamma_stat*"));
      statistics::constSet(set_temp, true);
    }
    
    // The actual fit command:
    RooNLLVar* varNLL
      = (RooNLLVar*)m_combPdf->createNLL(*m_workspace->data(dataName),
					 Extended(m_combPdf->canBeExtended()));
    RooFitResult *fitResult = statistics::minimize(varNLL, "", NULL, true);
  }
  
  double rMin = m_workspace->var(m_obsName)->getMin();
  double rMax = m_workspace->var(m_obsName)->getMax();
  int rBins = (int)((rMax - rMin) / (double)m_geVPerBin);
  
  // Plot the frame:
  RooPlot *rooPlot = m_workspace->var(m_obsName)
    ->frame(RooFit::Bins(rBins),RooFit::Range(rMin,rMax));
  rooPlot->SetYTitle(Form("Events / %d GeV", m_geVPerBin));
  rooPlot->SetXTitle("m_{#gamma#gamma} [GeV]");
  
  // Then add the data and PDF to the RooPlot:
  if (m_workspace->data(dataName)) {
    m_workspace->data(dataName)->plotOn(rooPlot);
  }
  else printer(Form("MassAnimation:makeSingleFrame: Missing dataset %s", 
		    dataName.Data()), true);
  
  if (m_workspace->pdf(m_pdfName)) {
    m_workspace->pdf(m_pdfName)->plotOn(rooPlot, RooFit::Range("fullRange"), 
				      RooFit::NormRange("fullRange"), 
				      RooFit::LineColor(kBlue+1));
  }
  else printer(Form("MassAnimation:makeSingleFrame: Missing PDF %s", 
		    m_pdfName.Data()), true);
  
  // Draw the RooPlot:
  rooPlot->Draw();
  // Special y-axis ranges for log-scale plots:
  gPad->SetLogy();
  
  // Plot text (ATLAS, sqrt(s), lumi, and time):
  TLatex l; l.SetNDC(); l.SetTextColor(kBlack);
  l.SetTextFont(72); l.SetTextSize(0.05); 
  l.DrawLatex(0.60, 0.88, "ATLAS");
  l.SetTextFont(42); l.SetTextSize(0.05); 
  l.DrawLatex(0.72, 0.88, m_config->getStr("ATLASLabel"));
  double frameLumi = (((double)(frame+1) / (double)m_nFrames) * 
		      m_config->getNum("AnalysisLuminosity") / 1000.0);
  l.DrawLatex(0.6, 0.82, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}", frameLumi));
  l.DrawLatex(0.6, 0.76, m_times[frame]);
  
  //---------- Pad 2: The subtraction plot ----------//
  pad2->cd();
  
  TGraphErrors* subData = plotComparison(m_workspace->data(dataName),
					 m_workspace->pdf(m_pdfName),
					 m_workspace->var(m_obsName), rBins);
  
  TH1F *medianHist = new TH1F("median", "median", rBins, rMin, rMax);
  for (int i_b = 1; i_b <= rBins; i_b++) medianHist->SetBinContent(i_b, 1.0);
  medianHist->GetYaxis()->SetRangeUser(0.0, 2.0);
  medianHist->SetLineColor(kBlue+1);
  medianHist->SetLineWidth(2);
  medianHist->GetYaxis()->SetTitle("Data / Fit");
  medianHist->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  medianHist->GetYaxis()->SetNdivisions(5);
  medianHist->GetXaxis()->SetTitleOffset(0.95);
  medianHist->GetYaxis()->SetTitleOffset(0.7);
  medianHist->GetXaxis()->SetTitleSize(0.1);
  medianHist->GetYaxis()->SetTitleSize(0.1);
  medianHist->GetXaxis()->SetLabelSize(0.1);
  medianHist->GetYaxis()->SetLabelSize(0.1);
  medianHist->Draw();
  subData->Draw("EPSAME");
  
  // Print and also save to file:
  can->Print(Form("%s/plot_mass_frame%d.eps", m_outputDir.Data(), frame));
  TFile *outputFile = new TFile(Form("%s/file_frame%d.root",
				     m_outputDir.Data(), frame), "RECREATE");
  
  can->Write();
  outputFile->Close();
  
  //delete outputFile;
  //delete can;
  //delete pad1;
  //delete pad2; 
}

/**
   -----------------------------------------------------------------------------
   Create a subtraction or ratio plot.
   @param data - The RooAbsData set for comparison.
   @param pdf - The PDF for comparison.
   @param observable - The mass observable for data and pdf. 
   @param xBins - The number of bins for the observable.
   @return - A TGraphErrors to plot.
*/
TGraphErrors* MassAnimation::plotComparison(RooAbsData *data, RooAbsPdf *pdf, 
					    RooRealVar *observable, 
					    double xBins, bool doRatio) {
  double minOrigin = observable->getMin();
  double maxOrigin = observable->getMax();
  double nEvents = data->sumEntries(Form("%s>%f&&%s<%f",
					 observable->GetName(), minOrigin,
					 observable->GetName(), maxOrigin));
  TH1F *originHist
    = (TH1F*)data->createHistogram("dataSub", *observable,
				   RooFit::Binning(xBins,minOrigin,maxOrigin));
  TGraphErrors *result = new TGraphErrors();
  double increment = (maxOrigin - minOrigin) / ((double)xBins);
  
  RooAbsReal* intTot
    = (RooAbsReal*)pdf->createIntegral(RooArgSet(*observable),
				       RooFit::NormSet(*observable), 
				       RooFit::Range("fullRange"));
  double valTot = intTot->getVal();
  
  int pointIndex = 0;
  for (double i_m = minOrigin; i_m < maxOrigin; i_m += increment) {
    //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
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
    double currWeight = doRatio ? 
      (currDataWeight / currPdfWeight) : (currDataWeight - currPdfWeight);
    double currError = doRatio ? 
      (originHist->GetBinError(pointIndex+1) / currPdfWeight) :
      originHist->GetBinError(pointIndex+1);
    if (doRatio && currPdfWeight == 0.0) {
      currWeight = 0.0;
      currError = 0.0;
    }
    result->SetPoint(pointIndex, currMass, currWeight);
    result->SetPointError(pointIndex, 0.0, currError);
    pointIndex++;
  }
  observable->setMin(minOrigin);
  observable->setMax(maxOrigin);
  delete originHist;
  return result;
}

/**
   -----------------------------------------------------------------------------
   Copy files from a slow resource (e.g. EOS) to the local disk for faster
   processing.
   @param fileNames - The original file names.
   @return - An updated list of file names.
*/
std::vector<TString> MassAnimation::makeLocalMxAODCopies(std::vector<TString>
							 fileNames) {
  std::cout << "MassAnimation: Making local copies of inputs."
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
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void MassAnimation::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Prints a progress bar to screen to provide elapsed time and remaining time
   information to the user. This is useful when processing large datasets. 
   @param index - The current event index.
   @param total - The total number of events.
*/
void MassAnimation::printProgressBar(int index, int total) {
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
   Remove any files that were copied over for speed.
   @param fileNames - The original file names.
*/
void MassAnimation::removeLocalMxAODCopies(std::vector<TString> fileNames) {
  std::cout << "MassAnimation: Removing local copies of inputs."
	    << std::endl;
  for (int i_f = 0; i_f < (int)fileNames.size(); i_f++) {
    system(Form("rm %s", fileNames[i_f].Data()));
  }
}

/**
   -----------------------------------------------------------------------------
   Set the number of frames to use for the GIF.
   @param nFrames - The number of frames.
*/
void MassAnimation::setNFrames(int nFrames) {
  m_nFrames= nFrames;
}

/**
   -----------------------------------------------------------------------------
   Set the output file location.
   @param directory - The directory for storing output files.
*/
void MassAnimation::setOutputDirectory(TString directory) {
  m_outputDir = directory;
  // Create output directory if it doesn't already exist:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
}

/**
   -----------------------------------------------------------------------------
   Convert a UNIX time into a legible date string. 
   @param timeInt a long integer corresponding to the UNIX time.
   @return - An ASC time string.
*/
TString MassAnimation::timeToString(time_t dateValue) {
  TString dateName = TString(asctime(gmtime(&dateValue)));
  TString days[7] = {"Fri ","Sat ","Sun ","Mon ","Tue ","Wed ","Thu "};
  for (int i_d = 0; i_d < 7; i_d++) dateName.ReplaceAll(days[i_d], "");
  TString year = dateName.Contains("2015") ? " 2015" : " 2016";
  dateName.Resize(6);
  dateName += year;
  dateName.ReplaceAll( "  ", " " );
  dateName.ReplaceAll( " 201", ", 201" );
  return dateName;
}
