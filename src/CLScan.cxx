////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CLScan.cxx                                                          //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/03/2016                                                          //
//                                                                            //
//  Performs a scan of the 95% CL using toy MC at various mass-width-cross-   //
//  section points.                                                           //
//                                                                            //
//  Macro options:                                                            //
//  - "New"        Calculate everything from scratch.                         //
//  - "FromFile"   Load CL values from file.                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CLScan.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the CLScan class.
   @param configFileName - The name of the analysis config file.
   @param options - Options for the CL scan.
*/
CLScan::CLScan(TString configFileName, TString options) {
  printer(Form("CLScan::CLScan(%s, %s)", configFileName.Data(), options.Data()),
	  false);
  
  // Load the config file:
  m_configFileName = configFileName;
  m_config = new Config(m_configFileName);
  m_options = options;
  
  // Set output directory:
  setInputDirectory(Form("%s/%s/GenericToys/single_files",
			 (m_config->getStr("MasterOutput")).Data(),
			 (m_config->getStr("JobName")).Data()));
  setOutputDirectory(Form("%s/%s/CLScan",
			  (m_config->getStr("MasterOutput")).Data(),
			  (m_config->getStr("JobName")).Data()));
  
  // Clear data:
  clearData();
  
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
}

/**
   -----------------------------------------------------------------------------
   Clear the class data.
*/
void CLScan::clearData() {
  m_massValues.clear();
  m_widthValues.clear();
  m_xsValues.clear();
  m_values95CL.clear();
}

/**
   -----------------------------------------------------------------------------
   Look at the available toy MC files to detect which mass, width, and cross-
   section values have associated pseudo-experiment ensembles.
   @param toyDirectory - The directory in which the toy MC files reside.
*/
void CLScan::detectMassWidthXSFiles(TString toyDirectory) {

  // First clear existing list:
  m_massValues.clear();
  m_widthValues.clear();
  m_xsValues.clear();
  
  // Make a list of toy MC files and loop over the list:
  system(Form("ls %s | tee tempToyList.txt", toyDirectory.Data()));
  std::ifstream toyFileList("tempToyList.txt");
  TString currToyName;
  while (toyFileList >> currToyName) {
    // Only look at toys made for scans:
    if (currToyName.Contains("ForScan")) {
      
      // Split up the file name into components:
      TObjArray *array = currToyName.Tokenize("_");
      for (int i_t = 0; i_t < array->GetEntries(); i_t++) {
	TString currElement = ((TObjString*)array->At(i_t))->GetString();
	// Token for mass value:
	if (currElement.Contains("mass")) {
	  TString massStr = currElement;
	  massStr.ReplaceAll("mass", "");
	  int massVal = massStr.Atoi();
	  if (!vectorContainsValue(m_massValues, massVal)) {
	    m_massValues.push_back(massVal);
	  }
	}
	// Token for width value:
	else if (currElement.Contains("width")) {
	  TString widthStr = currElement;
	  widthStr.ReplaceAll("width", "");
	  int widthVal = widthStr.Atoi();
	  if (!vectorContainsValue(m_widthValues, widthVal)) {
	    m_widthValues.push_back(widthVal);
	  }
	}
	// Token for cross-section value:
	else if (currElement.Contains("xs")) {
	  TString xsStr = currElement;
	  xsStr.ReplaceAll("xs", "");
	  xsStr.ReplaceAll(".root", "");
	  int xsVal = xsStr.Atoi();
	  if (!vectorContainsValue(m_xsValues, xsVal)) {
	    m_xsValues.push_back(xsVal);
	  }
	}
      }
    }
  }
  // Close and remove list of files:
  toyFileList.close();
  system("rm tempToyList.txt");
  
  // Then sort the vectors of unique mass, width, and cross-section values:
  std::sort(m_massValues.begin(), m_massValues.end());
  std::sort(m_widthValues.begin(), m_widthValues.end());
  std::sort(m_xsValues.begin(), m_xsValues.end());
}

/**
   -----------------------------------------------------------------------------
   Get the intersection point for the graph.
   @param graph - The graph for which intercept points will be found.
   @param valueToIntercept - The y-value to intercept.
   @return - The x-value of the intercept.
*/
double CLScan::getIntercept(TGraph *graph, double valueToIntercept) {
  
  // Loop over points in the graph to get search range:
  double rangeMin = 0.0; double rangeMax = 0.0;
  for (int i_p = 0; i_p < graph->GetN(); i_p++) {
    double xCurr; double yCurr;
    graph->GetPoint(i_p, xCurr, yCurr);
    if (i_p == 0) rangeMin = xCurr;
    if (i_p == (graph->GetN()-1)) rangeMax = xCurr;
  }
  
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

    if (currYValue > valueToIntercept) currXValue -= stepSize;
    else currXValue += stepSize;
  }
  
  // Print error message and return bad value if convergence not achieved:
  if (nIterations == maxIterations) {
    std::cout << "CLScan: ERROR! Intercept not found." << std::cout;
    return -999;
  } 
  return currXValue;
}

/**
   -----------------------------------------------------------------------------
   Get the limit which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limit (false for observed).
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @return - The 95% CL value of the cross-section.
*/
double CLScan::getLimit(int mass, int width, bool expected, int N) {
  TString key = expected ? Form("exp_mass%d_width%d_N%d", mass, width, N) : 
    Form("obs_mass%d_width%d_N%d", mass, width, N);
  if (m_values95CL.count(key) > 0) return m_values95CL[key];
  else printer(Form("CLScan: ERROR no value for key %s", key.Data()), true);
  return 0.0;
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void CLScan::printer(TString statement, bool isFatal) {
  if (m_config->getBool("Verbose") || isFatal) {
    std::cout << statement << std::endl;
  }
  if (isFatal) exit(0);
}

/**
   -----------------------------------------------------------------------------
   Plot the limits as a function of mass, for a given width slice.
   @param width - The width integer for the scan.
*/
void CLScan::scanMass(int width, bool makeNew) {
  
  // Arrays to store band information:
  double varValues[100] = {0};  
  double CLObs[100]     = {0};
  double CLExp_p2[100]  = {0};  
  double CLExp_p1[100]  = {0};
  double CLExp[100]     = {0};
  double CLExp_n1[100]  = {0};
  double CLExp_n2[100]  = {0};
  
  int nToyPoints = 0;

  // Loop over the mass points to load CL results
  for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
    if (singleCLScan(m_massValues[i_m], width, makeNew)) {
      CLObs[nToyPoints]    = getLimit(m_massValues[i_m], width, false, 0);
      CLExp_p2[nToyPoints] = getLimit(m_massValues[i_m], width, true, -2);
      CLExp_p1[nToyPoints] = getLimit(m_massValues[i_m], width, true, -1);
      CLExp[nToyPoints]    = getLimit(m_massValues[i_m], width, true, 0);
      CLExp_n1[nToyPoints] = getLimit(m_massValues[i_m], width, true, 1);
      CLExp_n2[nToyPoints]  = getLimit(m_massValues[i_m], width, true, 2);
      nToyPoints++;
    }
  }
  
  // Scan information:
  //std::vector<double> scanMXValues = config->getNumV("MXScanValues");
  
  //----------------------------------------//
  // Plot the results:
  double errExp_p2[100] = {0};  
  double errExp_p1[100] = {0};
  double errExp_n1[100] = {0};
  double errExp_n2[100] = {0};
  for (int i_t = 0; i_t < nToyPoints; i_t++) {
    errExp_p2[i_t] = fabs(CLExp_p2[i_t] - CLExp[i_t]);
    errExp_p1[i_t] = fabs(CLExp_p1[i_t] - CLExp[i_t]);
    errExp_n1[i_t] = fabs(CLExp_n1[i_t] - CLExp[i_t]);
    errExp_n2[i_t] = fabs(CLExp_n2[i_t] - CLExp[i_t]);
  }
  
  // Median expected and observed results:
  TGraph *gCLExp = new TGraph(nToyPoints, varValues, CLExp);
  TGraph *gCLObs = new TGraph(nToyPoints, varValues, CLObs);
  
  // Also plot the bands:
  TGraphAsymmErrors *gCLExp_2s
    = new TGraphAsymmErrors(nToyPoints, varValues, CLExp, 0, 0, 
 			    errExp_n2, errExp_p2);
  TGraphAsymmErrors *gCLExp_1s
    = new TGraphAsymmErrors(nToyPoints, varValues, CLExp, 0, 0, 
			    errExp_n1, errExp_p1);
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  gCLExp->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLObs->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLExp_2s->GetXaxis()->SetTitle("m_{X} [GeV]");
  gCLExp->GetYaxis()
    ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
  gCLObs->GetYaxis()
    ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
  gCLExp_2s->GetYaxis()
    ->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrowhh} [pb]");
  
  gCLExp->SetLineColor(kBlack);
  gCLObs->SetLineColor(kBlack);
  gCLExp->SetLineStyle(2);
  gCLObs->SetLineStyle(1);
  gCLExp->SetLineWidth(2);
  gCLObs->SetLineWidth(2);
  gCLExp_2s->SetFillColor(kYellow);
  gCLExp_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.61,0.68,0.89,0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!m_config->getBool("DoBlind")) leg.AddEntry(gCLObs,"Obs. limit","l");
  leg.AddEntry(gCLExp,"Exp. limit","l");
  leg.AddEntry(gCLExp_1s,"Exp. limit #pm1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_2s,"Exp. limit #pm2#sigma_{exp}","F");
  
  // Plotting options:
  gCLExp_2s->Draw("A3");
  gCLExp->Draw("Lsame");
  gCLExp_1s->Draw("3same");
  gCLExp->Draw("LSAME");
  if (!m_config->getBool("DoBlind")) gCLObs->Draw("LSAME");
  gPad->RedrawAxis();
  leg.Draw("SAME");
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.2, 0.87, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.87, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.2, 0.81, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			      (m_config->getNum("AnalysisLuminosity")/1000.0)));
  
  // Print the canvas:
  can->Print(Form("%s/limits_toy.eps", m_outputDir.Data()));
  
  // Delete pointers:
  printer(Form("CLScan: Finished mass scan for %d!", width), false);
  delete can;
  delete gCLObs;
  delete gCLExp;
  delete gCLExp_2s;
  delete gCLExp_1s;
}

/**
   -----------------------------------------------------------------------------
   Set the input file location.
   @param directory - The directory containing input files.
*/
void CLScan::setInputDirectory(TString directory) {
  m_inputDir = directory;
  // Also detect the input files in this new directory:
  detectMassWidthXSFiles(m_inputDir);
}

/**
   -----------------------------------------------------------------------------
   Set the limit which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limits (false for observed).
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @param value - The limit value to set.
*/
void CLScan::setLimit(int mass, int width, bool expected, int N, 
		      double limitValue) {
  if (expected) {
    m_values95CL[Form("exp_mass%d_width%d_N%d", mass, width, N)] = limitValue;
  }
  else {
    m_values95CL[Form("obs_mass%d_width%d_N%d", mass, width, N)] = limitValue;
  }
}

/**
   -----------------------------------------------------------------------------
   Set the output file location.
   @param directory - The directory for storing output files.
*/
void CLScan::setOutputDirectory(TString directory) {
  m_outputDir = directory;
  // Create output directory if it doesn't already exist:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
}

/**
   -----------------------------------------------------------------------------
   The main method scans the 95% CL for various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options: "New","FromFile","toy","asymptotic","NEvents"
   @param resMass - The resonance mass.
   @param makeNew - True if from toy MC, false if from text file storage.
   @return - True iff loaded successfully.
*/
bool CLScan::singleCLScan(int mass, int width, bool makeNew) {
  printer(Form("Form(singleCLScan(mass=%d, width=%d, makeNew=%d",
	       mass, width, (int)makeNew), false);
  
  // Arrays to store band information:
  int xsVals[100] = {0};
  double varValues[100] = {0};  
  double CLObs[100] = {0};
  double CLExp_p2[100] = {0};  
  double CLExp_p1[100] = {0};
  double CLExp[100] = {0};
  double CLExp_n1[100] = {0};
  double CLExp_n2[100] = {0};
  
  double qMuObs[100] = {0.0};
  double qMuExp[100] = {0.0};
  double qMuExp_p2[100] = {0.0};
  double qMuExp_p1[100] = {0.0};
  double qMuExp_n1[100] = {0.0};
  double qMuExp_n2[100] = {0.0};
  
  int nToyPoints = 0;
  
  // Scan information:
  //std::vector<double> scanValues = m_config->getNumV("CLScanValues");
  
  //----------------------------------------//
  // Open CL values from file:
  if (!makeNew) {
    
    // Open the saved CL values from toys:
    TString inputName = Form("%s/scan_CL_values.txt", m_outputDir.Data());
    std::ifstream inputFile(inputName);
    if (inputFile.is_open()) {
      while (inputFile >> varValues[nToyPoints] >> CLObs[nToyPoints] 
	     >> CLExp[nToyPoints] >> CLExp_p2[nToyPoints] 
	     >> CLExp_p1[nToyPoints] >> CLExp_n1[nToyPoints] 
	     >> CLExp_n2[nToyPoints]) {
	nToyPoints++;
      }
    }
    else {
      printer(Form("CLScan: ERROR opnening toy file %s", inputName.Data()), 
	      false);
    }
    inputFile.close();
  }
  
  //----------------------------------------//
  // Calculate value of qMu observed before processing toy files:
  else {
    // Open the workspace:
    TFile wsFile(m_config->getStr("WorkspaceFile"), "read");
    RooWorkspace *workspace
      = (RooWorkspace*)wsFile.Get(m_config->getStr("WorkspaceName"));
    
    // Instantiate the test statistic class for calculations and plots:
    TestStat *testStat = new TestStat(m_configFileName, "new", workspace);
    
    // Save values for plotting again!
    std::ofstream outFile(Form("%s/scan_CL_values.txt", m_outputDir.Data()));
    
    // Loop over cross-section:
    for (int i_x = 0; i_x < (int)m_xsValues.size(); i_x++) {
      double crossSection = ((double)m_xsValues[i_x])/1000.0;
      std::cout << "CLScan: cross-section = " << crossSection << std::endl;
      
      // Also force the mass and width to be constant always:
      testStat->setParam(m_config->getStr("PoIForMass"), 
			 (double)mass, true);
      testStat->setParam(m_config->getStr("PoIForWidth"), 
			 (((double)width)/100.0), true);
      
      // map of names and values of pois to set for fit.
      std::map<TString,double> mapPoIMu1; mapPoIMu1.clear();
      mapPoIMu1[m_config->getStr("PoIForNormalization")] = crossSection;
      mapPoIMu1[m_config->getStr("PoIForMass")] = (double)mass;
      mapPoIMu1[m_config->getStr("PoIForWidth")] = (((double)width)/100.0);
      
      // Perform the mu=1 fit and mu-free fit (necessary for qmu calculation):
      double nllObsMu1
	= testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
			      1, true, mapPoIMu1, false);
      double nllObsMuFree
	= testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
			      1, false, mapPoIMu1, false);
      
      // Get profiled signal strength from the mu-free fit:
      std::map<std::string,double> poiFromFit = testStat->getPoIs();
      double obsXSValue = poiFromFit[(std::string)m_config->getStr("PoIForNormalization")];
      double muHatValue = obsXSValue / crossSection;
      
      qMuObs[nToyPoints]
	= testStat->getQMuFromNLL(nllObsMu1, nllObsMuFree, muHatValue, 1);
      xsVals[nToyPoints] = m_xsValues[i_x];
      varValues[nToyPoints] = crossSection;
      nToyPoints++;
    }
    delete testStat;
    delete workspace;
    wsFile.Close();
    
    //----------------------------------------//
    // Process the toy MC files to obtain the CL values:
    for (int i_t = 0; i_t < nToyPoints; i_t++) {
      
      // Load the tool to analyze toys.
      // NOTE: this was moved outside the loop above because of interference
      // between the DHToyAnalysis class and DHTestStat, which is called
      // in DHToyAnalysis...
      ToyAnalysis *toyAna = new ToyAnalysis(m_configFileName, "None");
      toyAna->setOutputDir(m_outputDir);
      std::vector<TString> fitTypes; fitTypes.clear();
      fitTypes.push_back("0");
      fitTypes.push_back("1");
      fitTypes.push_back("Free");
      toyAna->setFitTypes(fitTypes);
      toyAna->loadToy(0, Form("%s/toy_mu0*ForScan_mass%d_width%d_xs%d.root",
			      m_outputDir.Data(), mass, width, xsVals[i_t]));
      toyAna->loadToy(1, Form("%s/toy_mu1*ForScan_mass%d_width%d_xs%d.root",
			      m_outputDir.Data(), mass, width, xsVals[i_t]));
      if (!(toyAna->areInputFilesOK())) {
	printer("CLScan: ERROR with toy scan option.", true);
      }
      
      // Calculate the expected qmu:
      qMuExp[i_t] = toyAna->calculateBkgQMuForN(0);
      qMuExp_p2[i_t] = toyAna->calculateBkgQMuForN(2.0);
      qMuExp_p1[i_t] = toyAna->calculateBkgQMuForN(1.0);
      qMuExp_n1[i_t] = toyAna->calculateBkgQMuForN(-1.0);
      qMuExp_n2[i_t] = toyAna->calculateBkgQMuForN(-2.0);
      
      CLObs[i_t] = toyAna->calculateCLFromToy(qMuObs[i_t]);
      CLExp[i_t] = toyAna->calculateCLFromToy(qMuExp[i_t]);
      CLExp_p2[i_t] = toyAna->calculateCLFromToy(qMuExp_p2[i_t]);
      CLExp_p1[i_t] = toyAna->calculateCLFromToy(qMuExp_p1[i_t]);
      CLExp_n1[i_t] = toyAna->calculateCLFromToy(qMuExp_n1[i_t]);
      CLExp_n2[i_t] = toyAna->calculateCLFromToy(qMuExp_n2[i_t]);
      
      // Write CL values to file:
      outFile << varValues[i_t] << " " << CLObs[i_t] << " " << CLExp[i_t] << " "
	      << CLExp_p2[i_t] << " " << CLExp_p1[i_t] << " " 
	      << CLExp_n1[i_t] << " " << CLExp_n2[i_t] << std::endl;
      delete toyAna;
    }
    
    // Close the files that save CL data:
    outFile.close();
  }
  
  //----------------------------------------//
  // Plot the CL scan results:
  double errExp_p2[100] = {0};  
  double errExp_p1[100] = {0};
  double errExp_n1[100] = {0};
  double errExp_n2[100] = {0};
  for (int i_t = 0; i_t < nToyPoints; i_t++) {
    errExp_p2[i_t] = fabs(CLExp_p2[i_t] - CLExp[i_t]);
    errExp_p1[i_t] = fabs(CLExp_p1[i_t] - CLExp[i_t]);
    errExp_n1[i_t] = fabs(CLExp_n1[i_t] - CLExp[i_t]);
    errExp_n2[i_t] = fabs(CLExp_n2[i_t] - CLExp[i_t]);
  }
  
  // Median expected and observed results:
  TGraph *gCLObs = new TGraph(nToyPoints, varValues, CLObs);
  TGraph *gCLExp = new TGraph(nToyPoints, varValues, CLExp);
  TGraph *gCLExp_p1 = new TGraph(nToyPoints, varValues, CLExp_p1);
  TGraph *gCLExp_p2 = new TGraph(nToyPoints, varValues, CLExp_p2);
  TGraph *gCLExp_n1 = new TGraph(nToyPoints, varValues, CLExp_n1);
  TGraph *gCLExp_n2 = new TGraph(nToyPoints, varValues, CLExp_n2);
  
  // Also plot the bands:
  TGraphAsymmErrors *gCLExp_2s 
    = new TGraphAsymmErrors(nToyPoints, varValues, CLExp, 0, 0, 
			    errExp_n2, errExp_p2);
  TGraphAsymmErrors *gCLExp_1s
    = new TGraphAsymmErrors(nToyPoints, varValues, CLExp, 0, 0, 
			    errExp_n1, errExp_p1);
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  gCLExp->GetXaxis()->SetTitle("#sigma#timesBR(G*#rightarrow#gamma#gamma [fb]");
  gCLObs->GetXaxis()->SetTitle("#sigma#timesBR(G*#rightarrow#gamma#gamma [fb]");
  gCLExp->GetYaxis()->SetTitle("CL_{s} value");
  gCLObs->GetYaxis()->SetTitle("CL_{s} value");
  gCLExp->SetLineColor(kBlack);
  gCLObs->SetLineColor(kBlack);
  gCLExp->SetLineStyle(2);
  gCLObs->SetLineStyle(1);
  gCLExp->SetLineWidth(2);
  gCLObs->SetLineWidth(2);
  gCLExp->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp_2s->SetFillColor(kYellow);
  gCLExp_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.64, 0.20, 0.89, 0.38);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!m_config->getBool("DoBlind")) leg.AddEntry(gCLObs, "Obs. CL_{s}", "l");
  leg.AddEntry(gCLExp, "Exp. CL_{s}", "l");
  leg.AddEntry(gCLExp_1s, "Exp. CL_{s} #pm1#sigma_{exp}", "F");
  leg.AddEntry(gCLExp_2s, "Exp. CL_{s} #pm2#sigma_{exp}", "F");
  
  // Plotting options:
  gCLExp->Draw("AL");
  gCLExp_2s->Draw("3same");
  gCLExp_1s->Draw("3same");
  gCLExp->Draw("LSAME");
  if (!m_config->getBool("DoBlind")) gCLObs->Draw("LSAME");
  gPad->RedrawAxis();
  leg.Draw("SAME");
  
  // 95% CL Line
  TLine *line = new TLine();
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  line->DrawLine(gCLExp->GetXaxis()->GetXmin(), 0.95,
		 gCLExp->GetXaxis()->GetXmax(), 0.95);
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.64, 0.48, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.76, 0.48, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.64, 0.42, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			       (m_config->getNum("AnalysisLuminosity")/1000.)));
  
  // Print the canvas:
  can->Print(Form("%s/scan_95CL_mass%d_width%d.eps",
		  m_outputDir.Data(), mass, width));
  
  // Get the actual limit values:
  double observedCL = getIntercept(gCLObs, 0.95);
  double expectedCL = getIntercept(gCLExp, 0.95);
  double expectedCL_p1 = getIntercept(gCLExp_p1, 0.95);
  double expectedCL_p2 = getIntercept(gCLExp_p2, 0.95);
  double expectedCL_n1 = getIntercept(gCLExp_n1, 0.95);
  double expectedCL_n2 = getIntercept(gCLExp_n2, 0.95);  
  
  // Store the limits for this mass and width:
  setLimit(mass, width, false, 0, observedCL);
  setLimit(mass, width, true, 0, expectedCL);
  setLimit(mass, width, true, 2, expectedCL_p1);
  setLimit(mass, width, true, 1, expectedCL_p2);
  setLimit(mass, width, true, -1, expectedCL_n1);
  setLimit(mass, width, true, -2, expectedCL_n2);
  
  // Then print to screen:
  std::cout << "\nCLScan: Results" << std::endl;
  std::cout << "\tobs. CL " << observedCL << std::endl;
  std::cout << "\texp. CL " << expectedCL << std::endl;
  std::cout << "\texp. CL +1 " << expectedCL_p1 << std::endl;
  std::cout << "\texp. CL +2 " << expectedCL_p2 << std::endl;
  std::cout << "\texp. CL -1 " << expectedCL_n1 << std::endl;
  std::cout << "\texp. CL -2 " << expectedCL_n2 << std::endl;
  
  // Delete pointers, close files, return:
  printer(Form("CLScan: Finished mass %d width %d!\n", mass, width), false);
  
  delete line;
  delete can;
  delete gCLObs;
  delete gCLExp;
  delete gCLExp_2s;
  delete gCLExp_1s;
  
  // Return true iff. calculation was successful:
  if (nToyPoints < 2 || observedCL < 0 || expectedCL < 0) return false;
  else return true;
}

/**
   -----------------------------------------------------------------------------
   Check if a value belongs to a vector.
   @param theVector - The vector that might contain the value.
   @param theValue - The value to check for membership in the vector.
   @return - True iff the vector already contains the value.
*/
bool CLScan::vectorContainsValue(std::vector<int> theVector, int theValue) {
  for (int i_v = 0; i_v < (int)theVector.size(); i_v++) {
    if (theVector[i_v] == theValue) return true;
  }
  return false;
}
