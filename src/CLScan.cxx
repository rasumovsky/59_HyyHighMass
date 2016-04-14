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
////////////////////////////////////////////////////////////////////////////////

#include "CLScan.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the CLScan class.
   @param configFileName - The name of the analysis config file.
   @param options - Options for the CL scan.
*/
CLScan::CLScan(TString configFileName, TString options) {
    
  // Load the config file:
  m_configFileName = configFileName;
  m_config = new Config(m_configFileName);
  m_options = options;

  printer(Form("CLScan::CLScan(%s, %s)", configFileName.Data(), options.Data()),
	  false);

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
  m_valuesP0.clear();
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
   Get the p0 which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limit (false for observed).
   @return - The p0 value.
*/
double CLScan::getP0(int mass, int width, bool expected) {
  TString key = expected ? Form("exp_mass%d_width%d", mass, width) : 
    Form("obs_mass%d_width%d", mass, width);
  if (m_valuesP0.count(key) > 0) return m_valuesP0[key];
  else printer(Form("CLScan: ERROR no value for key %s", key.Data()), true);
  return 0.0;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the mass values stored in the scan tool.
*/
std::vector<int> CLScan::listMasses() {
  return m_massValues;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the width values stored in the scan tool.
*/
std::vector<int> CLScan::listWidths() {
  return m_widthValues;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the cross-section values stored in the scan tool.
*/
std::vector<int> CLScan::listXS() {
  return m_xsValues;
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
   @param makeNew - True iff. calculating things from toys directly.
*/
void CLScan::scanMassLimit(int width, bool makeNew) {
  
  TString limitFileName = Form("%s/limit_values_width%d.txt",
			       m_outputDir.Data(), width);
  
  // Arrays to store band information:
  double varValues[100] = {0};  
  double CLObs[100]     = {0};
  double CLExp_p2[100]  = {0};  
  double CLExp_p1[100]  = {0};
  double CLExp[100]     = {0};
  double CLExp_n1[100]  = {0};
  double CLExp_n2[100]  = {0};
  int nToyPoints = 0;
  
  //----------------------------------------//
  // Access results either by loading toy MC or stored limit computations:
  // Loop over the mass points to load CL results
  if (makeNew || m_options.Contains("ForceScan")) {
    printer("CLScan: Calculating limits vs. mass from scratch", false);
    
    std::ofstream limitFileOut(limitFileName);
    for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
      if (singleCLScan(m_massValues[i_m], width, makeNew)) {
	varValues[nToyPoints] = m_massValues[i_m];
	CLObs[nToyPoints] = getLimit(m_massValues[i_m], width, false, 0);
	CLExp_p2[nToyPoints] = getLimit(m_massValues[i_m], width, true, -2);
	CLExp_p1[nToyPoints] = getLimit(m_massValues[i_m], width, true, -1);
	CLExp[nToyPoints] = getLimit(m_massValues[i_m], width, true, 0);
	CLExp_n1[nToyPoints] = getLimit(m_massValues[i_m], width, true, 1);
	CLExp_n2[nToyPoints] = getLimit(m_massValues[i_m], width, true, 2);
	limitFileOut << m_massValues[i_m] << " " << CLObs[nToyPoints] << " " 
		  << CLExp_p2[nToyPoints] << " " << CLExp_p1[nToyPoints] << " " 
		  << CLExp[nToyPoints] << " " << CLExp_n1[nToyPoints] << " " 
		  << CLExp_n2[nToyPoints] << " " << std::endl;
	nToyPoints++;
      }
    }
    limitFileOut.close();
  }
  // Or load from text file:
  else {
    printer("CLScan: Loading limits vs. mass from file.", false);
    
    double currMass = 0.0;
    std::ifstream limitFileIn(limitFileName);
    if (limitFileIn.is_open()) {
      while (limitFileIn >> currMass >> CLObs[nToyPoints]
	     >> CLExp_p2[nToyPoints] >> CLExp_p1[nToyPoints] 
	     >> CLExp[nToyPoints] >> CLExp_n1[nToyPoints] 
	     >> CLExp_n2[nToyPoints]) {
	varValues[nToyPoints] = currMass;
        setLimit(currMass, width, false, 0, CLObs[nToyPoints]);
	setLimit(currMass, width, true, -2, CLExp_p2[nToyPoints]);
	setLimit(currMass, width, true, -1, CLExp_p1[nToyPoints]);
	setLimit(currMass, width, true, 0, CLExp[nToyPoints]);
	setLimit(currMass, width, true, 1, CLExp_n1[nToyPoints]);
	setLimit(currMass, width, true, 2, CLExp_n2[nToyPoints]);
	
	std::cout << currMass << " " << CLObs[nToyPoints] << " " 
		  << CLExp_p2[nToyPoints] << " " << CLExp_p1[nToyPoints] << " " 
		  << CLExp[nToyPoints] << " " << CLExp_n1[nToyPoints] << " " 
		  << CLExp_n2[nToyPoints] << " " << std::endl;
	nToyPoints++;
      }
    }
    limitFileIn.close();
  }
  
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
  gCLExp->GetXaxis()->SetTitle("m_{G*} [GeV]");
  gCLObs->GetXaxis()->SetTitle("m_{G*} [GeV]");
  gCLExp_2s->GetXaxis()->SetTitle("m_{G*} [GeV]");
  gCLExp->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [pb]");
  gCLObs->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [pb]");
  gCLExp_2s->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [pb]");
  
  gCLExp->SetLineColor(kBlack);
  gCLObs->SetLineColor(kBlack);
  gCLExp->SetLineStyle(2);
  gCLObs->SetLineStyle(1);
  gCLExp->SetLineWidth(2);
  gCLObs->SetLineWidth(2);
  gCLExp_2s->SetFillColor(kYellow);
  gCLExp_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.61, 0.68, 0.89, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!m_config->getBool("DoBlind")) {
    leg.AddEntry(gCLObs,"Observed #it{CL_{s}} limit","l");
  }
  leg.AddEntry(gCLExp,"Expected #it{CL_{s}} limit","l");
  leg.AddEntry(gCLExp_1s,"Expected #pm 1#sigma_{exp}","F");
  leg.AddEntry(gCLExp_2s,"Expected #pm 2#sigma_{exp}","F");
  
  // Plotting options:
  gCLExp_2s->GetXaxis()->SetRangeUser(m_massValues[0], 
				      m_massValues[m_massValues.size()-1]);
  gCLExp_2s->GetYaxis()->SetRangeUser(0.1, 1000);
  
  gPad->SetLogy();
  gCLExp_2s->Draw("A3");
  gCLExp->Draw("Lsame");
  gCLExp_1s->Draw("3same");
  gCLExp->Draw("LSAME");
  if (!m_config->getBool("DoBlind")) gCLObs->Draw("PSAME");
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
  t.DrawLatex(0.2, 0.75, "Spin-2 Selection");
  t.DrawLatex(0.2, 0.69,
	      Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{PI}=%2.2f",
		   ((double)width)/100.0));
  
  // Print the canvas:
  can->Print(Form("%s/limits_width%d.eps", m_outputDir.Data(), width));
  gPad->SetLogy(false);
  
  // Save graphs to file:
  TFile *outLimitFile = new TFile(Form("%s/limit_graphs_width%d.root",
				       m_outputDir.Data(), width), "RECREATE");
  gCLObs->Write();
  gCLExp->Write();
  gCLExp_2s->Write();
  gCLExp_1s->Write();
  can->Write();
  outLimitFile->Close();
  
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
   Plot the p0 value as a function of mass, for a given width slice.
   @param width - The width integer for the scan.
   @param makeNew - True iff. calculating things from toys directly.
*/
void CLScan::scanMassP0(int width, bool makeNew) {
  
  TString p0FileName = Form("%s/p0_values_width%d.txt",
			    m_outputDir.Data(), width);
  
  // Arrays to store band information:
  double varValues[100] = {0};  
  double p0Obs[100] = {0};
  double p0Exp[100] = {0};
  int nToyPoints = 0;
  
  //----------------------------------------//
  // Access results either by loading toy MC or stored limit computations:
  // Loop over the mass points to load CL results
  if (makeNew || m_options.Contains("ForceScan")) {
    printer("CLScan: Calculating limits vs. mass from scratch", false);
    
    std::ofstream p0FileOut(p0FileName);
    for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
      if (singleP0Test(m_massValues[i_m], width, makeNew)) {
	varValues[nToyPoints] = m_massValues[i_m];
	p0Obs[nToyPoints] = getP0(m_massValues[i_m], width, false);
	//p0Exp[nToyPoints] = getP0(m_massValues[i_m], width, true);
	p0Exp[nToyPoints] = 0.5;// CURRENTLY NOT IMPLEMENTED
	// ISSUE IS: WHAT XSECTION TO EXPECT IN TOY ANALYSIS?
	p0FileOut << m_massValues[i_m] << " " << p0Obs[nToyPoints] << " " 
		 << p0Exp[nToyPoints] << std::endl;
	nToyPoints++;
      }
    }
    p0FileOut.close();
  }
  // Or load from text file:
  else {
    printer("CLScan: Loading p0 vs. mass from file.", false);
    
    double currMass = 0.0;
    std::ifstream p0FileIn(p0FileName);
    if (p0FileIn.is_open()) {
      while (p0FileIn >> currMass >> p0Obs[nToyPoints] >> p0Exp[nToyPoints]) {
	varValues[nToyPoints] = currMass;
	setP0(currMass, width, false, p0Obs[nToyPoints]);
	setP0(currMass, width, true, p0Exp[nToyPoints]);
	std::cout << currMass << " " << p0Obs[nToyPoints] << " " 
		  << p0Exp[nToyPoints] << std::endl;
	nToyPoints++;
      }
    }
    p0FileIn.close();
  }
  
  //----------------------------------------//
  // Plot the results:
    
  // Median expected and observed results:
  //TGraph *gP0Exp = new TGraph(nToyPoints, varValues, p0Exp);
  TGraph *gP0Obs = new TGraph(nToyPoints, varValues, p0Obs);
    
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  /*
  gP0Exp->GetXaxis()->SetTitle("m_{G*} [GeV]");
  gP0Exp->GetYaxis()->SetTitle("p_{0}");
  gP0Exp->SetLineColor(kBlack);
  gP0Exp->SetLineStyle(2);
  gP0Exp->SetLineWidth(2);
  */
  
  gP0Obs->GetXaxis()->SetTitle("m_{G*} [GeV]");
  gP0Obs->GetYaxis()->SetTitle("p_{0}");  
  gP0Obs->SetLineColor(kBlack);
  gP0Obs->SetLineStyle(1);
  gP0Obs->SetLineWidth(2);
  
  // Legend:
  /*
  TLegend leg(0.61, 0.68, 0.89, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (m_config->getBool("DoBlind")) leg.AddEntry(gP0Exp,"Expected p_{0}","P");
  else leg.AddEntry(gP0Obs,"Observed p_{0}","P");
  */

  // Plotting options:
  gP0Obs->GetYaxis()->SetRangeUser(0.0000001, 1.0);
  //gP0Exp->GetYaxis()->SetRangeUser(0.0000001, 1.0);
  
  gPad->SetLogy();
  //if (m_config->getBool("DoBlind")) gP0Exp->Draw("AP");
  gP0Obs->Draw("AP");
  //leg.Draw("SAME");
  
  // Significance lines and text:
  TLatex sigma; sigma.SetTextColor(kBlack);
  sigma.SetTextFont(42); sigma.SetTextSize(0.05);
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  double sigmaVals[5] = {0.158655, 0.02275, 0.001349, 0.000032, 0.0000002867};
  for (int i_s = 0; i_s < 4; i_s++) {
    line->DrawLine(gP0Obs->GetXaxis()->GetXmin(), sigmaVals[i_s],
		   gP0Obs->GetXaxis()->GetXmax(), sigmaVals[i_s]);
    sigma.DrawLatex(0.9, 1.1*sigmaVals[i_s], Form("%d#sigma",i_s+1));
  }
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.2, 0.27, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.27, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.2, 0.21, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			      (m_config->getNum("AnalysisLuminosity")/1000.0)));
  
  t.DrawLatex(0.50, 0.27, "Spin-2 Selection");
  t.DrawLatex(0.50, 0.21,
	      Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{PI}=%2.2f",
		   ((double)width)/100.0));
  
  // Print the canvas:
  can->Print(Form("%s/p0_width%d.eps", m_outputDir.Data(), width));
  
  // Save p0 graphs to file:
  TFile *outP0File = new TFile(Form("%s/p0_graphs_width%d.root",
				    m_outputDir.Data(), width), "RECREATE");
  gP0Obs->Write();
  can->Write();
  outP0File->Close();
  
  gPad->SetLogy(false);
  
  // Delete pointers:
  printer(Form("CLScan: Finished p0 mass scan for %d!", width), false);
  delete can;
  delete gP0Obs;
  //delete gP0Exp;
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
   @param limitValue - The limit value to set.
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
   Set the p0 which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limits (false for observed).
   @param p0Value - The p0 value to set.
*/
void CLScan::setP0(int mass, int width, bool expected, double p0Value) {
  if (expected) m_valuesP0[Form("exp_mass%d_width%d", mass, width)] = p0Value;
  else m_valuesP0[Form("obs_mass%d_width%d", mass, width)] = p0Value;
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
   The method finds the 95% CL by scanning various signal cross-sections.
   @param configFile - The analysis configuration file.
   @param options - Job options: "New","FromFile","toy","asymptotic","NEvents"
   @param resMass - The resonance mass.
   @param makeNew - True if from toy MC, false if from text file storage.
   @return - True iff loaded successfully.
*/
bool CLScan::singleCLScan(int mass, int width, bool makeNew) {
  printer(Form("singleCLScan(mass=%d, width=%d, makeNew=%d",
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
  
  //----------------------------------------//
  // Open CL values from file:
  if (!makeNew) {
    
    // Open the saved CL values from toys:
    TString inputName = Form("%s/scan_CL_values_mass%d_width%d.txt",
			     m_outputDir.Data(), mass, width);
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
    testStat->setNominalSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
    
    // Set the PoI ranges for this study:
    std::vector<TString> listPoI = m_config->getStrV("WorkspacePoIs");
    for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
      std::vector<double> currRange
	= m_config->getNumV(Form("ScanPoIRange_%s", (listPoI[i_p]).Data()));
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
      testStat->theWorkspace()
	->loadSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
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
    
    // Save values for plotting again!
    std::ofstream outFile(Form("%s/scan_CL_values_mass%d_width%d.txt",
			       m_outputDir.Data(), mass, width));
    
    // Loop over cross-section:
    printer("CLScan: Loop over cross-sections to get qMuObs", false);
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
    printer("CLScan: Processing toy MC in loop over cross-sections", false);
    for (int i_t = 0; i_t < nToyPoints; i_t++) {
      
      // Load the tool to analyze toys.
      // NOTE: this was moved outside the loop above because of interference
      // between the ToyAnalysis class and TestStat, which is called
      // in DHToyAnalysis...
      ToyAnalysis *toyAna = new ToyAnalysis(m_configFileName, "None");
      toyAna->setOutputDir(Form("%s/ToyPlots_mass%d_width%d_xs%d",
				m_outputDir.Data(), mass, width, xsVals[i_t]));
      std::vector<TString> fitTypes; fitTypes.clear();
      fitTypes.push_back("0");
      fitTypes.push_back("1");
      fitTypes.push_back("Free");
      toyAna->setFitTypes(fitTypes);
      toyAna->loadToy(0, Form("%s/toy_mu0*ForScan_mass%d_width%d_xs%d.root",
			      m_inputDir.Data(), mass, width, xsVals[i_t]));
      toyAna->loadToy(1, Form("%s/toy_mu1*ForScan_mass%d_width%d_xs%d.root",
			      m_inputDir.Data(), mass, width, xsVals[i_t]));
      if (!(toyAna->areInputFilesOK())) {
	printer("CLScan: ERROR with toy scan option.", true);
      }
      
      if (m_config->getBool("PlotToysForScan")) {
	
	// Plot the toy MC nuisance parameter, global observables, and PoI:
	std::vector<TString> namesGlobs = toyAna->getNamesGlobalObservables();
	std::vector<TString> namesNuis = toyAna->getNamesNuisanceParameters();
	std::vector<TString> namesPars = toyAna->getNamesPoI();
	for (int i_g = 0; i_g < (int)namesGlobs.size(); i_g++) {
	  if (!(namesGlobs[i_g]).Contains("gamma_stat_channel_bin")) {
	    toyAna->plotHist(namesGlobs[i_g], 0);// Mu=0 toy data
	    toyAna->plotHist(namesGlobs[i_g], 1);// Mu=1 data
	  }
	}
	for (int i_n = 0; i_n < (int)namesNuis.size(); i_n++) {
	  if (!(namesNuis[i_n]).Contains("gamma_stat_channel_bin")) {
	    toyAna->plotHist(namesNuis[i_n], 0);
	    toyAna->plotHist(namesNuis[i_n], 1);
	  }
	}
	for (int i_p = 0; i_p < (int)namesPars.size(); i_p++) {
	  toyAna->plotHist(namesPars[i_p], 0);
	  toyAna->plotHist(namesPars[i_p], 1);
	}
	
	// Then plot statistics:
	toyAna->plotProfiledMu();
	toyAna->plotTestStat("QMu");
	toyAna->plotTestStat("Q0");
	toyAna->plotTestStatComparison("QMu");
	toyAna->plotTestStatComparison("Q0");
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
  gCLExp->GetXaxis()
    ->SetTitle("#sigma_{G*}#timesBR(G*#rightarrow#gamma#gamma [fb]");
  gCLObs->GetXaxis()
    ->SetTitle("#sigma_{G*}#timesBR(G*#rightarrow#gamma#gamma [fb]");
  gCLExp->GetYaxis()->SetTitle("#it{CL_{s}} value");
  gCLObs->GetYaxis()->SetTitle("#it{CL_{s}} value");
  gCLExp->SetLineColor(kBlack);
  gCLObs->SetLineColor(kBlack);
  gCLObs->SetMarkerColor(kBlack);
  gCLObs->SetMarkerStyle(21);  
  gCLExp->SetLineStyle(2);
  gCLObs->SetLineStyle(1);
  gCLExp->SetLineWidth(2);
  gCLObs->SetLineWidth(2);
  gCLExp->GetYaxis()->SetRangeUser(0.0, 1.0);
  gCLExp->GetXaxis()
    ->SetRangeUser(((double)m_xsValues[0])/1000.0,
		   ((double)m_xsValues[(int)m_xsValues.size()-1])/1000.0);
  gPad->SetLogx();
  gCLExp_2s->SetFillColor(kYellow);
  gCLExp_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.64, 0.20, 0.89, 0.38);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!m_config->getBool("DoBlind")) {
    leg.AddEntry(gCLObs, "Observed #it{CL_{s}}", "l");
  }
  leg.AddEntry(gCLExp, "Expected #it{CL_{s}}", "l");
  leg.AddEntry(gCLExp_1s, "Expected #pm 1#sigma_{exp}", "F");
  leg.AddEntry(gCLExp_2s, "Expected #pm 2#sigma_{exp}", "F");
  
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
  gPad->SetLogx(false);
  
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
  setLimit(mass, width, true, 1, expectedCL_p1);
  setLimit(mass, width, true, 2, expectedCL_p2);
  setLimit(mass, width, true, -1, expectedCL_n1);
  setLimit(mass, width, true, -2, expectedCL_n2);
  
  // Then print to screen:
  std::cout << "\nCLScan: Results" << std::endl;
  std::cout << "\tobserved CL = " << observedCL << std::endl;
  std::cout << "\texpected CL = " << expectedCL << std::endl;
  std::cout << "\texpected CL +1 = " << expectedCL_p1 << std::endl;
  std::cout << "\texpected CL +2 = " << expectedCL_p2 << std::endl;
  std::cout << "\texpected CL -1 = " << expectedCL_n1 << std::endl;
  std::cout << "\texpected CL -2 = " << expectedCL_n2 << std::endl;
  
  // Delete pointers, close files, return:
  printer(Form("CLScan::singleCLScan finished mass %d width %d\n", mass, width),
	  false);
  
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
   Calculate the p0 for a given mass and width hypothesis.
   @param configFile - The analysis configuration file.
   @param options - Job options: "New","FromFile","toy","asymptotic","NEvents"
   @param resMass - The resonance mass.
   @param makeNew - True if from toy MC, false if from text file storage.
   @return - True iff loaded successfully.
*/
bool CLScan::singleP0Test(int mass, int width, bool makeNew) {
  printer(Form("singleP0Test(mass=%d, width=%d, makeNew=%d",
	       mass, width, (int)makeNew), false);
  
  //----------------------------------------//
  // Step 1: get observed q0 (no expected!)
  printer("CLScan: Calculating observed q0 for p0 test", false);
  
  // Open the workspace:
  TFile wsFile(m_config->getStr("WorkspaceFile"), "read");
  RooWorkspace *workspace
    = (RooWorkspace*)wsFile.Get(m_config->getStr("WorkspaceName"));
  
  // Instantiate the test statistic class for calculations and plots:
  TestStat *testStat = new TestStat(m_configFileName, "new", workspace);
  testStat->setNominalSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
  
  // Set the PoI ranges for this study:
  std::vector<TString> listPoI = m_config->getStrV("WorkspacePoIs");
  for (int i_p = 0; i_p < (int)listPoI.size(); i_p++) {
    std::vector<double> currRange
      = m_config->getNumV(Form("ScanPoIRange_%s", (listPoI[i_p]).Data()));
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
    testStat->theWorkspace()
      ->loadSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
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
  
  // Also force the mass and width to be constant always:
  testStat->setParam(m_config->getStr("PoIForMass"), (double)mass, true);
  testStat->setParam(m_config->getStr("PoIForWidth"),(((double)width)/100.0), 
		     true);
  
  // map of names and values of pois to set for fit.
  std::map<TString,double> mapPoIMu0; mapPoIMu0.clear();
  mapPoIMu0[m_config->getStr("PoIForNormalization")] = 0.0;//cross-section
  mapPoIMu0[m_config->getStr("PoIForMass")] = (double)mass;
  mapPoIMu0[m_config->getStr("PoIForWidth")] = (((double)width)/100.0);
  
  // Perform the mu=0 fit and mu-free fit (necessary for qmu calculation):
  double nllObsMu0 = testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
					 0, true, mapPoIMu0, false);
  double nllObsMuFree =testStat->getFitNLL(m_config->getStr("WorkspaceObsData"),
					   1, false, mapPoIMu0, false);
  
  // Get profiled signal strength from the mu-free fit:
  std::map<std::string,double> poiFromFit = testStat->getPoIs();
  double obsXSValue
    = poiFromFit[(std::string)m_config->getStr("PoIForNormalization")];
  double q0Observed
    = testStat->getQ0FromNLL(nllObsMu0, nllObsMuFree, obsXSValue);
  
  // Delete pointers and close files:  
  delete testStat;
  delete workspace;
  wsFile.Close();
  
  //----------------------------------------//
  // Process the toy MC files to obtain the CL values:
  printer("CLScan: Computing p0 from toy MC", false);
  
  // Load the tool to analyze toys.
  // NOTE: this was moved  because of intereference between the ToyAnalysis 
  // class and TestStat, which is also called in ToyAnalysis...
  TString toyAnaOptions
    = (m_options.Contains("ImportSamp")) ? "ImportSamp" : "None";
  ToyAnalysis *toyAna = new ToyAnalysis(m_configFileName, toyAnaOptions);
  toyAna->setOutputDir(Form("%s/ToyPlots_mass%d_width%d_xs0",
			    m_outputDir.Data(), mass, width));
  std::vector<TString> fitTypes; fitTypes.clear();
  fitTypes.push_back("0");
  fitTypes.push_back("1");
  fitTypes.push_back("Free");
  toyAna->setFitTypes(fitTypes);
  toyAna->loadToy(0, Form("%s/toy_mu0*ForScan_mass%d_width%d_xs*.root",
			  m_inputDir.Data(), mass, width));
  if (!(toyAna->areInputFilesOK())) {
    printer("CLScan: ERROR with toy scan option.", true);
  }
  
  if (m_config->getBool("PlotToysForScan")) {

    // Plot the toy MC nuisance parameter, global observables, and PoI:
    std::vector<TString> namesGlobs = toyAna->getNamesGlobalObservables();
    std::vector<TString> namesNuis = toyAna->getNamesNuisanceParameters();
    std::vector<TString> namesPars = toyAna->getNamesPoI();
    for (int i_g = 0; i_g < (int)namesGlobs.size(); i_g++) {
      if (!(namesGlobs[i_g]).Contains("gamma_stat_channel_bin")) {
	toyAna->plotHist(namesGlobs[i_g], 0);// Mu=0 toy data
      }
    }
    for (int i_n = 0; i_n < (int)namesNuis.size(); i_n++) {
      if (!(namesNuis[i_n]).Contains("gamma_stat_channel_bin")) {
	toyAna->plotHist(namesNuis[i_n], 0);
      }
    }
    for (int i_p = 0; i_p < (int)namesPars.size(); i_p++) {
      toyAna->plotHist(namesPars[i_p], 0);
    }
    
    // Then plot statistics:
    toyAna->plotTestStat("Q0");
    toyAna->plotTestStatComparison("Q0");
  }
  
  double observedP0 = toyAna->calculateP0FromToy(q0Observed);
  //double expectedP0 = toyAna->calculateP0FromToy(q0Expected);
  delete toyAna;
  
  // Store the p0 for this mass and width:
  setP0(mass, width, false, observedP0);
  //setP0(mass, width, true, expectedP0);
  
  // Then print to screen:
  std::cout << "\nCLScan: P0 Results" << std::endl;
  std::cout << "\tobserved p0 = " << observedP0 << std::endl;
  //std::cout << "\texpected p0 = " << expectedP0 << std::endl;
  
  // Delete pointers, close files, return:
  printer(Form("CLScan::singleP0Test finished mass %d width %d\n", mass, width),
	  false);
  
  // Return true iff. calculation was successful:
  if (toyAna->areInputFilesOK() && observedP0 > 0.0 && observedP0 < 1.0) {
    return true;
  }
  else return false;
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
