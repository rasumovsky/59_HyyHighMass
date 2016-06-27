////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: StatScan.cxx                                                        //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 29/03/2016                                                          //
//                                                                            //
//  Performs a scan of the 95% CL using toy MC at various mass-width-cross-   //
//  section points.                                                           //
//                                                                            //
//  Options:                                                                  //
//  - ImportSamp: implements importance sampling for p0 calculation.          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "StatScan.h"

/**
   -----------------------------------------------------------------------------
   Constructor for the StatScan class.
   @param configFileName - The name of the analysis config file.
   @param options - Options for the CL scan.
*/
StatScan::StatScan(TString configFileName, TString options) {
  
  // Load the config file:
  m_configFileName = configFileName;
  m_config = new Config(m_configFileName);
  m_options = options;

  printer(Form("StatScan::StatScan(%s, %s)", configFileName.Data(), 
	       options.Data()), false);
  
  // Set output directory:
  setInputDirectory(Form("%s/%s/GenericToys/single_files",
			 (m_config->getStr("MasterOutput")).Data(),
			 (m_config->getStr("JobName")).Data()));
  setOutputDirectory(Form("%s/%s/StatScan",
			  (m_config->getStr("MasterOutput")).Data(),
			  (m_config->getStr("JobName")).Data()));
  
  // Clear data:
  clearData();
  
  // Use observed dataset nominally:
  defineDataForFitting(m_config->getStr("WorkspaceObsData"), NULL);
      
  // Set ATLAS style template:
  CommonFunc::SetAtlasStyle();
}

/**
   -----------------------------------------------------------------------------
   Clear the class data.
*/
void StatScan::clearData() {
  m_massValues.clear();
  m_widthValues.clear();
  m_xsValues.clear();
  m_valuesCL.clear();
  m_valuesLimit.clear();
  m_valuesP0.clear();
}

/**
   -----------------------------------------------------------------------------
   @param dataNameForFits - The name of the dataset to use for fits.
   @param dataForFits - A pointer to the datset to use for fits.
*/
void StatScan::defineDataForFitting(TString dataNameForFits,
				    RooAbsData *dataForFits) {
  m_dataNameForFits = dataNameForFits;
  m_dataToFit = dataForFits;
}

/**
   -----------------------------------------------------------------------------
   Look at the available toy MC files to detect which mass, width, and cross-
   section values have associated pseudo-experiment ensembles.
   @param toyDirectory - The directory in which the toy MC files reside.
*/
void StatScan::detectMassWidthXSFiles(TString toyDirectory) {

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
    if ((currToyName.Contains("ForScan") && currToyName.Contains(".root")) ||
	(currToyName.Contains("scan_CL_values") && 
	 currToyName.Contains(".txt"))) {
      
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
	  widthStr.ReplaceAll(".txt", "");
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
double StatScan::getIntercept(TGraph *graph, double valueToIntercept) {
  
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
    std::cout << "StatScan: ERROR! Intercept not found." << std::cout;
    return -999;
  } 
  return currXValue;
}

/**
   -----------------------------------------------------------------------------
   Get the CL which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param crossSection - The crossSection integer for the point of interest.
   @param expected - True for expected limit (false for observed).
   @param asymptotic - True iff asymptotic result is desired.
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @return - The 95% CL value of the cross-section.
*/
double StatScan::getCL(int mass, int width, int crossSection, bool expected,
		       bool asymptotic, int N) {
  TString key = Form("mass%d_width%d_xs%d_exp%d_asym%d_N%d", mass, width,
		     crossSection, (int)expected, (int)asymptotic, N);
  if (m_valuesCL.count(key) == 0) {
    printer(Form("StatScan:getCL ERROR no value for key %s", key.Data()), true);
  }
  return m_valuesCL[key];
}

/**
   -----------------------------------------------------------------------------
   Get the limit which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limit (false for observed).
   @param asymptotic - True iff asymptotic result is desired.
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @return - The 95% CL value of the cross-section.
*/
double StatScan::getLimit(int mass, int width, bool expected, bool asymptotic,
			  int N) {
  TString key = Form("mass%d_width%d_exp%d_asym%d_N%d", mass, width,
		     (int)expected, (int)asymptotic, N);
  if (m_valuesLimit.count(key) == 0) {
    printer(Form("StatScan:getLimit ERROR no value for key %s", key.Data()),
	    true);
  }
  return m_valuesLimit[key];
}

/**
   -----------------------------------------------------------------------------
   Get the p0 which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limit (false for observed).
   @param asymptotic - True iff asymptotic result is desired.
   @return - The p0 value.
*/
double StatScan::getP0(int mass, int width, bool expected, bool asymptotic) {
  TString key = Form("mass%d_width%d_exp%d_asym%d", mass, width,
		     (int)expected, (int)asymptotic);
  if (m_valuesP0.count(key) == 0) {
    printer(Form("StatScan:getP0 ERROR no value for key %s", key.Data()), true);
  }
  return m_valuesP0[key];
}

/**
   -----------------------------------------------------------------------------
   Get a list of the mass values stored in the scan tool.
*/
std::vector<int> StatScan::listMasses() {
  return m_massValues;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the width values stored in the scan tool.
*/
std::vector<int> StatScan::listWidths() {
  return m_widthValues;
}

/**
   -----------------------------------------------------------------------------
   Get a list of the cross-section values stored in the scan tool.
*/
std::vector<int> StatScan::listXS() {
  return m_xsValues;
}

/**
   -----------------------------------------------------------------------------
   Prints a statement (if verbose) and exits (if fatal).
   @param statement - The statement to print.
   @param isFatal - True iff. this should trigger an exit command.
*/
void StatScan::printer(TString statement, bool isFatal) {
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
   @param asymptotic - True iff. asymptotic results are desired. 
   @param doTilde - True iff. using qMuTilde instead of qMu.
*/
void StatScan::scanMassLimit(int width, bool makeNew, bool asymptotic,
			     bool doTilde) {
  printer(Form("StatScan::scanMassLimit(%d, %d, %d, %d)",
	       width, (int)makeNew, (int)asymptotic, (int)doTilde), false);
  
  TString limitFileName = asymptotic ? 
    Form("%s/asymptotic_limit_values_width%d.txt", m_outputDir.Data(), width) :
    Form("%s/toy_limit_values_width%d.txt", m_outputDir.Data(), width);
  
  // Arrays to store limit graph information (medians and error bands):
  int nPlotPoint = 0;
  double observableValues[1000] = {0};  
  double limitObs[1000] = {0};
  double limitExp_p2[1000] = {0};  
  double limitExp_p1[1000] = {0};
  double limitExp[1000] = {0};
  double limitExp_n1[1000] = {0};
  double limitExp_n2[1000] = {0};

  // Load from file here. 
  if (asymptotic && !makeNew) {
    std::ifstream limitFileIn(limitFileName);
    while (limitFileIn >> m_massValues[nPlotPoint] >> limitObs[nPlotPoint] 
	   >> limitExp_n2[nPlotPoint] >> limitExp_n1[nPlotPoint]
	   >> limitExp[nPlotPoint] >> limitExp_p1[nPlotPoint] 
	   >> limitExp_p2[nPlotPoint]) {
      // Store the limits for this mass and width:
      setLimit(m_massValues[nPlotPoint], width, false, asymptotic, 0, 
	       limitObs[nPlotPoint]);
      setLimit(m_massValues[nPlotPoint], width, true, asymptotic, -2, 
	       limitExp_n2[nPlotPoint]);
      setLimit(m_massValues[nPlotPoint], width, true, asymptotic, -1, 
	       limitExp_n1[nPlotPoint]);
      setLimit(m_massValues[nPlotPoint], width, true, asymptotic, 0, 
	       limitExp[nPlotPoint]);
      setLimit(m_massValues[nPlotPoint], width, true, asymptotic, 1, 
	       limitExp_p1[nPlotPoint]);
      setLimit(m_massValues[nPlotPoint], width, true, asymptotic, 2, 
	       limitExp_p2[nPlotPoint]);
      nPlotPoint++;
    }
  }
  
  //----------------------------------------//
  // Loop over the mass points to load or calculate limit results:
  std::ofstream limitFileOut(limitFileName);
  for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
    if ((asymptotic && singleLimitTest(m_massValues[i_m], width, doTilde)) || 
	(!asymptotic && singleCLScan(m_massValues[i_m],width,makeNew,false))) {
      
      observableValues[nPlotPoint] = m_massValues[i_m];
         
      limitObs[nPlotPoint]
	= getLimit(m_massValues[i_m], width, false, asymptotic, 0);
      limitExp_n2[nPlotPoint] 
	= getLimit(m_massValues[i_m], width, true, asymptotic, -2);
      limitExp_n1[nPlotPoint]
	= getLimit(m_massValues[i_m], width, true, asymptotic, -1);
      limitExp[nPlotPoint]
	= getLimit(m_massValues[i_m], width, true, asymptotic, 0);
      limitExp_p1[nPlotPoint]
	= getLimit(m_massValues[i_m], width, true, asymptotic, 1);
      limitExp_p2[nPlotPoint]
	= getLimit(m_massValues[i_m], width, true, asymptotic, 2);

      limitFileOut << m_massValues[i_m] << " " 
		   << limitObs[nPlotPoint] << " " 
		   << limitExp_n2[nPlotPoint] << " " 
		   << limitExp_n1[nPlotPoint] << " "
		   << limitExp[nPlotPoint] << " " 
		   << limitExp_p1[nPlotPoint] << " " 
		   << limitExp_p2[nPlotPoint] << " " << std::endl;
      nPlotPoint++;
    }
  }
  limitFileOut.close();

  if (nPlotPoint >= 1000) printer("StatScan: Array bound exceeded",true);

  //----------------------------------------//
  // Plot the results:
  double errExp_p2[1000] = {0};  
  double errExp_p1[1000] = {0};
  double errExp_n1[1000] = {0};
  double errExp_n2[1000] = {0};
  for (int i_t = 0; i_t < nPlotPoint; i_t++) {
    errExp_p2[i_t] = fabs(limitExp_p2[i_t] - limitExp[i_t]);
    errExp_p1[i_t] = fabs(limitExp_p1[i_t] - limitExp[i_t]);
    errExp_n1[i_t] = fabs(limitExp_n1[i_t] - limitExp[i_t]);
    errExp_n2[i_t] = fabs(limitExp_n2[i_t] - limitExp[i_t]);
  }
  
  // Graphs for the median expected and observed results:
  TGraph *gLimitExp = new TGraph(nPlotPoint, observableValues, limitExp);
  TGraph *gLimitObs = new TGraph(nPlotPoint, observableValues, limitObs);
  
  // Graphs for the expected error bands:
  TGraphAsymmErrors *gLimitExp_2s 
    = new TGraphAsymmErrors(nPlotPoint, observableValues, limitExp, 0, 0, 
 			    errExp_n2, errExp_p2);
  TGraphAsymmErrors *gLimitExp_1s
    = new TGraphAsymmErrors(nPlotPoint, observableValues, limitExp, 0, 0, 
			    errExp_n1, errExp_p1);
  
  // Canvas and graph formatting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  if ((m_config->getStr("AnalysisType")).Contains("Scalar")) {
    gLimitExp->GetXaxis()->SetTitle("m_{X} [GeV]");
    gLimitObs->GetXaxis()->SetTitle("m_{X} [GeV]");
    gLimitExp_2s->GetXaxis()->SetTitle("m_{X} [GeV]");
    gLimitExp->GetYaxis()->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrow#gamma#gamma} [fb]");
    gLimitObs->GetYaxis()->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrow#gamma#gamma} [fb]");
    gLimitExp_2s->GetYaxis()->SetTitle("95% CL limit on #sigma_{X}#timesBR_{X#rightarrow#gamma#gamma} [fb]");
  }
  else {
    gLimitExp->GetXaxis()->SetTitle("m_{G*} [GeV]");
    gLimitObs->GetXaxis()->SetTitle("m_{G*} [GeV]");
    gLimitExp_2s->GetXaxis()->SetTitle("m_{G*} [GeV]");
    gLimitExp->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [fb]");
    gLimitObs->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [fb]");
    gLimitExp_2s->GetYaxis()->SetTitle("95% CL limit on #sigma_{G*}#timesBR_{G*#rightarrow#gamma#gamma} [fb]");
  }
  
  gLimitExp->SetLineColor(kBlack);
  gLimitExp->SetLineStyle(2);
  gLimitExp->SetLineWidth(2);
  
  gLimitObs->SetLineColor(kBlack);
  gLimitObs->SetLineStyle(1);
  gLimitObs->SetLineWidth(2);
  gLimitObs->SetMarkerColor(kBlack);
  gLimitObs->SetMarkerStyle(21);
  
  gLimitExp_2s->SetFillColor(kYellow);
  gLimitExp_1s->SetFillColor(kGreen);
  
  // Legend:
  TLegend leg(0.61, 0.68, 0.89, 0.91);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextSize(0.04);
  if (!m_config->getBool("DoBlind")) {
    leg.AddEntry(gLimitObs,"Observed #it{CL_{s}} limit","P");
  }
  leg.AddEntry(gLimitExp,"Expected #it{CL_{s}} limit","l");
  leg.AddEntry(gLimitExp_1s,"Expected #pm 1#sigma_{exp}","F");
  leg.AddEntry(gLimitExp_2s,"Expected #pm 2#sigma_{exp}","F");
  
  // Plotting options:
  gLimitExp_2s->GetXaxis()->SetRangeUser(m_massValues[0], 
					 m_massValues[m_massValues.size()-1]);
  gLimitExp_2s->GetYaxis()->SetRangeUser(1, 1000);
  
  gPad->SetLogy();
  gLimitExp_2s->Draw("A3");
  gLimitExp->Draw("Lsame");
  gLimitExp_1s->Draw("3same");
  gLimitExp->Draw("LSAME");
  if (!m_config->getBool("DoBlind")) gLimitObs->Draw("PSAME");
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
  if ((m_config->getStr("AnalysisType")).Contains("Scalar")) {
    t.DrawLatex(0.2, 0.75, "Spin-0 Selection");
    t.DrawLatex(0.2, 0.69,
		Form("G*#rightarrow#gamma#gamma, #Gamma/m_{X}=%2.2f",
		     ((double)width)/1000.0));
  }
  else {
    t.DrawLatex(0.2, 0.75, "Spin-2 Selection");
    t.DrawLatex(0.2, 0.69,
		Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{PI}=%2.2f",
		     ((double)width)/1000.0));
  }
  
  // Print the canvas:
  can->Print(Form("%s/limits_width%d.eps", m_outputDir.Data(), width));
  gPad->SetLogy(false);
  
  // Save graphs to file:
  TString graphFileName = asymptotic ? 
    Form("%s/asymptotic_limit_graphs_width%d.root", m_outputDir.Data(), width) :
    Form("%s/toy_limit_graphs_width%d.root", m_outputDir.Data(), width);
  TFile *outLimitFile = new TFile(graphFileName, "RECREATE");
  gLimitObs->Write();
  gLimitExp->Write();
  gLimitExp_2s->Write();
  gLimitExp_1s->Write();
  can->Write();
  outLimitFile->Close();
  
  // Delete pointers:
  printer(Form("StatScan: Finished mass scan for %d!", width), false);
  delete can;
  delete gLimitObs;
  delete gLimitExp;
  delete gLimitExp_2s;
  delete gLimitExp_1s;
}

/**
   -----------------------------------------------------------------------------
   Plot the p0 value as a function of mass, for a given width slice.
   @param width - The width integer for the scan.
   @param makeNew - True iff. calculating things from toys directly.
   @param asymptotic - True iff. asymptotic results are desired. 
*/
void StatScan::scanMassP0(int width, bool makeNew, bool asymptotic) {
  printer(Form("StatScan::scanMassP0(%d, %d, %d)",
	       width, (int)makeNew, (int)asymptotic), false);
  
  // Arrays to store p0 graph information:
  int nPlotPoint = 0;
  double observableValues[1000] = {0};  
  double p0Obs[1000] = {0};
  double p0Exp[1000] = {0};
  
  //----------------------------------------//
  // Loop over the mass points to load or calculate p0:
  printer("StatScan: Retrieving p0...", false);
  for (int i_m = 0; i_m < (int)m_massValues.size(); i_m++) {
    double currentXS = 1.0;
    if (singleP0Test(m_massValues[i_m], width, currentXS, makeNew, asymptotic)){
      observableValues[nPlotPoint] = m_massValues[i_m];
      p0Obs[nPlotPoint] = getP0(m_massValues[i_m], width, false, asymptotic);
      p0Exp[nPlotPoint] = getP0(m_massValues[i_m], width, true, asymptotic);
      nPlotPoint++;
    }
  }
  if (nPlotPoint >= 1000) printer("StatScan: Array bound exceeded",true);

  //----------------------------------------//
  // Plot the results:
  
  // Median expected and observed results:
  TGraph *gP0Exp = new TGraph(nPlotPoint, observableValues, p0Exp);
  TGraph *gP0Obs = new TGraph(nPlotPoint, observableValues, p0Obs);
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  if ((m_config->getStr("AnalysisType")).Contains("Scalar")) {
    gP0Exp->GetXaxis()->SetTitle("m_{X} [GeV]");
    gP0Obs->GetXaxis()->SetTitle("m_{X} [GeV]");
  }
  else {
    gP0Exp->GetXaxis()->SetTitle("m_{G*} [GeV]");
    gP0Obs->GetXaxis()->SetTitle("m_{G*} [GeV]");
  }
  gP0Exp->GetYaxis()->SetTitle("p_{0}");
  gP0Obs->GetYaxis()->SetTitle("p_{0}");  
  gP0Exp->SetLineColor(kBlack);
  gP0Obs->SetLineColor(kBlack);
  gP0Exp->SetLineStyle(2);
  gP0Obs->SetLineStyle(1);
  gP0Exp->SetLineWidth(2);
  gP0Obs->SetLineWidth(2);
  
  // Legend:
  TLegend leg(0.65, 0.27, 0.89, 0.33);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.05);
  leg.AddEntry(gP0Obs,"Observed p_{0}","L");
  //leg.AddEntry(gP0Exp,"Expected p_{0}","L");
  
  // Plotting options:
  gP0Obs->GetYaxis()->SetRangeUser(0.0000001, 1.0);
    
  gPad->SetLogy();
  gP0Obs->Draw("AL");
  //  gP0Exp->Draw("LPSAME");
  leg.Draw("SAME");
  
  // Significance lines and text:
  TLatex sigma; sigma.SetTextColor(kRed+1);
  sigma.SetTextFont(42); sigma.SetTextSize(0.04);
  TLine *line = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(1);
  line->SetLineColor(kRed+1);
  double sigmaVals[6] = {0.5,0.15865,0.02275,0.001349,0.000032,0.0000002867};
  for (int i_s = 0; i_s < 5; i_s++) {
    double sigmaXPos = gP0Obs->GetXaxis()->GetXmax()
      - (0.07*(gP0Obs->GetXaxis()->GetXmax() - gP0Obs->GetXaxis()->GetXmin()));
    line->DrawLine(gP0Obs->GetXaxis()->GetXmin(), sigmaVals[i_s],
		   gP0Obs->GetXaxis()->GetXmax(), sigmaVals[i_s]);
    sigma.DrawLatex(sigmaXPos, 1.1*sigmaVals[i_s], Form("%d#sigma",i_s));
  }
  
  // Print ATLAS text on the plot:    
  TLatex t; t.SetNDC(); t.SetTextColor(kBlack);
  t.SetTextFont(72); t.SetTextSize(0.05);
  t.DrawLatex(0.2, 0.33, "ATLAS");
  t.SetTextFont(42); t.SetTextSize(0.05);
  t.DrawLatex(0.32, 0.33, m_config->getStr("ATLASLabel"));
  t.DrawLatex(0.2, 0.27, Form("#sqrt{s} = 13 TeV, %2.1f fb^{-1}",
			      (m_config->getNum("AnalysisLuminosity")/1000.0)));
  if ((m_config->getStr("AnalysisType")).Contains("Scalar")) {
    t.DrawLatex(0.2, 0.21, "Spin-0 Selection");
    if (width == 0) {
      t.DrawLatex(0.65, 0.21, "X#rightarrow#gamma#gamma, NWA");
    }
    else {
      t.DrawLatex(0.65, 0.21,
		  Form("X#rightarrow#gamma#gamma, #Gamma/m_{X}=%2.2f",
		       ((double)width)/1000.0));
    }
  }
  else if ((m_config->getStr("AnalysisType")).Contains("GravitonLoose")) {
    t.DrawLatex(0.2, 0.21, "Spin-2 Loose Iso.");
    t.DrawLatex(0.65, 0.21,
		Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{PI}=%2.2f",
		     ((double)width)/1000.0));
  }
  else {
    t.DrawLatex(0.2, 0.21, "Spin-2 Selection");
    t.DrawLatex(0.65, 0.21,
		Form("G*#rightarrow#gamma#gamma, #it{k}/#bar{M}_{PI}=%2.2f",
		     ((double)width)/1000.0));
  }

  gP0Obs->Draw("LSAME");

  // Print the canvas:
  can->Print(Form("%s/p0_width%d.eps", m_outputDir.Data(), width));
  
  // Save p0 graphs to file:
  TString graphFileName = asymptotic ? 
    Form("%s/asymptotic_p0_graphs_width%d.root", m_outputDir.Data(), width) :
    Form("%s/toy_p0_graphs_width%d.root", m_outputDir.Data(), width);
  TFile *outP0File = new TFile(graphFileName, "RECREATE");
  gP0Obs->Write();
  can->Write();
  outP0File->Close();
  
  gPad->SetLogy(false);
  
  // Delete pointers:
  printer(Form("StatScan: Finished p0 mass scan for %d!", width), false);
  delete can;
  delete gP0Obs;
  delete gP0Exp;
}

/**
   -----------------------------------------------------------------------------
   Set the input file location.
   @param directory - The directory containing input files.
*/
void StatScan::setInputDirectory(TString directory) {
  m_inputDir = directory;
  // Also detect the input files in this new directory:
  detectMassWidthXSFiles(m_inputDir);
}

/**
   -----------------------------------------------------------------------------
   Set the CL which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param crossSection - The cross-section integer for the point of interest.
   @param expected - True for expected limits (false for observed).
   @param asymptotic - True iff. asymptotic value being set.
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @param CLValue - The CL value to set.
*/
void StatScan::setCL(int mass, int width, int crossSection, bool expected,
		     bool asymptotic, int N, double CLValue) {
  TString key = Form("mass%d_width%d_xs%d_exp%d_asym%d_N%d", mass, width,
		     crossSection, (int)expected, (int)asymptotic, N);
  m_valuesCL[key] = CLValue;
  printer(Form("StatScan::setCL(%s)=%f", key.Data(), m_valuesCL[key]), false);
}

/**
   -----------------------------------------------------------------------------
   Set the limit which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limits (false for observed).
   @param asymptotic - True iff. asymptotic value being set.
   @param N - For bands, use +2,+1,-1,-2. For medians, use 0.
   @param limitValue - The limit value to set.
*/
void StatScan::setLimit(int mass, int width, bool expected, bool asymptotic, 
			int N, double limitValue) {
  TString key = Form("mass%d_width%d_exp%d_asym%d_N%d", mass, width,
		     (int)expected, (int)asymptotic, N);
  m_valuesLimit[key] = limitValue;
  printer(Form("StatScan::setLimit(%s)=%f", key.Data(), m_valuesLimit[key]), 
	  false);
}

/**
   -----------------------------------------------------------------------------
   Set the p0 which was computed for a given mass, width, and sigma.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param expected - True for expected limits (false for observed).
   @param asymptotic - True iff. asymptotic value being set.
   @param p0Value - The p0 value to set.
*/
void StatScan::setP0(int mass, int width, bool expected, bool asymptotic, 
		     double p0Value) {
  TString key = Form("mass%d_width%d_exp%d_asym%d", mass, width,
		     (int)expected, (int)asymptotic);
  m_valuesP0[key] = p0Value;
  //printer(Form("StatScan::setP0(%s)=%f", key.Data(), m_valuesP0[key]), false);
}

/**
   -----------------------------------------------------------------------------
   Set the output file location.
   @param directory - The directory for storing output files.
*/
void StatScan::setOutputDirectory(TString directory) {
  m_outputDir = directory;
  // Create output directory if it doesn't already exist:
  system(Form("mkdir -vp %s", m_outputDir.Data()));
}

/**
   -----------------------------------------------------------------------------
   The method finds the 95% CL for a signal by scanning over cross-sections.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param makeNew - True if from toy MC, false if from text file storage.
   @param asymptotic - True iff. asymptotic results are desired. 
   @return - True iff loaded successfully.
*/
bool StatScan::singleCLScan(int mass, int width, bool makeNew, bool asymptotic){
  printer(Form("StatScan::singleCLScan(mass=%d, width=%d, makeNew=%d",
	       mass, width, (int)makeNew), false);
  
  // Arrays to store band information:
  int nPlotPoint = 0;
  double varValues[1000] = {0};  
  double CLObs[1000] = {0};
  double CLExp_p2[1000] = {0};  
  double CLExp_p1[1000] = {0};
  double CLExp[1000] = {0};
  double CLExp_n1[1000] = {0};
  double CLExp_n2[1000] = {0};
  
  //----------------------------------------//
  // Scan CL for various cross-sections:
  printer("StatScan::singleCLScan: Loop over cross-sections to get CL", false);
  for (int i_x = 0; i_x < (int)m_xsValues.size(); i_x++) {
    if (singleCLTest(mass, width, m_xsValues[i_x], makeNew, asymptotic)) {
      varValues[nPlotPoint] = ((double)m_xsValues[i_x])/1000.0;
      CLObs[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], false, asymptotic, 0);
      CLExp_n2[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], true, asymptotic, -2);
      CLExp_n1[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], true, asymptotic, -1);
      CLExp[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], true, asymptotic, 0);
      CLExp_p1[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], true, asymptotic, 1);
      CLExp_p2[nPlotPoint]
	= getCL(mass, width, m_xsValues[i_x], true, asymptotic, 2);
      nPlotPoint++;
    }
    else {
      printer(Form("StatScan::singleCLScan: ERROR for cross-section %d", 
		   m_xsValues[i_x]), false);
    }
  }
  
  // Check that array is large enough:
  if (nPlotPoint >= 1000) printer("StatScan: Array bound exceeded",true);
  
  //----------------------------------------//
  // Plot the CL scan results:
  double errExp_p2[1000] = {0};  
  double errExp_p1[1000] = {0};
  double errExp_n1[1000] = {0};
  double errExp_n2[1000] = {0};
  for (int i_t = 0; i_t < nPlotPoint; i_t++) {
    errExp_p2[i_t] = fabs(CLExp_p2[i_t] - CLExp[i_t]);
    errExp_p1[i_t] = fabs(CLExp_p1[i_t] - CLExp[i_t]);
    errExp_n1[i_t] = fabs(CLExp_n1[i_t] - CLExp[i_t]);
    errExp_n2[i_t] = fabs(CLExp_n2[i_t] - CLExp[i_t]);
  }
  
  // Median expected and observed results:
  TGraph *gCLObs = new TGraph(nPlotPoint, varValues, CLObs);
  TGraph *gCLExp = new TGraph(nPlotPoint, varValues, CLExp);
  TGraph *gCLExp_p1 = new TGraph(nPlotPoint, varValues, CLExp_p1);
  TGraph *gCLExp_p2 = new TGraph(nPlotPoint, varValues, CLExp_p2);
  TGraph *gCLExp_n1 = new TGraph(nPlotPoint, varValues, CLExp_n1);
  TGraph *gCLExp_n2 = new TGraph(nPlotPoint, varValues, CLExp_n2);
  
  // Also plot the bands:
  TGraphAsymmErrors *gCLExp_2s 
    = new TGraphAsymmErrors(nPlotPoint, varValues, CLExp, 0, 0, 
			    errExp_n2, errExp_p2);
  TGraphAsymmErrors *gCLExp_1s
    = new TGraphAsymmErrors(nPlotPoint, varValues, CLExp, 0, 0, 
			    errExp_n1, errExp_p1);
  
  // Start plotting:
  TCanvas *can = new TCanvas("can","can");
  can->cd();
  
  // Toy graph formatting:
  if ((m_config->getStr("AnalysisType")).Contains("Scalar")) {
    gCLExp->GetXaxis()
      ->SetTitle("#sigma_{X}#timesBR(X#rightarrow#gamma#gamma [fb]");
    gCLObs->GetXaxis()
      ->SetTitle("#sigma_{X}#timesBR(X#rightarrow#gamma#gamma [fb]");
  }
  else {
    gCLExp->GetXaxis()
      ->SetTitle("#sigma_{G*}#timesBR(G*#rightarrow#gamma#gamma [fb]");
    gCLObs->GetXaxis()
      ->SetTitle("#sigma_{G*}#timesBR(G*#rightarrow#gamma#gamma [fb]");
  }
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
  double observedLimit = getIntercept(gCLObs, 0.95);
  double expectedLimit_n2 = getIntercept(gCLExp_n2, 0.95);  
  double expectedLimit_n1 = getIntercept(gCLExp_n1, 0.95);
  double expectedLimit = getIntercept(gCLExp, 0.95);
  double expectedLimit_p1 = getIntercept(gCLExp_p1, 0.95);
  double expectedLimit_p2 = getIntercept(gCLExp_p2, 0.95);
  
  // Store the limits for this mass and width (reverse bands for this):
  setLimit(mass, width, false, asymptotic, 0, observedLimit);
  setLimit(mass, width, true, asymptotic, 2, expectedLimit_n2);
  setLimit(mass, width, true, asymptotic, 1, expectedLimit_n1);
  setLimit(mass, width, true, asymptotic, 0, expectedLimit);
  setLimit(mass, width, true, asymptotic, -1, expectedLimit_p1);
  setLimit(mass, width, true, asymptotic, -2, expectedLimit_p2);
  
  // Then print to screen:
  std::cout << "\nCLScan: Results" << std::endl;
  std::cout << "\tobserved CL = " << observedLimit << std::endl;
  std::cout << "\texpected CL = " << expectedLimit << std::endl;
  std::cout << "\texpected CL +1 = " << expectedLimit_p1 << std::endl;
  std::cout << "\texpected CL +2 = " << expectedLimit_p2 << std::endl;
  std::cout << "\texpected CL -1 = " << expectedLimit_n1 << std::endl;
  std::cout << "\texpected CL -2 = " << expectedLimit_n2 << std::endl;
  
  // Delete pointers, close files, return:
  printer(Form("StatScan::singleCLScan finished mass %d width %d\n", 
	       mass, width), false);
  
  delete line;
  delete can;
  delete gCLObs;
  delete gCLExp;
  delete gCLExp_2s;
  delete gCLExp_1s;
  
  // Return true iff. calculation was successful:
  if (nPlotPoint < 2 || observedLimit < 0 || expectedLimit < 0) return false;
  else return true;
}

/**
   -----------------------------------------------------------------------------
   The method finds the CL for a given signal cross-section.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param crossSection - The cross-section integer for the point of interest.
   @param makeNew - True if from toy MC, false if from text file storage.
   @param asymptotic - True iff. asymptotic results are desired. 
   @param doTilde - True iff. using qMuTilde instead of qMu.
   @return - True iff loaded successfully.
*/
bool StatScan::singleCLTest(int mass, int width, int crossSection,
			    bool makeNew, bool asymptotic, bool doTilde) {
  printer(Form("StatScan::singleCLTest(mass=%d, width=%d, makeNew=%d)",
	       mass, width, (int)makeNew), false);

  bool successful = true;
  
  if (!asymptotic) {
    printer("singleCLTest: ERROR! Toys need to be updated for CL with qmutilde",
	    true);
  }

  // Store band information:
  double CLObs = 0.0;
  double CLExp_p2 = 0.0;
  double CLExp_p1 = 0.0;
  double CLExp = 0.0;
  double CLExp_n1 = 0.0;
  double CLExp_n2 = 0.0;
  
  TString textFileNameCL = asymptotic ? 
    Form("%s/asymptotic_CL_values_mass%d_width%d_xs%d.txt", m_outputDir.Data(),
	 mass, width, crossSection) :
    Form("%s/toy_CL_values_mass%d_width%d_xs%d.txt", m_outputDir.Data(),
	 mass, width, crossSection);
  
  //----------------------------------------//
  // Open CL values from file:
  if (!makeNew) {
    printer(Form("StatScan::singleCLTest: Loading CL from file %s",
		 textFileNameCL.Data()), false);
    
    // Open the saved CL values from toys:
    std::ifstream inputFile(textFileNameCL);
    if (inputFile.is_open()) {
      while (inputFile >> CLObs >> CLExp_n2 >> CLExp_n1 >> CLExp >> CLExp_p1
	     >> CLExp_p2) {
	std::cout << "observed CL = " << CLObs << std::endl;
	std::cout << "expected CL -2sigma = " << CLExp_n2 << std::endl;
	std::cout << "expected CL -1sigma = " << CLExp_n1 << std::endl;
	std::cout << "expected CL = " << CLExp << std::endl;
	std::cout << "expected CL +1sigma = " << CLExp_p1 << std::endl;
	std::cout << "expected CL +2sigma = " << CLExp_p2 << std::endl;
      }
    }
    else {
      printer(Form("StatScan: ERROR opnening text file %s",
		   textFileNameCL.Data()), false);
      successful = false;
    }
    inputFile.close();
  }
  
  //----------------------------------------//
  // Calculate new CL values:
  else {
    printer("StatScan::singleCLTest: Calculating CL from scratch", false);
    
    // Save values for plotting again!
    std::ofstream outFile(textFileNameCL);

    // Set the dataset to fit:
    //TString datasetToFit = m_config->getStr("WorkspaceObsData");
    
    double xSectionDouble = ((double)crossSection)/1000.0;
    
    // Load the workspace:
    // First check for local copy, then go to central copy.
    TString workspaceFileName = m_config->getStr("WorkspaceFile");
    TObjArray *array = workspaceFileName.Tokenize("/");
    
    workspaceFileName
      = ((TObjString*)array->At(array->GetEntries()-1))->GetString();
    TFile *wsFile = new TFile(workspaceFileName);
    if (!wsFile->IsZombie()) {
      printer(Form("StatScan: Loaded ws from %s", workspaceFileName.Data()),
	      false);
    }
    else {
      workspaceFileName = m_config->getStr("WorkspaceFile");
      printer(Form("StatScan: Load ws from %s", workspaceFileName.Data()),
	      false);
      wsFile = new TFile(workspaceFileName, "read");
    }
    RooWorkspace *workspace
      = (RooWorkspace*)wsFile->Get(m_config->getStr("WorkspaceName"));
    
    // Instantiate the test statistic class for calculations and plots:
    TestStat *testStat = new TestStat(m_configFileName, "new", workspace);
    testStat->setNominalSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
    if (!testStat->theWorkspace()->data(m_dataNameForFits) && m_dataToFit) {
      testStat->theWorkspace()->import(*m_dataToFit);
    }
    
    // Strategy settings (if defined):
    if (m_config->isDefined("FitOptions")) {
      testStat->setFitOptions(m_config->getStr("FitOptions"));
    }
    
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
	printer(Form("StatScan: Workspace has no variable %s",
		     listPoI[i_p].Data()), true);
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
    
    // Map of names and values of PoIs to set for fit:
    std::map<TString,double> mapPoI; mapPoI.clear();
    mapPoI[m_config->getStr("PoIForNormalization")] = xSectionDouble;
    mapPoI[m_config->getStr("PoIForMass")] = (double)mass;
    mapPoI[m_config->getStr("PoIForWidth")] = (((double)width)/1000.0);
    
    
    //----------------------------------------//
    // Get the asymptotic CL results:
    if (asymptotic) {
      std::vector<double> asymptoticCLValues
	= testStat->asymptoticCL(mapPoI, m_dataNameForFits,
				 m_config->getStr("WorkspaceSnapshotMu1"),
				 m_config->getStr("PoIForNormalization"),
				 doTilde);
      CLObs = asymptoticCLValues[0];
      CLExp_n2 = asymptoticCLValues[1];
      CLExp_n1 = asymptoticCLValues[2];
      CLExp = asymptoticCLValues[3];
      CLExp_p1 = asymptoticCLValues[4];
      CLExp_p2 = asymptoticCLValues[5];
      
      successful = testStat->fitsAllConverged();
      
      delete testStat;
      delete workspace;
    }
    
    //----------------------------------------//
    // Get CL results using toy MC:
    else {
      printer("StatScan: Calculating observed qMu...", false);
      
      // Also force the mass and width to be constant always:
      testStat->setParam(m_config->getStr("PoIForMass"), 
			 (double)mass, true);
      testStat->setParam(m_config->getStr("PoIForWidth"), 
			 (((double)width)/1000.0), true);
      
      // Perform the mu=1 fit and mu-free fit (necessary for qmu calculation):
      double nllObsMu1 
	= testStat->getFitNLL(m_dataNameForFits, 1, true, mapPoI, false);
      double nllObsMuFree
	= testStat->getFitNLL(m_dataNameForFits, 1, false, mapPoI, false);
      
      // Get profiled signal strength from the mu-free fit:
      std::map<std::string,double> poiFromFit = testStat->getPoIs();
      double obsXSValue 
	= poiFromFit[(std::string)m_config->getStr("PoIForNormalization")];
      double muHat = obsXSValue / xSectionDouble;
      double muForQMu = 1.0;
      double qMuObs = testStat->getQMuFromNLL(nllObsMu1, nllObsMuFree,
					      muHat, muForQMu);
      // If using qMuTilde:
      if (doTilde) {
	double originNorm = mapPoI[m_config->getStr("PoIForNormalization")];
	mapPoI[m_config->getStr("PoIForNormalization")] = 0.0;
	double nllObsMu0
	  = testStat->getFitNLL(m_dataNameForFits, 0.0, true, mapPoI, false);
	mapPoI[m_config->getStr("PoIForNormalization")] = originNorm;
	qMuObs = testStat->getQMuTildeFromNLL(nllObsMu1, nllObsMu0, 
					      nllObsMuFree, muHat, muForQMu);
      }
      
      delete testStat;
      delete workspace;
      wsFile->Close();
      
      //----------------------------------------//
      // Process the toy MC files to obtain the CL values:
      printer("StatScan: Processing toy MC for qMu->CL", false);
      
      // Load the tool to analyze toys.
      // NOTE: interference between the ToyAnalysis class and TestStat
      ToyAnalysis *toyAna = new ToyAnalysis(m_configFileName, "None");
      toyAna->setOutputDir(Form("%s/ToyPlots_mass%d_width%d_xs%d",
				m_outputDir.Data(), mass, width, crossSection));
      std::vector<TString> fitTypes; fitTypes.clear();
      fitTypes.push_back("0");
      fitTypes.push_back("1");
      fitTypes.push_back("Free");
      toyAna->setFitTypes(fitTypes);
      toyAna->loadToy(0, Form("%s/toy_mu0*ForScan_mass%d_width%d_xs%d.root",
			      m_inputDir.Data(), mass, width, crossSection));
      toyAna->loadToy(1, Form("%s/toy_mu1*ForScan_mass%d_width%d_xs%d.root",
			      m_inputDir.Data(), mass, width, crossSection));
      if (toyAna->areInputFilesOK()) {
	
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
	  toyAna->plotTestStat("QMuTilde");
	  toyAna->plotTestStat("QMu");
	  toyAna->plotTestStat("Q0");
	  toyAna->plotTestStatComparison("QMuTilde");
	  toyAna->plotTestStatComparison("QMu");
	  toyAna->plotTestStatComparison("Q0");
	}
	
	// Calculate the expected qMu:
	double qMuExp_n2 = toyAna->calculateBkgQMuForN(-2.0);
	double qMuExp_n1 = toyAna->calculateBkgQMuForN(-1.0);
	double qMuExp = toyAna->calculateBkgQMuForN(0);      
	double qMuExp_p1 = toyAna->calculateBkgQMuForN(1.0);
	double qMuExp_p2 = toyAna->calculateBkgQMuForN(2.0);
	
	// Finally, calculate observed and expected CL:
	CLObs = toyAna->calculateCLFromToy(qMuObs);
	CLExp_n2 = toyAna->calculateCLFromToy(qMuExp_n2);
	CLExp_n1 = toyAna->calculateCLFromToy(qMuExp_n1);
	CLExp = toyAna->calculateCLFromToy(qMuExp);
	CLExp_p1 = toyAna->calculateCLFromToy(qMuExp_p1);
	CLExp_p2 = toyAna->calculateCLFromToy(qMuExp_p2);
	
	// Write CL values to file:
	outFile << CLObs << " " << CLExp_n2 << " " << CLExp_n1 << " " 
		<< CLExp << " " << CLExp_p1 << " " << CLExp_p2 << std::endl;
      }
      else printer("StatScan: ERROR with toy scan option.", false);
      
      successful = toyAna->areInputFilesOK();
      delete toyAna;
    }
    
    // Close the files that save CL data:
    outFile.close();
  }
  
  if (successful) {
    setCL(mass, width, crossSection, false, asymptotic, 0, CLObs);
    setCL(mass, width, crossSection, true, asymptotic, -2, CLExp_n2);
    setCL(mass, width, crossSection, true, asymptotic, -1, CLExp_n1);
    setCL(mass, width, crossSection, true, asymptotic, 0, CLExp);
    setCL(mass, width, crossSection, true, asymptotic, 1, CLExp_p1);
    setCL(mass, width, crossSection, true, asymptotic, 2, CLExp_p2);
    
    // Then print to screen:
    std::cout << "\nStatScan: CL Results:" << std::endl;
    std::cout << "\tobserved CL = " << CLObs << std::endl;
    std::cout << "\texpected CL -2sigma = " << CLExp_n2 << std::endl;
    std::cout << "\texpected CL -1sigma = " << CLExp_n1 << std::endl;
    std::cout << "\texpected CL = " << CLExp << std::endl;
    std::cout << "\texpected CL +1sigma = " << CLExp_p1 << std::endl;
    std::cout << "\texpected CL +2sigma = " << CLExp_p2 << std::endl;
  }
  return successful;
}

/**
   -----------------------------------------------------------------------------
   This method uses Aaron Armbruster's macro to find the 95%CLs exclusion
   with asymptotics.
   @param mass - The mass integer for the point of interest.
   @param width - The width integer for the point of interest.
   @param doTilde - True iff. using qMuTilde instead of qMu.
   @return - True iff calculated successfully.
*/
bool StatScan::singleLimitTest(int mass, int width, bool doTilde) {
  printer(Form("StatScan::singleLimitTest(mass=%d, width=%d)", mass, width),
	  false);
  
  // Options for AsymptoticsCLs.cxx macro:
  TString asymptoticsOptions = "SetVal_NoTemplateStat";
  if (doTilde) asymptoticsOptions.Append("_DoTilde");
  
  system(Form("./bin/AsymptoticsCLs %s %s %s=%d %s=%f", 
	      m_configFileName.Data(), asymptoticsOptions.Data(),
	      (m_config->getStr("PoIForMass")).Data(), mass,
	      (m_config->getStr("PoIForWidth")).Data(), 
	      ((double)width/1000.0)));
  
  TString asymptoticsCLsDir = Form("%s/%s/AsymptoticsCls", 
				   (m_config->getStr("MasterOutput")).Data(),
				   (m_config->getStr("JobName")).Data());
  std::ifstream limitInput(Form("%s/text_limits__%s_%2.2f_%s_%2.2f.txt", 
				asymptoticsCLsDir.Data(),
				(m_config->getStr("PoIForMass")).Data(), 
				((double)mass),
				(m_config->getStr("PoIForWidth")).Data(), 
				((double)width/1000.0)));
  std::map<TString,double> asymptoticLimits; asymptoticLimits.clear();
  if (limitInput.is_open()) {
    while (limitInput >> asymptoticLimits["Obs"]
	   >> asymptoticLimits["ExpN0"] >> asymptoticLimits["ExpN2"]
	   >> asymptoticLimits["ExpN1"] >> asymptoticLimits["ExpN-1"]
	   >> asymptoticLimits["ExpN-2"]) {
      setLimit(mass, width, false, true, 0, asymptoticLimits["Obs"]);
      for (int i_n = -2; i_n <= 2; i_n++) {
	setLimit(mass, width, true, true, i_n, 
		 asymptoticLimits[Form("ExpN%d",i_n)]);
      }
    }
    limitInput.close();
    return true;
  }
  else return false;
}

/**
   -----------------------------------------------------------------------------
   Calculate the p0 for a given mass and width hypothesis.
   @param mass - Integer representing the resonance mass.
   @param width - Integer representing the resonance width * 1000.
   @param crossSection - Integer representing the cross-section * 1000.
   @param makeNew - True if from toy MC, false if from text file storage.
   @param asymptotic - True iff. asymptotic results are desired. 
   @return - True iff loaded successfully.
*/
bool StatScan::singleP0Test(int mass, int width, int crossSection, 
			    bool makeNew, bool asymptotic){
  printer(Form("StatScan::singleP0Test(mass=%d, width=%d, New=%d, asympt=%d)",
	       mass, width, (int)makeNew, (int)asymptotic), false);
  
  bool successful = true;
  double observedP0 = 1.0;
  double expectedP0 = 1.0;
  
  TString textFileNameP0 = asymptotic ? 
    Form("%s/asymptotic_p0_values_mass%d_width%d.txt", m_outputDir.Data(),
	 mass, width) :
    Form("%s/toy_p0_values_mass%d_width%d.txt", m_outputDir.Data(),
	 mass, width);
  
  //----------------------------------------//
  // Open p0 values from file:
  if (!makeNew) {
    //printer(Form("StatScan::singleP0Test: Loading p0 from file %s",
    //		 textFileNameP0.Data()), false);
    
    // Open the saved p0 values from toys:
    std::ifstream inputFile(textFileNameP0);
    if (inputFile.is_open()) {
      while (inputFile >> observedP0 >> expectedP0) {
	//std::cout << "observed p0 = " << observedP0 << std::endl;
	//std::cout << "expected p0 = " << expectedP0 << std::endl;
      }
    }
    else {
      printer(Form("StatScan: ERROR opening text file %s",
		   textFileNameP0.Data()), false);
      successful = false;
    }
    inputFile.close();
  }
  
  //----------------------------------------//
  // Calculate new p0 value.
  else {
    printer("StatScan::singleP0Test: Calculating p0 from scratch", false);
    
    std::ofstream outputFile(textFileNameP0);    
    
    // Set the dataset to fit:
    //TString datasetToFit = m_config->getStr("WorkspaceObsData");
    
    // Load the workspace:
    // First check for local copy, then go to central copy.
    TString workspaceFileName = m_config->getStr("WorkspaceFile");
    TObjArray *array = workspaceFileName.Tokenize("/");
    
    workspaceFileName
      = ((TObjString*)array->At(array->GetEntries()-1))->GetString();
    TFile *wsFile = new TFile(workspaceFileName);
    if (!wsFile->IsZombie()) {
      printer(Form("StatScan: Loaded ws from %s", workspaceFileName.Data()),
	      false);
    }
    else {
      workspaceFileName = m_config->getStr("WorkspaceFile");
      printer(Form("StatScan: Load ws from %s", workspaceFileName.Data()),
	      false);
      wsFile = new TFile(workspaceFileName, "read");
    }
    RooWorkspace *workspace
      = (RooWorkspace*)wsFile->Get(m_config->getStr("WorkspaceName"));
    
    // Instantiate the test statistic class for calculations and plots:
    TestStat *testStat = new TestStat(m_configFileName, "new", workspace);
    testStat->setNominalSnapshot(m_config->getStr("WorkspaceSnapshotMu1"));
    if (!testStat->theWorkspace()->data(m_dataNameForFits) && m_dataToFit) {
      testStat->theWorkspace()->import(*m_dataToFit);
    }
    
    // Strategy settings (if defined):
    if (m_config->isDefined("FitOptions")) {
      testStat->setFitOptions(m_config->getStr("FitOptions"));
    }
  
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
	printer(Form("StatScan: Workspace has no variable %s",
		     listPoI[i_p].Data()), true);
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
    
    // Map of names and values of PoIs to set for fit:
    std::map<TString,double> mapPoI; mapPoI.clear();
    mapPoI[m_config->getStr("PoIForNormalization")]
      = ((double)crossSection/1000.0);
    mapPoI[m_config->getStr("PoIForMass")] = (double)mass;
    
    if ((m_config->getStr("AnalysisType")).Contains("Scalar") &&
	!(m_config->getStr("PoIForWidth")).EqualTo("GoM")) {
      if (width == 0) {
	mapPoI[m_config->getStr("PoIForWidth")] = 0.004;
      }
      else {
	mapPoI[m_config->getStr("PoIForWidth")] = 0.001 * (double)(mass*width);
      }
    }
    else {
      mapPoI[m_config->getStr("PoIForWidth")] = (((double)width)/1000.0);
    }
    
    //----------------------------------------//
    // Get the asymptotic p0 results:
    if (asymptotic) {
      testStat->useTwoSidedTestStat(m_config->getBool("UseTwoSided"));
      std::vector<double> asymptoticP0Values
 	= testStat->asymptoticP0(mapPoI, m_dataNameForFits,
				 m_config->getStr("WorkspaceSnapshotMu1"),
				 m_config->getStr("PoIForNormalization"));
      observedP0 = asymptoticP0Values[0];
      expectedP0 = asymptoticP0Values[1];
      successful = testStat->fitsAllConverged();
      delete testStat;
      delete workspace;
    }
    
    //----------------------------------------//
    // Get p0 results using toy MC:
    else {
      // Step 1: get observed q0 (no expected!)
      
      // Also force the mass and width to be constant always:
      testStat->setParam(m_config->getStr("PoIForMass"), (double)mass, true);
      
      
      if ((m_config->getStr("AnalysisType")).Contains("Scalar") &&
	  !(m_config->getStr("PoIForWidth")).EqualTo("GoM")) {
	if (width == 0) {
	  testStat->setParam(m_config->getStr("PoIForWidth"), 0.004, true);
	}
	else {
	  testStat->setParam(m_config->getStr("PoIForWidth"), 
			     0.001 * (double)(mass*width), true);
	}
      }
      else {
	testStat->setParam(m_config->getStr("PoIForWidth"), 
			   (((double)width)/1000.0), true);
      }
      
      // Perform the mu=0 fit, make sure normalization PoI is set to zero:
      mapPoI[m_config->getStr("PoIForNormalization")] = 0.0;
      double nllObsMu0 = testStat->getFitNLL(m_dataNameForFits, 0, true, 
					     mapPoI, true);
      
      // Perform mu-free fit (necessary for q0 calculation) and get signal:
      double nllObsMuFree = testStat->getFitNLL(m_dataNameForFits, 1, false,
						mapPoI, true);
      std::map<std::string,double> poiFromFit = testStat->getPoIs();
      double obsXSValue
	= poiFromFit[(std::string)m_config->getStr("PoIForNormalization")];
      
      // Calculate q0 and asymptotic p0:
      double q0Observed
	= testStat->getQ0FromNLL(nllObsMu0, nllObsMuFree, obsXSValue);
      double p0Asymptotic = testStat->getP0FromQ0(q0Observed);
      
      // Delete pointers and close files:  
      delete testStat;
      delete workspace;
      wsFile->Close();
      
      //----------------------------------------//
      // Process the toy MC files to obtain the CL values:
      printer("StatScan: Computing p0 from toy MC", false);
      
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
	printer("StatScan: ERROR with toy scan option.", true);
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
      
      observedP0 = toyAna->calculateP0FromToy(q0Observed);
      double q0Expected = 0.0;//toyAna->calculateExpQ0();

      expectedP0 = toyAna->calculateP0FromToy(q0Expected);
      
      successful = (toyAna->areInputFilesOK() && observedP0 > 0.0 &&
		    observedP0 < 1.0);
      delete toyAna;
    }
    
    // Save the p0 for this mass and width:
    outputFile << observedP0 << " " << expectedP0 << std::endl;
    outputFile.close();
  }
  
  // Store the p0 for this mass and width in the class:
  setP0(mass, width, false, asymptotic, observedP0);
  setP0(mass, width, true, asymptotic, expectedP0);
  
  // Then print to screen:
  std::cout << "StatScan: P0 Results:" << std::endl;
  std::cout << "\tobserved p0 = " << observedP0 << std::endl;
  std::cout << "\texpected p0 = " << expectedP0 << "\n" << std::endl;
  
  // Delete pointers, close files, return:
  printer(Form("StatScan::singleP0Test finished mass %d width %d\n", 
	       mass, width), false);
  
  return successful;
}

/**
   -----------------------------------------------------------------------------
   Check if a value belongs to a vector.
   @param theVector - The vector that might contain the value.
   @param theValue - The value to check for membership in the vector.
   @return - True iff the vector already contains the value.
*/
bool StatScan::vectorContainsValue(std::vector<int> theVector, int theValue) {
  for (int i_v = 0; i_v < (int)theVector.size(); i_v++) {
    if (theVector[i_v] == theValue) return true;
  }
  return false;
}

/**
   -----------------------------------------------------------------------------
   Set the list of mass points to use. The vector will be sorted automatically.
   @param massValues - The new vector of mass values to use.
*/
void StatScan::useTheseMasses(std::vector<int> massValues) {
  m_massValues = massValues;
  std::sort(m_massValues.begin(), m_massValues.end());
}

/**
   -----------------------------------------------------------------------------
   Set the list of width points to use. The vector will be sorted automatically.
   @param widthValues - The new vector of width values to use.
*/
void StatScan::useTheseWidths(std::vector<int> widthValues) {
  m_widthValues = widthValues;
  std::sort(m_widthValues.begin(), m_widthValues.end());
}

/**
   -----------------------------------------------------------------------------
   Set the list of cross-section points to use. The vector will be sorted 
   automatically.
   @param xsValues - The new vector of cross-section values to use.
*/
void StatScan::useTheseXS(std::vector<int> xsValues) {
  m_xsValues = xsValues;
  std::sort(m_xsValues.begin(), m_xsValues.end());
}
