////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: CommonHead.h                                                        //
//  Creator: Hongtao Yang                                                     //
//                                                                            //
//  This file contains most of the standard ROOT includes necessary for the   //
//  H->gg + DM search with 13 TeV data.                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef COMMONHEAD_H
#define COMMONHEAD_H

/// c++ headers
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cassert>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <list>
#include <map>
#include <cstdlib>
#include <cmath>

/// ROOT headers
#include <TApplication.h>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGaxis.h>
#include <TGenPhaseSpace.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <THStack.h>
#include "TImage.h"
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TPolyLine3D.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorT.h>
#include <TVersionCheck.h>

/// ROOT math headers
#include <Math/QuantFuncMathCore.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <Math/Factory.h>

/// TMVA headers
#include <TMVA/Tools.h>
#include <TMVA/Factory.h>
#include <TMVA/Reader.h>
#include <TMVA/MethodCuts.h>

#endif
