////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: HMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 27/02/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "Config.h"
#include "MassAnimation.h"
#include "TestStat.h"
#include "ToyAnalysis.h"
#include "Workspace.h"

bool m_isFirstJob;

Config *m_config;

void submitGIFViaBsub(TString exeConfigFile, TString exeOption, int exeMinFrame,
		      int exeMaxFrame);

void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     int exeSeed, int exeToysPerJob);

void submitStatViaBsub(TString exeConfigFile, TString exeOption, int exeWidth,
		       int exeMassMin, int exeMassMax, int exeMassStep);
