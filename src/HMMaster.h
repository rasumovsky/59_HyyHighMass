////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: HMMaster.h                                                          //
//                                                                            //
//  Created: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 25/02/2016                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "Config.h"
#include "TestStat.h"
#include "ToyAnalysis.h"

bool m_isFirstJob;

Config *m_config;

void submitTSViaBsub(TString exeConfigFile, TString exeOption,
		     TString exeSignal);

void submitMLViaBsub(TString exeConfigFile, TString exeOption, 
		     TString exeSignal);

void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     int exeSeed, int exeToysPerJob, int resonanceMass);
