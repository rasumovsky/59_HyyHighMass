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
#include "TestStat.h"
#include "ToyAnalysis.h"

bool m_isFirstJob;

Config *m_config;

void submitPEViaBsub(TString exeConfigFile, TString exeOption,
		     int exeSeed, int exeToysPerJob);
