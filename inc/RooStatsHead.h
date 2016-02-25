////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: RooStatsHead.h                                                      //
//  Creator: Hongtao Yang                                                     //
//                                                                            //
//  All of the includes necessary for RooStats and HistFactory.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ROOSTATHEAD_H
#define ROOSTATHEAD_H
/// RooStat headers

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/NeymanConstruction.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestResult.h"

#include "RooStats/NumberCountingUtils.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/SPlot.h"

#include "RooStats/FrequentistCalculator.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
//#include "RooStats/ToyMCImportanceSampler.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/LinInterpVar.h"
#include "RooStats/HistFactory/HistoToWorkspaceFactoryFast.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/ParamHistFunc.h"

//#include "RooStats/HistFactory/Systematics.h"

#endif
