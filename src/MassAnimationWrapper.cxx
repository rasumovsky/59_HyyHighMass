////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Name: MassAnimationWrapper.cxx                                            //
//                                                                            //
//  Creator: Andrew Hard                                                      //
//  Email: ahard@cern.ch                                                      //
//  Date: 18/04/2016                                                          //
//                                                                            //
//  A wrapper for the mass animation class for cluster jobs that contain sub- //
//  sets of the frames.                                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CommonFunc.h"
#include "CommonHead.h"
#include "Config.h"
#include "MassAnimation.h"
#include <time.h>

/**
   -----------------------------------------------------------------------------
   The main method plots the p0 or 95CL values as a function of mass.
   @param configFile - The analysis configuration file.
   @param options - Job options. Can be "toy" or "asymptotic" or "both"
*/
int main(int argc, char **argv) {
  
  // Check that arguments are provided.
  if (argc < 5) {
    std::cout << "\nUsage: " << argv[0] 
	      << " <configFile> <options> <minFrame> <maxFrame>" << std::endl;
    exit(0);
  }
  
  // Clock the animation:
  clock_t time;
  time = clock();
  
  TString configFile = argv[1];
  TString options = argv[2];
  int minFrame = atoi(argv[3]);
  int maxFrame = atoi(argv[4]);
  
  // Load the analysis configuration file:
  Config *config = new Config(configFile);

  MassAnimation *animation = new MassAnimation(configFile, options);
  animation->setNFrames(config->getInt("AnimationFrames"));
  animation->getDataForFrames();
  for (int i_f = minFrame; i_f < maxFrame; i_f++) {
    animation->makeSingleFrame(i_f);
  }
  
  // Clock the animation:
  time = clock() - time;
  if (config->getBool("Verbose")) {
    std::cout << "MassAnimationWrapper: Finished frames " << minFrame << " - "
	      << maxFrame << std::endl;
    printf("\tRequired %d clock cycles (%f seconds).\n\n",
	   (int)time, ((float)time/CLOCKS_PER_SEC));
  }
  return 0;
}
