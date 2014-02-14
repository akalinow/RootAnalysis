/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : Artur Kalinowski
 Version     :
 Copyright   : Your copyright notice
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/libRootAnalysis/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include <iostream>

#include "TreeAnalyzer.h"
#include "OTFAnalyzer.h"
#include "EventProxyOTF.h"

#include "TFile.h"
#include "TStopwatch.h"

int main(int argc, char ** argv) {

	std::string cfgFileName = "cfg.py";

	  if(argc<2){
	    std::cout<<"Usage: readEvents cfg.py"<<std::endl;
	    return 1;
	  }
	  else cfgFileName = argv[1];


	std::cout<<"Start"<<std::endl;
	TStopwatch timer;
	timer.Start();
	  //----------------------------------------------------------
	 std::vector<Analyzer*> myAnalyzers;
	 EventProxyOTF *myEvent = new EventProxyOTF();

	  myAnalyzers.push_back(new OTFAnalyzer("OTFAnalyzer"));

	  TreeAnalyzer *tree = new TreeAnalyzer("TreeAnalyzer",cfgFileName, myEvent);
	  tree->init(myAnalyzers);
	  int nEventsAnalysed = tree->loop();
	  tree->finalize();

	  timer.Stop();
	  Double_t rtime = timer.RealTime();
	  Double_t ctime = timer.CpuTime();
	  printf("Analysed events: %d \n",nEventsAnalysed);
	  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
	  printf("%4.2f events / RealTime second .\n", nEventsAnalysed/rtime);
	  printf("%4.2f events / CpuTime second .\n", nEventsAnalysed/ctime);

	  tree->scaleHistograms();
	  for(unsigned int i=0;i<myAnalyzers.size();++i) delete myAnalyzers[i];
	  delete tree;
	  delete myEvent;

	  std::cout<<"Done"<<std::endl;
	  return 0;


}
