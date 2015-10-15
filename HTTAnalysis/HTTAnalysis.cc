/*
 ============================================================================
 Name        : HTTAnalysis.cc
 Author      : Artur Kalinowski
 Version     :
 Copyright   : GPL
 ============================================================================
 */

#include <iostream>

#include "TreeAnalyzer.h"

#include "EventProxyHTT.h"
#include "HTTAnalyzer.h"


#include "TFile.h"
#include "TStopwatch.h"

int main(int argc, char ** argv) {

	std::string cfgFileName = "cfg.ini";

	  if(argc<2){
	    std::cout<<"Usage: readEvents cfg.init"<<std::endl;
	    return 1;
	  }
	  else cfgFileName = argv[1];


	std::cout<<"Start"<<std::endl;
	TStopwatch timer;
	timer.Start();
	  //----------------------------------------------------------
	 std::vector<Analyzer*> myAnalyzers;
	 EventProxyHTT *myEvent = new EventProxyHTT();

	 myAnalyzers.push_back(new HTTAnalyzer("HTTAnalyzer"));

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
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
