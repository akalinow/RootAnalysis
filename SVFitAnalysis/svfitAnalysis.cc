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
#include "svfitAnalyzer.h"
#include "MLAnalyzer.h"
#include "MLObjectMessenger.h"

#include "TFile.h"
#include "TStopwatch.h"

#include "boost/functional/hash.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"

#include "TROOT.h"
#include "TObject.h"


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

	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(cfgFileName, pt);
	std::string processName = pt.get<std::string>("TreeAnalyzer.processName","Test");
	unsigned noOfThreads = pt.get("TreeAnalyzer.threads",1);

	//Tell Root we want to be multi-threaded
	ROOT::EnableThreadSafety();
	//When threading, also have to keep ROOT from logging all TObjects into a list
	TObject::SetObjectStat(false);

	//----------------------------------------------------------
	 std::vector<Analyzer*> myAnalyzers;
	 EventProxyHTT *myEvent = new EventProxyHTT();

	 std::string decayModeName;
	 if(processName=="AnalysisMuTau") decayModeName = "MuTau";
	 else if(processName=="AnalysisTauTau") decayModeName = "TauTau";
	 else if(processName=="AnalysisMuMu") decayModeName = "MuMu";
	 else{
	   std::cout<<"Incorrect process name: "<<processName<<std::endl;
	   return 1;
	 }

	 if(processName.find("Analysis")!=std::string::npos) {
	   myAnalyzers.push_back(new svfitAnalyzer("svfitAnalyzer",decayModeName));
				if(noOfThreads==1) 
					myAnalyzers.push_back(new MLAnalyzer("MLAnalyzer",decayModeName));
   }
	 else{
	   std::cout<<"Incorrect process name: "<<processName<<std::endl;
	   return 1;
	 }

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
