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
#include "HTTWeightsMaker.h"
#include "HTTAnalyzer.h"
#include "HTTSynchNTuple.h"
#include "HTTSynchNTupleTT.h"
#include "HTTSynchNTupleMM.h"

#include "TFile.h"
#include "TStopwatch.h"

#include "boost/functional/hash.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/tokenizer.hpp"

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
	
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	std::string processName = pt.get<std::string>("TreeAnalyzer.processName","Test");

	//Tell Root we want to be multi-threaded
	ROOT::EnableThreadSafety();
	//When threading, also have to keep ROOT from logging all TObjects into a list
	TObject::SetObjectStat(false);
	
	//----------------------------------------------------------
	 std::vector<Analyzer*> myAnalyzers;
	 EventProxyHTT *myEvent = new EventProxyHTT();

	 if(processName=="Weights" || processName=="PU") myAnalyzers.push_back(new HTTWeightsMaker("HTTWeightsMaker"));
	 else 
	 if(processName.find("SynchNTupleMT")!=std::string::npos) myAnalyzers.push_back(new HTTSynchNTuple("HTTSynchNTuple"));
	 if(processName.find("SynchNTupleTT") myAnalyzers.push_back(new HTTSynchNTupleTT("HTTSynchNTupleTT"));
	 if(processName.find("SynchNTupleMM") myAnalyzers.push_back(new HTTSynchNTupleMM("HTTSynchNTupleMM"));
	 else myAnalyzers.push_back(new HTTAnalyzer("HTTAnalyzer"));

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
