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
#include "MLAnalyzer.h"
#include "MLObjectMessenger.h"

#include "HTTSynchNTuple.h"

#include "TFile.h"
#include "TStopwatch.h"

#include "boost/functional/hash.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"

#include "TROOT.h"
#include "TObject.h"


/*! 
*	Prints the explanatory string of an exception. If the exception is nested,
*	recurses to print the explanatory of the exception it holds
*/
void print_exception(const std::exception& e, int level =  0)
{
    std::cerr << std::string(level, ' ') << "exception: " << e.what() << '\n';
    try
    {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e)
    {
        print_exception(e, level+1);
    } 
    catch(...) {}
}

int main(int argc, char ** argv)
{
  try
    {
      std::string cfgFileName = "cfg.ini";

      if(argc<2){
	std::cout<<"Usage: readEvents cfg.init"<<std::endl;
	return 1;
      }
      else 
	{
	  cfgFileName = argv[1];
	}


      std::cout<<"Start"<<std::endl;
      TStopwatch timer;
      timer.Start();

      boost::property_tree::ptree pt;
      boost::property_tree::ini_parser::read_ini(cfgFileName, pt);
      std::string processName = pt.get<std::string>("TreeAnalyzer.processName","Test");
      unsigned noOfThreads = pt.get("TreeAnalyzer.threads",1);
      if(noOfThreads != 1)
	std::cerr<<"[WARNING] Number of threads != 1. TTree output is available only in single thread mode."<<std::endl;
      //Tell Root we want to be multi-threaded
      ROOT::EnableThreadSafety();
      ROOT::DisableImplicitMT();
      //When threading, also have to keep ROOT from logging all TObjects into a list
      TObject::SetObjectStat(false);

      //----------------------------------------------------------
      std::vector<Analyzer*> myAnalyzers;
      EventProxyHTT *myEvent = new EventProxyHTT();
      ObjectMessenger* OMess = new ObjectMessenger("ObjectMessenger created in HTTAnalysis.cc");

      std::string decayModeName;
      if(processName.find("MuTau")!=std::string::npos) decayModeName = "MuTau";
      else if(processName.find("TT")!=std::string::npos) decayModeName = "TauTau";
      else if(processName.find("MM")!=std::string::npos) decayModeName = "MuMu";     
      else{
	std::cout<<"Incorrect process name: "<<processName<<std::endl;
	return 1;
      }
            
      if(processName.find("Analysis")!=std::string::npos){
	myAnalyzers.push_back(new HTTAnalyzer("HTTAnalyzer",decayModeName));
	// MLAnalyzer's purpose is to flush data to TTree, but this kind of output is disabled in multithread mode
	if(processName=="MLAnalysisMuTau" and noOfThreads==1) {
	  delete OMess;
	  OMess = new MLObjectMessenger("MLObjectMessenger created in HTTAnalysis.cc");
	  myAnalyzers.push_back(new MLAnalyzer("MLAnalyzer",decayModeName));
	}
      }
      else if(processName.find("Synch")!=std::string::npos)
	myAnalyzers.push_back(new HTTSynchNTuple("SynchNTuple",decayModeName));
      else{
	std::cout<<"Incorrect process name: "<<processName<<std::endl;
	return 1;
      }

      TreeAnalyzer *tree = new TreeAnalyzer("TreeAnalyzer",cfgFileName, myEvent);	 	
      tree->setObjectMessenger(OMess);
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
  catch(const std::exception& e)
    {
      print_exception(e);
      return EXIT_FAILURE;
    }
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
