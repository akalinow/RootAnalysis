#include <iostream>
#include <string>

#include "TFile.h"
#include "TStopwatch.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "PFAnalyses/CommonTools/interface/FWLiteTreeAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteTriggerAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteJetVetoAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteDiTauAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/FWLiteQCDAnalyzer.h"


int main(int argc, char ** argv){


  std::string cfgFileName = "ps.cfg";

  if(argc<2){
    std::cout<<"Usage: readEvents cfg.py"<<std::endl;
    return 1;
  }
  else cfgFileName = argv[1];

  std::cout<<"Start"<<std::endl;
  TStopwatch timer;
  timer.Start();
  //----------------------------------------------------------
  AutoLibraryLoader::enable();

  std::vector<FWLiteAnalyzer*> myAnalyzers;

  //myAnalyzers.push_back(new FWLiteDiTauAnalyzer("DiTau"));
  //myAnalyzers.push_back(new FWLiteTriggerAnalyzer("TriggerAnalyzer"));
  //myAnalyzers.push_back(new FWLiteJetVetoAnalyzer("JetVeto"));
 
  myAnalyzers.push_back(new FWLiteTriggerAnalyzer("TriggerAnalyzer"));
  myAnalyzers.push_back(new FWLiteQCDAnalyzer("QCDAnalysis"));

  
  FWLiteTreeAnalyzer *tree = new FWLiteTreeAnalyzer("TreeAnalyzer",cfgFileName);
  tree->init(myAnalyzers);
  int nEventsAnalysed = tree->loop();
  tree->finalize();
  //----------------------------------------------------------
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

  std::cout<<"Done"<<std::endl;
  return 0;
}
