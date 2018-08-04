#include <sstream>
#include <fstream>
#include <iterator>
#include <string>
#include <exception>
#include <stdexcept>

#include <omp.h>
#include "TreeAnalyzer.h"
#include "SummaryAnalyzer.h"
#include "ObjectMessenger.h"
#include "EventProxyBase.h"

#include "EventProxyBase.h"
#include "AnalysisHistograms.h"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/ini_parser.hpp"
#include "boost/tokenizer.hpp"
#include "boost/functional/hash.hpp"

#include "TFile.h"
#include "TH1D.h"
#include "TProofOutputFile.h"
#include "TTree.h"
#include "TChain.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::setObjectMessenger(ObjectMessenger* mess)
{
  myObjMessenger_=mess;
}

TreeAnalyzer::TreeAnalyzer(const std::string & aName,
			   const std::string & cfgFileName,
			   EventProxyBase *aProxy,
			   TProofOutputFile * proofFile) : Analyzer(aName){

  ///Analysis control
  nEventsToAnalyze_ = 0;

  cfgFileName_ = cfgFileName;

  parseCfg(cfgFileName_);

  if(proofFile) {
    proofFile->SetOutputFileName((filePath_+"/RootAnalysis_"+sampleName_+".root").c_str());
    store_ = proofFile->OpenFile("RECREATE");
  }
  else{
    // Create histogram store
    std::string fullPath = filePath_+"/RootAnalysis_"+sampleName_+".root";
    store_ = new TFile(fullPath.c_str(),"RECREATE");
  }

  ///Histogram with processing statistics. Necessary for the PROOF based analysis
  hStats_ = new TH1D("hStats","Various statistics",21,-0.5,20.5);
  hStats_->GetXaxis()->SetBinLabel(1,"Xsection");
  hStats_->GetXaxis()->SetBinLabel(2,"external scaling factor");
  hStats_->GetXaxis()->SetBinLabel(3,"number of runs");
  hStats_->GetXaxis()->SetBinLabel(4,"number of parrarel nodes");
  hStats_->GetXaxis()->SetBinLabel(5,"number of event analyzed");
  hStats_->GetXaxis()->SetBinLabel(6,"number of event skipped");
  hStats_->GetXaxis()->SetBinLabel(7,"generator preselection eff.");
  hStats_->GetXaxis()->SetBinLabel(8,"number of events processed at RECO/AOD");
  hStats_->GetXaxis()->SetBinLabel(9,"number of events saved from RECO/AOD");
  hStats_->GetXaxis()->SetBinLabel(10,"reco preselection eff.");
  TDirectory *aDir = store_->mkdir("Statistics");
  hStats_->SetDirectory(aDir);

  myStrSelections_ = new pat::strbitset();

  myObjMessenger_ = new ObjectMessenger("Default ObjMessenger Created In TreeAnalyzer.h");

  myProxy_ = aProxy;
  mySummary_ = 0;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TreeAnalyzer::~TreeAnalyzer(){

  std::cout<<"TreeAnalyzer::~TreeAnalyzer() Begin"<<std::endl;

  if(mySummary_) delete mySummary_;
  if(store_) {
    store_->Write();
    delete store_;
  //TODO: delete myObjMessenger ???
  }

  std::cout<<"TreeAnalyzer::~TreeAnalyzer() Done"<<std::endl;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::scaleHistograms(){

  float weight = 1.0;

  for(unsigned int i=0; i<myAnalyzers_.size(); ++i) {
    std::string name = myAnalyzers_[i]->name();
    TDirectory* summary = (TDirectory*)store_->Get(name.c_str());
    if(!summary) {
      std::cout<<"Histogram directory for analyzer: "<<name.c_str()
	       <<" not found!"<<std::endl;
      continue;
    }
    TList *list = summary->GetList();
    TIter next(list);
    TObject *obj = 0;
    while ((obj = next())) {
      if(obj->IsA()->InheritsFrom("TH1")) {
	TH1 *h = (TH1*)summary->Get(obj->GetName());
	if(h) h->Scale(weight);
      }
      if(obj->IsA()->InheritsFrom("TDirectory")) {
	TDirectory* aDir = (TDirectory*)summary->Get(obj->GetName());
	TList *listSubDir = aDir->GetList();
	TIter next2(listSubDir);
	TObject *obj2 = 0;
	while ((obj2 = next2())) {
	  if(obj2->IsA()->InheritsFrom("TH1")) {
	    TH1 *h1 = (TH1*)aDir->Get(obj2->GetName());
	    if(h1) h1->Scale(weight);
	  }
	}
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::parseCfg(const std::string & cfgFileName){

  eventWeight_ = 1.0;

  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(cfgFileName, pt);

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(", ");
  std::string str = pt.get<std::string>("TreeAnalyzer.inputFiles");
  std::string dataPath = pt.get<std::string>("TreeAnalyzer.dataPath","");
  tokenizer tokens(str, sep);
  std::cout<<"Reading files: "<<std::endl;
  std::string fileName;
  for (auto it: tokens) {
    fileName = dataPath+"/"+it;
    std::cout<<fileName<<std::endl;
    fileNames_.push_back(fileName);
  }

  filePath_ = pt.get<std::string>("TreeAnalyzer.outputPath");
  sampleName_ = pt.get<std::string>("TreeAnalyzer.processName","Test");

  nEventsToAnalyze_ = pt.get("TreeAnalyzer.eventsToAnalyze",-1);
  nEventsToPrint_ = pt.get("TreeAnalyzer.eventsToPrint",100);
  nThreads_ = pt.get("TreeAnalyzer.threads",1);

  if(nThreads_>AnalysisHistograms::maxThreads-1) {
    std::cout<<"Number of threads exceeds maximum value of "
	     << AnalysisHistograms::maxThreads-1
	     <<" will use the maximum value"<<std::endl;
    nThreads_ = AnalysisHistograms::maxThreads-1;
  }
  omp_set_num_threads(nThreads_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::init(std::vector<Analyzer*> myAnalyzers){

  myProxy_->init(fileNames_);
  myAnalyzers_ = myAnalyzers;

  if(nThreads_==1) {
    mySummary_ = new SummaryAnalyzer("Summary");
    myAnalyzers_.push_back(mySummary_);
  }

  try
  {
    for(unsigned int i=0; i<myAnalyzers_.size(); ++i) 
    {
      std::string analyzerName = myAnalyzers_[i]->name();
      TDirectory *analyzerDir = store_->mkdir(analyzerName.c_str());
      myAnalyzers_[i]->initialize(analyzerDir, myStrSelections_);
    }
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN TreeAnalyzer::init! PROBABLY SOMETHING WRONG WITH THE NAME OF ANALYZER!"));
  }



  for(int iThread=0; iThread<omp_get_max_threads(); ++iThread) {
    myProxiesThread_[iThread] = myProxy_->clone();
    myProxiesThread_[iThread]->init(fileNames_);
    for(unsigned int iAnalyzer=0; iAnalyzer<myAnalyzers_.size(); ++iAnalyzer) {
      if(iThread==0) myAnalyzersThreads_[iThread].push_back(myAnalyzers_[iAnalyzer]);
      else myAnalyzersThreads_[iThread].push_back(myAnalyzers_[iAnalyzer]->clone());
    }
  }

  ///Tree making does not work with multithread.
  if(nThreads_==1) {
    for(unsigned int i=0; i<myAnalyzers_.size(); ++i) {
      myAnalyzers_[i]->addBranch(mySummary_->getTree());
      myAnalyzers_[i]->addCutHistos(mySummary_->getHistoList());
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::finalize(){
  for(unsigned int i=0; i<myAnalyzers_.size(); ++i) myAnalyzers_[i]->finalize();
}
//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
int TreeAnalyzer::loop(){

  std::cout<<"Events total: "<<myProxy_->size()<<std::endl;
  unsigned int printoutStep = myProxy_->size()*0.1;
  if(!printoutStep) printoutStep = 1;
  TH1::AddDirectory(kFALSE);
  nEventsAnalyzed_ = 0;
  nEventsSkipped_ = 0;
  if(nEventsToAnalyze_<0 || nEventsToAnalyze_>myProxy_->size()) nEventsToAnalyze_ = myProxy_->size();
 
  myProxiesThread_[0]->toBegin();

  unsigned int eventCount[nThreads_] {0};

#pragma omp parallel for schedule(dynamic,500)
  for(int aEvent=0; aEvent<nEventsToAnalyze_; ++aEvent) {
    if(aEvent< nEventsToPrint_ || aEvent%printoutStep==0)
      std::cout<<"Events analyzed: "<<aEvent<<"/"<<nEventsToAnalyze_
	       <<" ("<<(float)aEvent/nEventsToAnalyze_<<")"
	       <<" thread: "<<omp_get_thread_num()<<std::endl;
    analyze(myProxiesThread_[omp_get_thread_num()]->toN(aEvent));
    ++eventCount[omp_get_thread_num()];
  }

  for(auto it: eventCount) nEventsAnalyzed_+=it;

  return nEventsAnalyzed_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool TreeAnalyzer::analyze(const EventProxyBase& iEvent){

  bool decision = true;
  for(unsigned int i=0; i<myAnalyzers_.size(); ++i) {
    if(!decision) break;
    decision &= myAnalyzersThreads_[omp_get_thread_num()][i]->analyze(iEvent,myObjMessenger_);
  }

  return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::clear(){
  myStrSelections_->set(false);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
