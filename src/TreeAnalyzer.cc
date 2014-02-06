#include <sstream>
#include <fstream>
#include <iterator>
#include <string>

#include "TreeAnalyzer.h"
//#include "SummaryAnalyzer.h"
#include "ObjectMessenger.h"

#include "boost/functional/hash.hpp"

#include "TFile.h"
#include "TH1D.h"
#include "TProofOutputFile.h"
#include "TTree.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TreeAnalyzer::TreeAnalyzer(const std::string & aName,
				       const std::string & cfgFileName, 
				       TProofOutputFile * proofFile)
  :Analyzer(aName){

  ///Analysis control
  nEventsToAnalyze_ = 0;

  cfgFileName_ = cfgFileName;

  parseCfg(cfgFileName_);

  if(proofFile){   
    proofFile->SetOutputFileName((filePath_+"/PFAnalysis_"+sampleName_+".root").c_str());
    store_ = new TFileService(proofFile->OpenFile("RECREATE"));
  }
  else{ 
    // Create histogram store
    store_ = new TFileService((filePath_+"/PFAnalysis_"+sampleName_+".root").c_str());
  }

  ///Histogram with processing statistics. Necessary for the PROOF based analysis
  hStats_ = store_->mkdir("Statistics").make<TH1D>("hStats","Various statistics",21,-0.5,20.5);
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

  myStrSelections_ = new pat::strbitset(); 

  myObjMessenger_ = new ObjectMessenger("ObjMessenger");

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TreeAnalyzer::~TreeAnalyzer(){

  std::cout<<"TreeAnalyzer::~TreeAnalyzer() Begin"<<std::endl;

  //delete mySummary_;
  delete store_;

  std::cout<<"TreeAnalyzer::~TreeAnalyzer() Done"<<std::endl;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::scaleHistograms(){

  float weight = 1.0;

  for(unsigned int i=0;i<myAnalyzers_.size();++i){
    std::string name = myAnalyzers_[i]->name();  
    TDirectoryFile* summary = (TDirectoryFile*)store_->file().Get(name.c_str());
    if(!summary){
      std::cout<<"Histogram directory for analyzer: "<<name.c_str()
	       <<" not found!"<<std::endl;
      continue;
    }
    TList *list = summary->GetList();
    TIter next(list);
    TObject *obj = 0;
    while ((obj = next())){
      if(obj->IsA()->InheritsFrom("TH1")){ 
	TH1 *h = (TH1*)summary->Get(obj->GetName());
	h->Scale(weight);
      }
      if(obj->IsA()->InheritsFrom("TDirectory")){ 
	TDirectory* aDir = (TDirectory*)summary->Get(obj->GetName());
	TList *listSubDir = aDir->GetList();
	TIter next2(listSubDir);
	TObject *obj2 = 0;
	while ((obj2 = next2())){
	  if(obj2->IsA()->InheritsFrom("TH1")){ 
	    TH1 *h1 = (TH1*)aDir->Get(obj2->GetName());
	    h1->Scale(weight);
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
  
  boost::hash<std::string> string_hash;
  std::ifstream in((filePath_+"/PFAnalysis_"+sampleName_+"_SelEvents.dat").c_str());
  std::string separator;
  int counter = 0;
  unsigned int run,ls,event;
  char text[500];
  while (in.good() && !in.eof()){
    std::string tmpString;
    in>>run>>ls>>event;
    sprintf(text,"%d:%d:%d",run,ls,event);
    tmpString.append(text);
    eventsToProcessHash_[tmpString] = true;    
  }
  in.close();

  std::cout<<"Event tag file name: "<<filePath_+"PFAnalysis_"+sampleName_+"_SelEvents.dat"<<std::endl;
  std::cout<<"Number of events to process with event list: "<<eventsToProcessHash_.size()<<std::endl
	   <<"If 0, all events will be processed."<<std::endl;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  TreeAnalyzer::init(std::vector<Analyzer*> myAnalyzers){
/*
  myAnalyzers_ = myAnalyzers;
  mySummary_ = new SummaryAnalyzer("Summary");
  myAnalyzers_.push_back(mySummary_);

  for(unsigned int i=0;i<myAnalyzers_.size();++i){ 
    myDirectories_.push_back(store_->mkdir(myAnalyzers_[i]->name()));
    myAnalyzers_[i]->initialize(parameterSet_->getParameter<edm::ParameterSet>(myAnalyzers_[i]->name()), 
			        myDirectories_[myDirectories_.size()-1],
				myStrSelections_);
  }

 for(unsigned int i=0;i<myAnalyzers_.size();++i){
   myAnalyzers_[i]->addBranch(mySummary_->getTree());  
   myAnalyzers_[i]->addCutHistos(mySummary_->getHistoList());  
 }

 myDirectories_.push_back(store_->mkdir("RunInfoAccounting"));
 runInfoAccounting_ = new RunInfoAccounting( myDirectories_.back() ,
					    "RunInfoAccounting" );
*/
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  TreeAnalyzer::finalize(){

  for(unsigned int i=0;i<myAnalyzers_.size();++i) myAnalyzers_[i]->finalize();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int TreeAnalyzer::loop(){
/*
  fwlite::ChainEvent ev(fileNames_);
  std::cout<<"Events total: "<<ev.size()<<std::endl;

  nEventsAnalyzed_ = 0;
  nEventsSkipped_ = 0;
  int eventPreviouslyPrinted=-1;
  TFile* currentFile = 0;
  ///////
   for(ev.toBegin();
       !ev.atEnd() && (nEventsToAnalyze_<0 || (nEventsAnalyzed_+nEventsSkipped_)<nEventsToAnalyze_); ++ev){
     
     //std::cout<<ev.id()<<std::endl;
     if((( nEventsAnalyzed_ < 10) ||
	 nEventsAnalyzed_%10000==0) &&  nEventsAnalyzed_ != eventPreviouslyPrinted ) {
       eventPreviouslyPrinted = nEventsAnalyzed_;
     }
     //analyze the Run data for every new run, and also a new file
     if(ev.id().run()!= currentRun_) ++nRunsAnalyzed_;
     analyze(ev,ev.getRun());
   }
   
   std::cout << "Events skipped: " << nEventsSkipped_ << std::endl ;
   return nEventsAnalyzed_;
   */
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool TreeAnalyzer::analyze(const EventBase& iEvent){

  clear();
    ///////
    for(unsigned int i=0;i<myAnalyzers_.size();++i){
      ///If analyzer returns false, skip to the last one, the Summary, unless filtering is disabled for this analyzer.
      if(!myAnalyzers_[i]->analyze(iEvent,myObjMessenger_) && myAnalyzers_[i]->filter() && myAnalyzers_.size()>1) i = myAnalyzers_.size()-2;
    }
    ///Clear all the analyzers, even if it was not called in this event.
    ///Important for proper TTree filling.
    for(unsigned int i=0;i<myAnalyzers_.size();++i) myAnalyzers_[i]->clear(); 
        
    myObjMessenger_->clear();
  
  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::clear(){

  myStrSelections_->set(false);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
