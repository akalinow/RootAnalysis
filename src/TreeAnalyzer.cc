#include <sstream>
#include <fstream>
#include <iterator>
#include <string>

//#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "TreeAnalyzer.h"
#include "SummaryAnalyzer.h"
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

  ///Cross section estimate 
  xSectionError_ = 0;
  xSectionWeight_ = 1.0;
  xSectionWMean_ = 1.0;
  preselectionEff_ = 1.0; 
  genPreselectionEff_ = 1.0;
  recoPreselectionEff_ = 1.0;
  initialRecoEvents_ = 0;
  finalRecoEvents_ = 0;
  ///
  eventWeight_ = 1.0;
  ///Analysis control
  nEventsToAnalyze_ = 0;
  nEventsAnalyzed_ = 0;
  nEventsSkipped_ = 0;
  currentRun_ = 0;
  nRunsAnalyzed_ = 0;
  runInfoAccounting_ = 0;

  cfgFileName_ = cfgFileName;

  parseCfg(cfgFileName_);

  if(proofFile){   
    proofFile->SetOutputFileName((filePath_+"/PFAnalysis_"+sampleName_+".root").c_str());
    store_ = new fwlite::TFileService(proofFile->OpenFile("RECREATE"));
  }
  else{ 
    // Create histogram store
    store_ = new fwlite::TFileService((filePath_+"/PFAnalysis_"+sampleName_+".root").c_str());
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

  std::cout<<"FWLiteTreeAnalyzer::~FWLiteTreeAnalyzer() Begin"<<std::endl;

  if(initialRecoEvents_) recoPreselectionEff_ = (float)finalRecoEvents_/initialRecoEvents_;
  
  hStats_->SetBinContent(1,xSectionWMean_);
  hStats_->SetBinContent(2,preselectionEff_);
  hStats_->SetBinContent(3,nRunsAnalyzed_);
  hStats_->SetBinContent(4,1.0);
  hStats_->SetBinContent(5,nEventsAnalyzed_);
  hStats_->SetBinContent(6,nEventsSkipped_);
  hStats_->SetBinContent(7,genPreselectionEff_);
  hStats_->SetBinContent(8,initialRecoEvents_);
  hStats_->SetBinContent(9,finalRecoEvents_);
  hStats_->SetBinContent(10,recoPreselectionEff_);

  std::cout.precision(6);
  std::cout<<"------------ Weighting report: "<<std::endl;
  std::cout<<"Cross section [pb]: "<<xSectionWMean_<<std::endl;
  std::cout<<"Generator level preselection eff: "<<genPreselectionEff_<<std::endl;
  std::cout<<"number of Runs analyzed: "<<nRunsAnalyzed_<<std::endl;
  std::cout<<"Number of events processed from RECO/AOD: "<<initialRecoEvents_<<std::endl;
  std::cout<<"Number of events saved at RECO/AOD: "<<finalRecoEvents_<<std::endl;  
  std::cout<<"Reco preselection efficiency: "<< recoPreselectionEff_<<std::endl;
  std::cout<<"External scaling: "<<preselectionEff_<<std::endl;
  if(preselectionEff_>0){
  std::cout<<"Full weight for histograms/cut counters: "
	   <<"xSection*genPreselectionEff*recoPreselectionEff*externalFactor/(nEventsAnalyzed+nEventsSkipped_) ="
	   <<xSectionWMean_*genPreselectionEff_*recoPreselectionEff_*preselectionEff_/(nEventsAnalyzed_+nEventsSkipped_)
	   <<std::endl;  
  }
  else{
    std::cout<<"Histograms are NOT prescaled due to negative external weight."<<std::endl;
  }  
  std::cout<<"The End -----------------------"<<std::endl;
  delete mySummary_;
  delete store_;
  delete runInfoAccounting_;

  std::cout<<"FWLiteTreeAnalyzer::~FWLiteTreeAnalyzer() Done"<<std::endl;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::scaleHistograms(){

  if(initialRecoEvents_) recoPreselectionEff_ = (float)finalRecoEvents_/initialRecoEvents_;

  float weight = xSectionWMean_*genPreselectionEff_*recoPreselectionEff_*preselectionEff_/(nEventsAnalyzed_+nEventsSkipped_);
  
  
  if(preselectionEff_<0) return;

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

  // Create the ParameterSet object from this configuration string.
  PythonProcessDesc builder(cfgFileName);
  parameterSet_ = builder.processDesc()->getProcessPSet();

  ///Get list of input files
  edm::ParameterSet source = parameterSet_->getParameter<edm::ParameterSet>("@main_input");
  fileNames_ = source.getUntrackedParameter<std::vector<std::string> >("fileNames");

  std::vector<edm::EventRange>  whichEventsToProcess = source.getUntrackedParameter<std::vector<edm::EventRange> >("eventsToProcess",std::vector<edm::EventRange>());
  std::vector<edm::LuminosityBlockRange>  whichLumisToProcess = source.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess",std::vector<edm::LuminosityBlockRange>());

  /// initialize the EventSkipperByID
  eventSkipper_ = edm::EventSkipperByID::create(source);

  ///Get number of events to process
  edm::ParameterSet maxEvents =  parameterSet_->getUntrackedParameter<edm::ParameterSet>("maxEvents");
  nEventsToAnalyze_ = maxEvents.getUntrackedParameter<int>("input");

  ///Get sample name (temporary, should be encoded in the POOL file)
  edm::ParameterSet metadata =  parameterSet_->getUntrackedParameter<edm::ParameterSet>("configurationMetadata");
  sampleName_ = metadata.getUntrackedParameter<std::string>("name");

  ///Temporary, should be encoded in the GenRunInfoProduct
  float preselectionEffTmp = metadata.getUntrackedParameter<double>("preselectionEff",-1.0);
  preselectionEff_ = (preselectionEff_ < 0.) ? 1.0 : preselectionEffTmp;
  
  filePath_ = metadata.getUntrackedParameter<std::string>("outputPath","./");
  
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

  myAnalyzers_ = myAnalyzers;
  mySummary_ = new FWLiteSummaryAnalyzer("Summary");
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

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  TreeAnalyzer::finalize(){

  xSectionWMean_ /= xSectionWeight_;
  for(unsigned int i=0;i<myAnalyzers_.size();++i) myAnalyzers_[i]->finalize();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int TreeAnalyzer::loop(){

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
     if(ev.getTFile() != currentFile || ev.id().run()!= currentRun_ || currentRun_==0){
       currentFile = ev.getTFile();
       currentRun_ = ev.id().run();
       processRunInfo(ev.getRun());
       //runInfoAccounting_->processRunInfo(ev.getRun());
     }
     analyze(ev,ev.getRun());
   }
   
   std::cout << "Events skipped: " << nEventsSkipped_ << std::endl ;
   return nEventsAnalyzed_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::processRunInfo(const edm::RunBase& aRun){

  for(unsigned int i=0;i<myAnalyzers_.size();++i){
    myAnalyzers_[i]->processRunInfo(aRun);
  } 

 //now try to find the MC weights
  edm::Handle<GenRunInfoProduct> gen;
  edm::InputTag genInfoLabel("generator");
  aRun.getByLabel(genInfoLabel,gen);
  if(!gen.isValid()){
    std::cout<<"Collection: "<<genInfoLabel<<" is missing."<<std::endl;   
  }
  else{
    xSectionWMean_ = gen->crossSection();
    genPreselectionEff_ = gen->filterEfficiency();
  }
  ///Count initial and final number of events used in reco file processing,
  ///e.g. PF2PAT production.       
  typedef MEtoEDM<double> MEtoEDMD;

  edm::Handle<MEtoEDMD> hist;
  edm::InputTag counterLabel("MEtoEDMConverter","MEtoEDMConverterRun","PAT");
  aRun.getByLabel(counterLabel,hist);
  if(!hist.isValid()){
    std::cout<<"Collection: "<<counterLabel<<" is missing."<<std::endl;   
  }
  else{    
    const MEtoEDMD::MEtoEdmObjectVector objects = hist->getMEtoEdmObject();
    for(MEtoEDMD::MEtoEdmObjectVector::const_iterator it = objects.begin(); it != objects.end(); ++it ){
      if(it->name == "initialEvents")initialRecoEvents_ += it->object;
      if(it->name == "finalEvents") finalRecoEvents_ += it->object;
    }    
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool TreeAnalyzer::analyze(const edm::EventBase& iEvent,
				 const edm::RunBase& iRun
				 ){

  clear();
  ////////
  // check that the event skipper exists, then apply it
  boost::hash<std::string> string_hash;
  char text[500];
  sprintf(text,"%d:%d:%d",iEvent.id().run(),
                          iEvent.id().luminosityBlock(),
	                  iEvent.id().event());
  std::string tmpString;
  tmpString.append(text);

  runInfoAccounting_->processRunInfo(iEvent,iRun);

  bool skipUsingList = eventsToProcessHash_.size() && !eventSkipper_.get() &&
    eventsToProcessHash_.find(tmpString)==eventsToProcessHash_.end();
  
  bool skipUsingRanges = eventSkipper_.get() && !eventsToProcessHash_.size() &&
    eventSkipper_->skipIt(iEvent.id().run(), iEvent.luminosityBlock(),iEvent.id().event());

  skipUsingRanges = false;

  if(!skipUsingList && !skipUsingRanges) {
    ///////
    for(unsigned int i=0;i<myAnalyzers_.size();++i){
      ///If analyzer returns false, skip to the last one, the Summary, unless filtering is disabled for this analyzer.
      if(!myAnalyzers_[i]->analyze(iEvent,myObjMessenger_) && myAnalyzers_[i]->filter() && myAnalyzers_.size()>1) i = myAnalyzers_.size()-2;
    }
    ///Clear all the analyzers, even if it was not called in this event.
    ///Important for proper TTree filling.
    for(unsigned int i=0;i<myAnalyzers_.size();++i) myAnalyzers_[i]->clear(); 
        
    myObjMessenger_->clear();
    ++nEventsAnalyzed_;
  }
  else ++nEventsSkipped_;
  

  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TreeAnalyzer::clear(){

  myStrSelections_->set(false);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
