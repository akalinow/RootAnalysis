#include <cstdlib> 

#include "PFAnalyses/CommonTools/interface/FWLiteTreeAnalyzer.h"

#include "PFAnalyses/VBFHTauTau/interface/FWLiteTriggerAnalyzer.h"
#include "PFAnalyses/VBFHTauTau/interface/JetVetoHistograms.h"

#include "FWCore/Utilities/interface/Algorithms.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include "DataFormats/FWLite/interface/Event.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
FWLiteTriggerAnalyzer::FWLiteTriggerAnalyzer(const std::string & aName)
:FWLiteAnalyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
FWLiteTriggerAnalyzer::~FWLiteTriggerAnalyzer(){

  //delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void FWLiteTriggerAnalyzer::initialize(const edm::ParameterSet& ps, 
				       const TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  triggerEventLabel_ = ps.getParameter<edm::InputTag>("triggerEventLabel");
  triggerItemNames_ = ps.getParameter<std::vector<std::string> >("triggerItemNames");
  selectionFlavours_ = ps.getUntrackedParameter<std::vector<std::string> >("selectionFlavours");

  mySelections_ = aSelections;

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void FWLiteTriggerAnalyzer::finalize(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void FWLiteTriggerAnalyzer::registerCuts(){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
    mySelections_->push_back("HLT"+selectionFlavours_[i]);
  
 ///Register tree variables
 for(unsigned int i=0;i<triggerItemNames_.size();++i){
   treeVariables_[triggerItemNames_[i]] = -999.0;
 }
 }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool FWLiteTriggerAnalyzer::analyze(const edm::EventBase& iEvent){

  clear();

  using namespace reco;

  eventWeight_ = 1.0;
  /////////////////
  ///Trigger event
  try{
    iEvent.getByLabel(triggerEventLabel_,triggerEvent_);
  }
  catch(...){
    std::cout<<"Trigger event collection label: "<<triggerEventLabel_<<" not found!";
    return false;
  };

  bool decision = false;
  for(unsigned int i=0;i<selectionFlavours_.size();++i){
    decision |= checkSelections(selectionFlavours_[i]); 
  }

  return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  FWLiteTriggerAnalyzer::checkSelections(const std::string & type){

  bool decision = false;
  unsigned int index = 0;
  /*
  ////  
  const pat::TriggerPathCollection *paths = triggerEvent_->paths();
  for(unsigned int i=0;i<paths->size();++i){
    const pat::TriggerPath triggerPath = (*paths)[i];
    if(triggerPath.name().find("Tau")==std::string::npos &&
       triggerPath.name().find("Mu")==std::string::npos) continue;
    std::cout<<"Trigger path name: "<<triggerPath.name()
	     <<" wasAccept? "<<triggerPath.wasAccept()
	     <<" prescale: "<<triggerPath.prescale()
	     <<std::endl;
  }
  ////
*/

  for(unsigned int i=0;i<triggerItemNames_.size();++i){
    const pat::TriggerPath *triggerPath = triggerEvent_->path(triggerItemNames_[i]);
    //std::cout<<"Item: "<<triggerItemNames_[i]<<" path: "<<triggerPath<<std::endl;
    //if(!triggerPath) std::cout<<"trigger path not found!"<<triggerItemNames_[i]<<std::endl;
    //if(triggerPath) std::cout<<"prescale= "<< triggerPath->prescale()<<" accept? "<<triggerPath->wasAccept()<<std::endl;
    ///For QCD stream we take only non isolated triggers now.
    ///Later wqe we take als oprescaled triggers
    if(type.find("QCD")!=std::string::npos && 
       triggerItemNames_[i].find("Iso")!=std::string::npos) continue;
    if(triggerPath && triggerPath->prescale()==1)decision|=triggerPath->wasAccept();    
    if(triggerPath) treeVariables_[triggerItemNames_[i]] = triggerPath->wasAccept()*triggerPath->prescale(); 
  }

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  FWLiteTriggerAnalyzer::addBranch(TTree *tree){ 

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void FWLiteTriggerAnalyzer::clear(){ 

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


