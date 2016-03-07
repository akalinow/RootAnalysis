#include <cstdlib> 
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>
#include <algorithm>
#include <iostream>

#include "TreeAnalyzer.h"
#include "HWvsEMULAnalyzer.h"
#include "EventProxyOTF.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HWvsEMULAnalyzer::HWvsEMULAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HWvsEMULAnalyzer::~HWvsEMULAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HWvsEMULAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new HWvsEMULHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HWvsEMULAnalyzer::clone() const{

  HWvsEMULAnalyzer* clone = new HWvsEMULAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HWvsEMULAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0); 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  HWvsEMULAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HWvsEMULAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;

  emulCandidates.clear();
  hwCandidates.clear();
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HWvsEMULAnalyzer::fillCandidates(){


  for(auto aCandidate: emulCandidates) myHistos_->fill1DHistogram("h1DiProcessorEMUL",aCandidate.iProcessor*(aCandidate.eta/fabs(aCandidate.eta)));    
  for(auto aCandidate: hwCandidates) myHistos_->fill1DHistogram("h1DiProcessorHW",aCandidate.iProcessor*(aCandidate.eta/fabs(aCandidate.eta)));

  for(auto aCandidate: hwCandidates){
    myHistos_->fill1DHistogram("h1DPhiHWPositive",aCandidate.phi*(aCandidate.eta>0));
    myHistos_->fill1DHistogram("h1DPhiHWNegative",aCandidate.phi*(aCandidate.eta<0));
    myHistos_->fill1DHistogram("h1DEtaHW",aCandidate.eta);
    myHistos_->fill1DHistogram("h1DPtHW",aCandidate.pt);
  }

  if(hwCandidates.size() && emulCandidates.size()==hwCandidates.size()){    
    for(unsigned int iCandidate=0;iCandidate<hwCandidates.size();++iCandidate){

      if(std::bitset<17>(emulCandidates[iCandidate].hits)!=
      	 std::bitset<17>(hwCandidates[iCandidate].hits>>1)) continue;
      
      float delta =  hwCandidates[iCandidate].pt - emulCandidates[iCandidate].pt;
      myHistos_->fill1DHistogram("h1DDeltaPt",delta);
      delta =  hwCandidates[iCandidate].phi - emulCandidates[iCandidate].phi;
      myHistos_->fill1DHistogram("h1DDeltaPhi",delta);
      delta =  hwCandidates[iCandidate].eta - emulCandidates[iCandidate].eta;
      myHistos_->fill1DHistogram("h1DDeltaEta",delta);
      delta =  hwCandidates[iCandidate].charge - emulCandidates[iCandidate].charge;
      myHistos_->fill1DHistogram("h1DDeltaCharge",delta);
      myHistos_->fill1DHistogram("h1DQualityHW",hwCandidates[iCandidate].q+3);
      myHistos_->fill1DHistogram("h1DQualityEMUL",emulCandidates[iCandidate].q);
      std::cout<<"emulCandidates[iCandidate].q: "
	       <<std::bitset<17>(emulCandidates[iCandidate].hits)<<" "
	       <<std::bitset<17>(hwCandidates[iCandidate].hits>>1)<<std::endl;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HWvsEMULAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.event;

  std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsOtf);
  
  std::copy_if(myL1Coll->begin(),myL1Coll->end(), std::back_inserter(emulCandidates), [](const L1Obj & aCand){return (aCand.type==10);});
  std::copy_if(myL1Coll->begin(),myL1Coll->end(), std::back_inserter(hwCandidates), [](const L1Obj & aCand){return (aCand.type==11);});

  fillCandidates();
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

