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



  
  for(auto aCandidate: emulCandidates)myHistos_->fill1DHistogram("h1DiProcessorEMUL",(aCandidate.iProcessor+0.99)*(aCandidate.eta/fabs(aCandidate.eta)));  
  for(auto aCandidate: hwCandidates) myHistos_->fill1DHistogram("h1DiProcessorHW",(aCandidate.iProcessor+0.99)*(aCandidate.eta/fabs(aCandidate.eta)));

  for(auto aCandidate: hwCandidates){      
      myHistos_->fill1DHistogram("h1DPhiHWPositive",aCandidate.phi*(aCandidate.eta>0) + 99*(aCandidate.eta<0));
      myHistos_->fill1DHistogram("h1DPhiHWNegative",aCandidate.phi*(aCandidate.eta<0) + 99*(aCandidate.eta>0));
      myHistos_->fill1DHistogram("h1DEtaHW",aCandidate.eta/2.61*240);      
      myHistos_->fill1DHistogram("h1DPtHW",aCandidate.pt);
  }

  for(auto aCandidate: emulCandidates){      
      myHistos_->fill1DHistogram("h1DPhiEMULPositive",aCandidate.phi*(aCandidate.eta>0) + 99*(aCandidate.eta<0));
      myHistos_->fill1DHistogram("h1DPhiEMULNegative",aCandidate.phi*(aCandidate.eta<0) + 99*(aCandidate.eta>0));
      myHistos_->fill1DHistogram("h1DEtaEMUL",aCandidate.eta/2.61*240);
      myHistos_->fill1DHistogram("h1DPtEMUL",aCandidate.pt);
  }
  /*
  if(emulCandidates.size()!=hwCandidates.size() && hwCandidates.size()<2 && emulCandidates.size()<2){
    std::cout<<"Event "
	     <<" HW size: "<<hwCandidates.size()
	     <<" EMUL size: "<<emulCandidates.size()
	     <<std::endl;
    for(auto aCandidate: hwCandidates){
      if(fabs(aCandidate.eta/2.61*240)<115) std::cout<<"HW iEta: "<<aCandidate.eta/2.61*240<<" PHI: "<<aCandidate.phi<<" PT: "<<aCandidate.pt<<std::endl;
	}

    for(auto aCandidate: emulCandidates){
      if(fabs(aCandidate.eta/2.61*240)<115) std::cout<<"EMUL iEta: "<<aCandidate.eta/2.61*240<<" PHI: "<<aCandidate.phi<<" PT: "<<aCandidate.pt<<std::endl;
    }    
  }
  */

  
  /*
  if(emulCandidates.size() && hwCandidates.size()){
    std::cout<<"Event: "<<eventNumber<<" EMUL cands: "<<emulCandidates.size()
	     <<" hits: "<<std::bitset<17>(emulCandidates[0].hits)
      	     <<" quality: "<<emulCandidates[0].q
	     <<" iProcessor: "<<emulCandidates[0].iProcessor
	     <<" HW cands: "<<hwCandidates.size()
	     <<" hits: "<<std::bitset<17>(hwCandidates[0].hits>>1)
      	     <<" quality: "<<hwCandidates[0].q+3
      	     <<" iProcessor: "<<hwCandidates[0].iProcessor
	     <<std::endl;
  }
  */

  
  if(hwCandidates.size() && emulCandidates.size()==hwCandidates.size()){    
    for(unsigned int iCandidate=0;iCandidate<hwCandidates.size();++iCandidate){

      myHistos_->fill1DHistogram("h1DPtERHW",hwCandidates[iCandidate].pt);
      myHistos_->fill1DHistogram("h1DPtEREMUL",emulCandidates[iCandidate].pt);

      //if(emulCandidates[iCandidate].iProcessor !=
      //	 hwCandidates[iCandidate].iProcessor) continue;
      
      float delta =  hwCandidates[iCandidate].pt - emulCandidates[iCandidate].pt;
      myHistos_->fill1DHistogram("h1DDeltaPt",delta);
      delta =  (hwCandidates[iCandidate].phi - emulCandidates[iCandidate].phi)*576.0/(2*M_PI);
      
      myHistos_->fill1DHistogram("h1DDeltaPhi",delta);
      delta =  (hwCandidates[iCandidate].eta - emulCandidates[iCandidate].eta)/2.61*240;
      myHistos_->fill1DHistogram("h1DDeltaEta",delta);
      delta =  hwCandidates[iCandidate].charge - emulCandidates[iCandidate].charge;
      myHistos_->fill1DHistogram("h1DDeltaCharge",delta);
      myHistos_->fill1DHistogram("h1DQualityHW",hwCandidates[iCandidate].q);
      myHistos_->fill1DHistogram("h1DQualityEMUL",emulCandidates[iCandidate].q);
      /*
      std::cout<<"Event: "<<eventNumber
	       <<" EMUL: "
	       <<std::bitset<17>(emulCandidates[iCandidate].hits)
	       <<" HW: "
	       <<std::bitset<17>(hwCandidates[iCandidate].hits>>1)<<std::endl;
      */
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

  //eventNumber = theEvent->event;
  //runNumber = theEvent->run;

  std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsOtf);
  
  std::copy_if(myL1Coll->begin(),myL1Coll->end(), std::back_inserter(emulCandidates), [](const L1Obj & aCand){return (aCand.type==10);});
  std::copy_if(myL1Coll->begin(),myL1Coll->end(), std::back_inserter(hwCandidates), [](const L1Obj & aCand){return (aCand.type==11);});

  fillCandidates();
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

