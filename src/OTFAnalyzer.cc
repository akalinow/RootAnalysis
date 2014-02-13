#include <cstdlib> 

#include "TreeAnalyzer.h"
#include "OTFAnalyzer.h"
#include "EventProxyOTF.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::OTFAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::~OTFAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OTFHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::finalize(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::registerCuts(){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
    mySelections_->push_back("HLT"+selectionFlavours_[i]);
  /*
 ///Register tree variables
 for(unsigned int i=0;i<triggerItemNames_.size();++i){
   treeVariables_[triggerItemNames_[i]] = -999.0;
 }
 */
 }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::fillTurnOnCurve(float & ptCut, const std::string & sysType){

	bool qualityCut = true;

std::vector<L1Obj> & myL1Coll = theEvent->l1ObjectsGmt;
std::string hName = "h2DPtGmt";

//if(sysType=="Rpc") myL1Coll = theEvent->l1ObjectsRpc;
if(sysType=="Otf") {
	myL1Coll = theEvent->l1ObjectsOtf;
	qualityCut = myL1Coll.size() && myL1Coll[0].q!=103 &&
		         myL1Coll[0].q!=104 && myL1Coll[0].q!=105 &&
		         myL1Coll[0].q%100>3;
	hName = "h2DPtOtf";
	}

   bool pass = (myL1Coll.size() && myL1Coll[0].pt>=ptCut && qualityCut) || !myL1Coll.size();

   myHistos_->fill2DHistogram(hName,pass,theEvent->pt);
   myHistos_->fill2DHistogram(hName,pass,theEvent->etaHit);
   myHistos_->fill2DHistogram(hName,pass,theEvent->phiHit);

}


bool OTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events;

  std::string sysType = "Gmt";
  float ptCut = 20;

  sysType = "Gmt";
  fillTurnOnCurve(ptCut,sysType);
  sysType = "Otf";
  fillTurnOnCurve(ptCut,sysType);



  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  OTFAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  OTFAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


