#include <cstdlib> 
#include <string>
#include <omp.h>

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
void OTFAnalyzer::finalize(){ myHistos_->finalizeHistograms(0,1.0); }
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
void OTFAnalyzer::fillTurnOnCurve(const int & iPtCut,
			                      const std::string & sysType,
			                      const std::string & selType){

  bool qualityCut = true;
  int ptCut = OTFHistograms::ptBins[iPtCut];
  
  std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsGmt);
  std::string hName = "h2DGmt"+selType;
  
  if(sysType=="Rpc"){
    myL1Coll = &(theEvent->l1ObjectsRpc);
    hName = "h2DRpc"+selType;
  }
  if(sysType=="Other"){
    myL1Coll = &(theEvent->l1ObjectsOther);
    hName = "h2DOther"+selType;
  }
  if(sysType=="Otf") {
    myL1Coll = &(theEvent->l1ObjectsOtf);
    qualityCut = myL1Coll->size() &&
      myL1Coll->operator[](0).q!=103 &&
      myL1Coll->operator[](0).q!=203 &&
      myL1Coll->operator[](0).q!=303 &&
      myL1Coll->operator[](0).q%100>3;
    hName = "h2DOtf"+selType;
  }
  
  bool pass = myL1Coll->size() && myL1Coll->operator[](0).pt>=ptCut && qualityCut;
  
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->pt, pass);
  
  ///Fill histos for eff vs eta/phi only for events at the plateau.
  if(selType.size()==0 && theEvent->pt<(ptCut + 20)) return;
  
  tmpName = hName+"EtaHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->etaHit, pass);
  
  tmpName = hName+"PhiHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phiHit, pass);
  
  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->eta, pass);
  
  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phi, pass);
}


void OTFAnalyzer::fillRateHisto(const std::string & sysType,
			                    const std::string & selType){

	bool qualityCut = true;
	int iCut = 2;
	float ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsGmt[iCut]];

	std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsGmt);
	std::string hName = "h2DRate"+selType+"Gmt";

	if(sysType=="Rpc"){
			myL1Coll = &(theEvent->l1ObjectsRpc);
			hName = "h2DRate"+selType+"Rpc";
	}
	if(sysType=="Other"){
			myL1Coll = &(theEvent->l1ObjectsOther);
			hName = "h2DRate"+selType+"Other";
	}
	if(sysType=="Otf") {
	  myL1Coll = &(theEvent->l1ObjectsOtf);
	  qualityCut = myL1Coll->size() && 
	    myL1Coll->operator[](0).q!=103 &&
	    myL1Coll->operator[](0).q!=203 && 
	    myL1Coll->operator[](0).q!=303 &&
	    myL1Coll->operator[](0).q%100>3;
	  hName = "h2DRate"+selType+"Otf";
	  ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsOtf[iCut]];

	}

	bool pass = myL1Coll->size() && qualityCut;
	float val = 0;
	if(pass) val = myL1Coll->operator[](0).pt;
        if(selType=="Tot") myHistos_->fill2DHistogram(hName,theEvent->pt,val);
	pass = myL1Coll->size() && qualityCut && (myL1Coll->operator[](0).pt>=ptCut);
	if(selType=="VsEta") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*theEvent->eta+(!pass)*99);
	

}

bool OTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events;

  if( !(theEvent->eta<1.24 && theEvent->eta>0.8) ) return true;
  //if(theEvent->eta<1.9) return true;
  //if(theEvent->eta>0.8) return true;

  std::string selType = "";
  std::string sysTypeGmt="Gmt";
  std::string sysTypeOtf="Otf";
  std::string sysTypeRpc="Rpc";
  std::string sysTypeOther="Other";

  
  //omp_set_num_threads(2);
  //pragma omp parallel for
  for(int iCut=0;iCut<22;++iCut){
	  if(iCut>0 && iCut<14) continue;
	  //std::cout<<"thread number: "<<omp_get_thread_num()<<std::endl;
	  fillTurnOnCurve(iCut,sysTypeGmt,selType);
	  fillTurnOnCurve(iCut,sysTypeOtf,selType);
	  fillTurnOnCurve(iCut,sysTypeRpc,selType);
	  fillTurnOnCurve(iCut,sysTypeOther,selType);
  }
  ////////////////
  int iCut = 2;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
	  float ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsGmt[iCut]];
	  if(iType==0) pass = theEvent->pt>(ptCut + 20);
	  if(iType==1) pass = theEvent->pt>ptCut && theEvent->pt<(ptCut+5);
	  if(iType==2) pass = theEvent->pt<10;
	  if(!pass) continue;

	  selType = std::string(TString::Format("Type%d",iType));
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeGmt,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeRpc,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeOther,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsOtf[iCut],sysTypeOtf,selType);

  }
  /////////////////
  fillRateHisto("Gmt","Tot");
  fillRateHisto("Otf","Tot");

  fillRateHisto("Gmt","VsEta");
  fillRateHisto("Otf","VsEta");


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

