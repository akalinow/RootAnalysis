#include <cstdlib> 
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include "TreeAnalyzer.h"
#include "PACAnalyzer.h"
#include "EventProxyPAC.h"


int iCandPAC = 0;
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
PACAnalyzer::PACAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
PACAnalyzer::~PACAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void PACAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new PACHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void PACAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0); 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void PACAnalyzer::registerCuts(){

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
bool PACAnalyzer::passQuality(std::vector<L1Obj> * myL1Coll,
		 const std::string & sysType, 
		 int iCand){

  if(sysType.find("Gmt")!=std::string::npos) return myL1Coll->size()>iCand;
  else if(sysType.find("Rpc")!=std::string::npos) return myL1Coll->size()>iCand;
  else return myL1Coll->size();
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void PACAnalyzer::fillTurnOnCurve(const int & iPtCut,
				  const std::string & sysType,
				  const std::string & selType){

  int ptCut = PACHistograms::ptBins[iPtCut];
  int iCand = 0;

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
  
  bool qualityCut = passQuality(myL1Coll,sysType,iCand);
  bool pass = myL1Coll->size() && myL1Coll->operator[](iCand).pt>=ptCut && qualityCut;
  
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


void PACAnalyzer::fillRateHisto(const std::string & sysType,
				const std::string & selType){

	int iCut = 2;
	int iCand = 0;
	float ptCut = PACHistograms::ptBins[PACHistograms::ptCutsGmt[iCut]];

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

	bool qualityCut = passQuality(myL1Coll,sysType,iCand);
	bool pass = myL1Coll->size() && qualityCut;
	float val = 0;
	if(pass) val = myL1Coll->operator[](iCand).pt;

        if(selType=="Tot") myHistos_->fill2DHistogram(hName,theEvent->pt,val);

	///Rate vs selected variable is plotted for given pt cut.
	pass = pass && (myL1Coll->operator[](iCand).pt>=ptCut);
	int q = 0;
	if(myL1Coll->size()) myL1Coll->operator[](iCand).q;

	if(selType=="VsEta") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*theEvent->eta+(!pass)*99);
	if(selType=="VsPt") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*theEvent->pt+(!pass)*(-100));
	if(selType=="VsQuality") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*q+(!pass)*(-10));	
}

bool PACAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();
 
  eventWeight_ = 1.0;
  /////////////////
  const EventProxyPAC & myEvent = static_cast<const EventProxyPAC&>(iEvent);
  theEvent = myEvent.events;

  //if(theEvent->eta<0.0 || theEvent->eta>2.1) return true;
   
  std::string selType = "";
  std::string sysTypeGmt="Gmt";
  std::string sysTypeRpc="Rpc";
  std::string sysTypeOther="Other";

  for(int iCut=0;iCut<22;++iCut){
	  if(iCut>0 && iCut<14) continue;
	  fillTurnOnCurve(iCut,sysTypeGmt,selType);
	  fillTurnOnCurve(iCut,sysTypeRpc,selType);
	  fillTurnOnCurve(iCut,sysTypeOther,selType);
  }


  ////////////////
  int iCut = 2;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
	  float ptCut = PACHistograms::ptBins[PACHistograms::ptCutsGmt[iCut]];
	  if(iType==0) pass = theEvent->pt>(ptCut + 20);
	  if(iType==1) pass = theEvent->pt>ptCut && theEvent->pt<(ptCut+5);
	  if(iType==2) pass = theEvent->pt<10;
	  if(!pass) continue;

	  selType = std::string(TString::Format("Type%d",iType));
	  fillTurnOnCurve(PACHistograms::ptCutsGmt[iCut],sysTypeGmt,selType);
	  fillTurnOnCurve(PACHistograms::ptCutsGmt[iCut],sysTypeRpc,selType);
	  fillTurnOnCurve(PACHistograms::ptCutsGmt[iCut],sysTypeOther,selType);
  }
  /////////////////

  fillRateHisto("Gmt","Tot");
  fillRateHisto("Gmt","VsEta");
  fillRateHisto("Gmt","VsPt");
  fillRateHisto("Gmt","VsQuality");

  fillRateHisto("Rpc","Tot");
  fillRateHisto("Rpc","VsEta");
  fillRateHisto("Rpc","VsPt");
  fillRateHisto("Rpc","VsQuality");

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  PACAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  PACAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void PACAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

