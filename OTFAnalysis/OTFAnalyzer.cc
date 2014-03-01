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
void OTFAnalyzer::fillTurnOnCurve(int & ptCut,
			                      const std::string & sysType,
			                      const std::string & selType){

	bool qualityCut = true;

std::vector<L1Obj> & myL1Coll = theEvent->l1ObjectsGmt;
std::string hName = "h2DGmt"+selType;

//if(sysType=="Rpc") myL1Coll = theEvent->l1ObjectsRpc;
//if(sysType=="Other") myL1Coll = theEvent->l1ObjectsOther;
if(sysType=="Otf") {
	myL1Coll = theEvent->l1ObjectsOtf;
	qualityCut = myL1Coll.size() && myL1Coll[0].q!=103 &&
		         myL1Coll[0].q!=104 && myL1Coll[0].q!=105 &&
		         myL1Coll[0].q%100>3;
	hName = "h2DOtf"+selType;
	}

   bool pass = myL1Coll.size() && myL1Coll[0].pt>=ptCut && qualityCut;


   std::string tmpName = hName+"Pt"+std::to_string(ptCut);
   myHistos_->fill2DHistogram(tmpName,theEvent->pt, pass);

   tmpName = hName+"EtaHit"+std::to_string(ptCut);
   myHistos_->fill2DHistogram(tmpName,theEvent->etaHit, pass);

   tmpName = hName+"PhiHit"+std::to_string(ptCut);
   myHistos_->fill2DHistogram(tmpName,theEvent->phiHit, pass);

   tmpName = hName+"EtaVx"+std::to_string(ptCut);
   myHistos_->fill2DHistogram(tmpName,theEvent->eta, pass);

   tmpName = hName+"PhiVx"+std::to_string(ptCut);
   myHistos_->fill2DHistogram(tmpName,theEvent->phi, pass);
}


bool OTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events;

  //std::string sysType = "Gmt";
  //int ptCut = 0;

  std::string selType = "";

  omp_set_num_threads(2);
#pragma omp parallel for
  for(int iCut=0;iCut<4;++iCut){

	  //std::cout<<"thread number: "<<omp_get_thread_num()<<std::endl;

	  int ptCut = OTFHistograms::ptCutsGmt[iCut];
	  std::string sysType = "Gmt";
	  fillTurnOnCurve(ptCut,sysType,selType);
	  sysType = "Otf";
	  fillTurnOnCurve(ptCut,sysType,selType);
	  ptCut = OTFHistograms::ptCutsOtf[iCut];
	  sysType = "Gmt";
	  fillTurnOnCurve(ptCut,sysType,selType);
	  sysType = "Otf";
	  fillTurnOnCurve(ptCut,sysType,selType);
  }
  ////////////////
  int iCut = 1;
  std::string selection = "";
  for(int iType=0;iType<3;++iType){
	  if(iType==0) selection = std::string(TString::Format("(pt>(%d + 20))",OTFHistograms::ptCutsGmt[iCut]).Data());
	  if(iType==1) selection = std::string(TString::Format("(pt>%d && pt<(%d+5))",OTFHistograms::ptCutsGmt[iCut],
			                               OTFHistograms::ptCutsGmt[iCut]).Data());
	  if(iType==2) selection = std::string(TString::Format("(pt<10)"));
	  selType = std::string(TString::Format("Type%d",iType));

	  int ptCut = OTFHistograms::ptCutsGmt[iCut];
	  std::string sysType = "Gmt";
	  fillTurnOnCurve(ptCut,sysType,selType);
	  ptCut = OTFHistograms::ptCutsOtf[iCut];
	  sysType = "Otf";
	  fillTurnOnCurve(ptCut,sysType,selType);
  }


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

