#include <cstdlib> 
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include "TreeAnalyzer.h"
#include "OTFDiMuonAnalyzer.h"
#include "EventProxyOTF.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFDiMuonAnalyzer::OTFDiMuonAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFDiMuonAnalyzer::~OTFDiMuonAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFDiMuonAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OTFHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFDiMuonAnalyzer::finalize(){ 

  myHistos_->finalizeDiMuonHistograms(0,1.0); 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFDiMuonAnalyzer::registerCuts(){

  for(unsigned int i=0;i<selectionFlavours_.size();++i){
    mySelections_->push_back("HLT"+selectionFlavours_[i]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OTFDiMuonAnalyzer::passQuality(std::vector<L1Obj> * myL1Coll,
		 const std::string & sysType, 
		 int iCand){

  if(sysType.find("Gmt")!=std::string::npos){
    return myL1Coll->size()>iCand &&
      myL1Coll->operator[](iCand).eta>0.8 &&
      myL1Coll->operator[](iCand).eta<1.3 &&
      true;
  }

  if(sysType.find("Otf")!=std::string::npos){
  return myL1Coll->size()>iCand &&       
    ///Barrel            
    myL1Coll->operator[](iCand).bx/100!=99840 &&
    myL1Coll->operator[](iCand).bx/100!=34304 &&
    myL1Coll->operator[](iCand).bx/100!=3075 &&
    myL1Coll->operator[](iCand).bx/100!=36928 &&
    ////
    myL1Coll->operator[](iCand).bx/100!=12300 &&     
    ///Endcap
    myL1Coll->operator[](iCand).bx/100!=98816 &&
    myL1Coll->operator[](iCand).bx/100!=98944 &&

    myL1Coll->operator[](iCand).bx/100!=33408 &&
    myL1Coll->operator[](iCand).bx/100!=66688 && 
    myL1Coll->operator[](iCand).bx/100!=66176 && 
    ///        
    myL1Coll->operator[](iCand).q>2 &&
    true;
  }
  else return myL1Coll->size();
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
unsigned int OTFDiMuonAnalyzer::findMuon(const unsigned int iMix){

  std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsOtf);

  L1Obj aCand;
  L1Obj originalMuon, secondMuon;
  L1Obj firstRecoMuon, secondRecoMuon;


  for(unsigned int iCandTmp=0;iCandTmp<myL1Coll->size();++iCandTmp){
    if(myL1Coll->operator[](iCandTmp).bx==0) originalMuon = myL1Coll->operator[](iCandTmp);
    else if(myL1Coll->operator[](iCandTmp).bx%2==1 && 
	    myL1Coll->operator[](iCandTmp).bx/2==iMix) secondMuon = myL1Coll->operator[](iCandTmp);
  }

  float deltaR1 = 0.05, deltaR2 = 0.05, tmpR = 999.0;
  ///Looking for first candidate in di muon event
  if(originalMuon.pt>0){
    for(unsigned int iCandTmp=0;iCandTmp<myL1Coll->size();++iCandTmp){
      if(myL1Coll->operator[](iCandTmp).bx &&
	 myL1Coll->operator[](iCandTmp).bx%2==0 && 
	 int(myL1Coll->operator[](iCandTmp).bx-0.5)/2==iMix) aCand = myL1Coll->operator[](iCandTmp);
      
      tmpR = sqrt(pow(aCand.phi-originalMuon.phi,2) + 
		  pow(aCand.eta-originalMuon.eta,2));
      if(aCand.pt>0 && tmpR<deltaR1 &&
	 aCand.charge==originalMuon.charge){
	deltaR1 = tmpR;
	firstRecoMuon = aCand;
      }	
      ///Looking for second candidate in di muon event     
      tmpR = sqrt(pow(aCand.phi-secondMuon.phi,2) + 
		  pow(aCand.eta-secondMuon.eta,2));
      if(tmpR<deltaR2 &&
	 aCand.charge==secondMuon.charge){
	deltaR2 = tmpR;
	secondRecoMuon = aCand;
      }			  
    }    
  }


  if(originalMuon.pt>10 && secondMuon.pt>10){
    float deltaPhi = originalMuon.phi - secondMuon.phi;  
    if(deltaPhi>M_PI) deltaPhi-=2*M_PI;
    if(deltaPhi<-M_PI) deltaPhi+=2*M_PI;

    float deltaR =  sqrt(pow(originalMuon.phi-secondMuon.phi,2) + 
			 pow(originalMuon.eta-secondMuon.eta,2));
   
    bool pass = firstRecoMuon.pt>10 && secondRecoMuon.pt>10;
    int ptCut = 0;
    std::string selType = "";
    std::string hName = "h2DOtf"+selType;

    if(fabs(deltaPhi)>0.4 && !pass){
      std::cout<<originalMuon<<std::endl;
      std::cout<<secondMuon<<std::endl;
      std::cout<<"DeltaPhi: "<<deltaPhi<<" deltaR1: "<<deltaR1<<" deltaR2: "<<deltaR2<<std::endl;
      std::cout<<firstRecoMuon<<std::endl;
      std::cout<<secondRecoMuon<<std::endl;
      std::cout<<"----"<<std::endl;
    }

    std::string tmpName = hName+"DeltaPhi"+std::to_string(ptCut);
    if(originalMuon.charge*secondMuon.charge==1) myHistos_->fill2DHistogram(tmpName,deltaPhi,pass);

    tmpName = hName+"PhiHit"+std::to_string(ptCut);
    if(fabs(deltaPhi)<0.1 && originalMuon.charge*secondMuon.charge==1) myHistos_->fill2DHistogram(tmpName,firstRecoMuon.phi,pass);
    ///
    ptCut = 14;
    tmpName = hName+"DeltaPhi"+std::to_string(ptCut);
    if(originalMuon.charge*secondMuon.charge==1) myHistos_->fill2DHistogram(tmpName,deltaPhi,pass);
 
    tmpName = hName+"PhiHit"+std::to_string(ptCut);
    if(fabs(deltaPhi)<0.1 && originalMuon.charge*secondMuon.charge==1) myHistos_->fill2DHistogram(tmpName,firstRecoMuon.phi,pass);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFDiMuonAnalyzer::fillTurnOnCurve(const int & iPtCut,
					const std::string & sysType,
					const std::string & selType,
					const unsigned int iMuon){

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
    hName = "h2DOtf"+selType;
  }
  hName+="iMuon"+std::to_string(iMuon);

  float deltaEta = 999;
  float tmp = 999;
  unsigned int iCand = 0;
  for(unsigned int iCandTmp=0;iCandTmp<myL1Coll->size();++iCandTmp){
    tmp = fabs(myL1Coll->operator[](iCandTmp).eta-theEvent->etaMuon(iMuon)); 
    if(tmp<deltaEta){
      deltaEta = tmp;
      iCand = iCandTmp;
    }
  }

  bool qualityCut = passQuality(myL1Coll,sysType,iCand);
  bool pass = myL1Coll->size() && myL1Coll->operator[](iCand).pt>=ptCut && qualityCut && deltaEta<0.15;

  std::string tmpName = "h1DDeltaEta"+sysType+"iMuon"+std::to_string(iMuon)+std::to_string(ptCut);
  myHistos_->fill1DHistogram(tmpName,deltaEta,1.0);
  //std::cout<<"tmpName: "<<tmpName;
  //exit(0);

  tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->ptMuon(iMuon), pass);
  
  ///Fill histos for eff vs eta/phi only for events at the plateau.
  //if(selType.size()==0 && theEvent->ptMuon(iMuon)<(ptCut + 20)) return;
  
  tmpName = hName+"EtaHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->etaMuon(iMuon), pass);
  
  tmpName = hName+"PhiHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phiMuon(iMuon), pass);
  
  tmpName = hName+"EtaVx"+std::to_string(ptCut);

  myHistos_->fill2DHistogram(tmpName,theEvent->etaMuon(iMuon), pass);
  
  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phiMuon(iMuon),pass);
}
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
bool OTFDiMuonAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events;

  if(theEvent->eta<0.83 || theEvent->eta>1.24) return true;

  for(unsigned int iMix=0;iMix<10;++iMix) findMuon(iMix);
  return true;

  std::string selType = "";
  std::string sysTypeGmt="Gmt";
  std::string sysTypeOtf="Otf";
  std::string sysTypeRpc="Rpc";
  std::string sysTypeOther="Other";

  for(unsigned int iMuon=0;iMuon<2;++iMuon){
    for(int iCut=0;iCut<22;++iCut){
      if(iCut>0 && iCut<14) continue;
      fillTurnOnCurve(iCut,sysTypeGmt,selType,iMuon);
      fillTurnOnCurve(iCut,sysTypeOtf,selType,iMuon);
      fillTurnOnCurve(iCut,sysTypeRpc,selType,iMuon);
      fillTurnOnCurve(iCut,sysTypeOther,selType,iMuon);
    }   
    ////////////////
    int iCut = 2;
    bool pass = false;
    for(int iType=0;iType<=3;++iType){
      float ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsGmt[iCut]];
      if(iType==0) pass = theEvent->ptMuon(iMuon)>(ptCut + 20);
      if(iType==1) pass = theEvent->ptMuon(iMuon)>ptCut && theEvent->ptMuon(iMuon)<(ptCut+5);
      if(iType==2) pass = theEvent->ptMuon(iMuon)<10;
      if(!pass) continue;

      selType = std::string(TString::Format("Type%d",iType));
      fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeGmt,selType,iMuon);
      fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeRpc,selType,iMuon);
      fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeOther,selType,iMuon);
      fillTurnOnCurve(OTFHistograms::ptCutsOtf[iCut],sysTypeOtf,selType,iMuon);	  
    }
  }
  /////////////////
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  OTFDiMuonAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  OTFDiMuonAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFDiMuonAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

