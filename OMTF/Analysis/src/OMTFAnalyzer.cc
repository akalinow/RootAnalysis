#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "OMTFAnalyzer.h"
#include "EventProxyOMTF.h"


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::OMTFAnalyzer(const std::string & aName):Analyzer(aName){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::~OMTFAnalyzer(){

  for(auto aItem: hitsPatterns){
    std::cout<<"hits: "<<aItem.first<<" "<<std::bitset<21>(aItem.first)<<" index: "<<aItem.second<<std::endl;
  }
  
  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::initialize(TDirectory* aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;
  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OMTFHistograms(aDir, selectionFlavours_);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* OMTFAnalyzer::clone() const{

  OMTFAnalyzer* clone = new OMTFAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::finalize(){ myHistos_->finalizeHistograms(); }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType){
  /*
hits: 12300 000000011000000001100 index: 0
hits: 1027 000000000010000000011 index: 1
hits: 3075 000000000110000000011 index: 2
hits: 8204 000000010000000001100 index: 3
hits: 4108 000000001000000001100 index: 4
hits: 14336 000000011100000000000 index: 5
hits: 14348 000000011100000001100 index: 6

hits: 11264 000000010110000000000 index: 1
hits: 6156 000000001100000001100 index: 2
hits: 2051 000000000100000000011 index: 5
hits: 10252 000000010100000001100 index: 7

hits: 14336 000000011100000000000 index: 0
hits: 16432 000000100000000110000 index: 9
hits: 13312 000000011010000000000 index: 3
hits: 15360 000000011110000000000 index: 5

hits: 13324 000000011010000001100 index: 0

hits: 30732 000000111100000001100 index: 2
hits: 28684 000000111000000001100 index: 3

hits: 9228 000000010010000001100 index: 0
hits: 30723 000000111100000000011 index: 1

*/
  if(sysType.find("OMTF")!=std::string::npos){
    return aL1Cand.type==L1Obj::OMTF_emu &&
    aL1Cand.q==12 && aL1Cand.bx==0;
    
    bool lowPtVeto =  aL1Cand.disc<220 && (aL1Cand.hits==1027 || aL1Cand.hits==3075 ||
					   aL1Cand.hits==4108 || aL1Cand.hits==8204 ||
					   aL1Cand.hits==12300 || aL1Cand.hits==14348 ||
					   aL1Cand.hits==11264 || aL1Cand.hits==6156 ||
					   aL1Cand.hits==2051 || aL1Cand.hits==10252 ||

					   aL1Cand.hits==14336 || aL1Cand.hits==16432 ||
					   aL1Cand.hits==13312 || aL1Cand.hits==15360 ||
					   aL1Cand.hits==13324 || aL1Cand.hits==30732 ||
					   aL1Cand.hits==28684					   
					   );
    lowPtVeto |= aL1Cand.disc<140;    
    return aL1Cand.type==L1Obj::OMTF_emu &&
    	   aL1Cand.q>0 && aL1Cand.bx==0 && !lowPtVeto;       
  }
  else if(sysType.find("kBMTF")!=std::string::npos){    
    return aL1Cand.type==L1Obj::BMTF && aL1Cand.q>0 && aL1Cand.bx==0;
  }
  else if(sysType.find("BMTF")!=std::string::npos){    
    return aL1Cand.type==L1Obj::EMTF && aL1Cand.q>0 && aL1Cand.bx==0;
  }
  else if(sysType.find("Vx")!=std::string::npos){
    return true;
  }

  return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillTurnOnCurve(const int & iPtCut,
				  const std::string & sysType,
				  const std::string & selType){

  //int important for histo name construction
  int ptCut = OMTFHistograms::ptBins[iPtCut];

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2DGmt"+selType;
  if(sysType=="OMTF") {   
    hName = "h2DOMTF"+selType;
  }
  if(sysType=="BMTF") {   
    hName = "h2DBMTF"+selType;
  }
  if(sysType=="kBMTF") {   
    hName = "h2DkBMTF"+selType;
  }

  ///Find best matching L1 candidate
  float deltaR = 0.2, tmpR = 999;
  L1Obj selectedCand;

  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType);
    if(!pass) continue;
    tmpR = pow(genMuMom.Eta()-aCand.etaValue(),2); //Only eta used, as phi is not propageted
    if(tmpR<deltaR){
      deltaR = tmpR;
      selectedCand = aCand;
    }
  }
  bool passPtCut = selectedCand.ptValue()>=ptCut;
      
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, genMuMom.Pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, genMuMom.Pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, genMuMom.Pt(), selectedCand.ptValue());

  if(genMuMom.Pt()<5 && selectedCand.ptValue()>=20 && sysType=="OMTF"){    
    myHistos_->fill1DHistogram("h1DLLH_Low", selectedCand.disc);
    if(hitsPatterns.find(selectedCand.hits)==hitsPatterns.end()){
      hitsPatterns[selectedCand.hits] = hitsPatterns.size();
    }
    int index = hitsPatterns[selectedCand.hits];
    myHistos_->fill1DHistogram("h1DHitsPattern_Low",index);
    myHistos_->fill1DHistogram("h1DDeltaEta_Low", selectedCand.etaValue());
    myHistos_->fill1DHistogram("h1DDeltaEta_Low_RefLayer", selectedCand.refLayer);    
  }
  if(genMuMom.Pt()>20 && selectedCand.ptValue()>=20 && sysType=="OMTF"){
    myHistos_->fill1DHistogram("h1DLLH_High", selectedCand.disc);
    if(hitsPatterns.find(selectedCand.hits)!=hitsPatterns.end()){
      int index = hitsPatterns[selectedCand.hits];
      myHistos_->fill1DHistogram("h1DHitsPattern_High",index);
    }
      myHistos_->fill1DHistogram("h1DDeltaEta_High", selectedCand.etaValue());
      myHistos_->fill1DHistogram("h1DDeltaEta_High_RefLayer", selectedCand.refLayer);    
  }
  
  ///Fill histos for eff vs eta/phi only for events at the plateau.
  if(selType.size()==0 && genMuMom.Pt()<(ptCut + 20)) return; 
  tmpName = hName+"EtaVx"+std::to_string(ptCut); 
  myHistos_->fill2DHistogram(tmpName, genMuMom.Eta(), passPtCut);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, genMuMom.Phi(), passPtCut);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillRateHisto(const std::string & sysType,
				 const std::string & selType){

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType);    
    if(pass) selectedCand = aCand;
  }

  float val = selectedCand.ptValue();
  if(selType=="Tot") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),val);

  bool pass = val>=30;
  if(selType=="VsEta") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),pass*genMuMom.Eta()+(!pass)*99);
  if(selType=="VsPt") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),pass*genMuMom.Pt()+(!pass)*(-100));
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillHistosForGenMuon(){

  //bool isOMTFAcceptance = fabs(genMuMom.Eta())>0.83 && fabs(genMuMom.Eta())<1.24;
  //if(!isOMTFAcceptance) return;
  
  bool isBMTFAcceptance = fabs(genMuMom.Eta())<0.83;
  if(!isBMTFAcceptance) return;

  fillRateHisto("Vx","Tot");
  fillRateHisto("OMTF","Tot");
  fillRateHisto("kBMTF","Tot");

  fillRateHisto("Vx","VsPt");
  fillRateHisto("OMTF","VsPt");
  fillRateHisto("kBMTF","VsPt");

  fillRateHisto("Vx","VsEta");
  fillRateHisto("OMTF","VsEta");
  fillRateHisto("kBMTF","VsEta");
      
  std::string selType = "";
  for(int iCut=0;iCut<22;++iCut){
    fillTurnOnCurve(iCut, "OMTF", selType);
    fillTurnOnCurve(iCut, "kBMTF", selType);
    fillTurnOnCurve(iCut, "BMTF", selType);
  }

  int iCut = 19;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = OMTFHistograms::ptBins[19];
    if(iType==0) pass = genMuMom.Pt()>ptCut + 20;
    if(iType==1) pass = genMuMom.Pt()>ptCut && genMuMom.Pt()<(ptCut+5);
    if(iType==2) pass = genMuMom.Pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(iCut, "OMTF", selType);
    fillTurnOnCurve(iCut, "BMTF", selType);
    fillTurnOnCurve(iCut, "kBMTF", selType);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  const EventProxyOMTF & myProxy = static_cast<const EventProxyOMTF&>(iEvent);

  myEventId = myProxy.getEventId();
  myGenObjColl = myProxy.getGenObjColl();
  myL1ObjColl = myProxy.getL1ObjColl();

  const std::vector<GenObj> genObjVec = myGenObjColl->data();  
  if(genObjVec.empty()) return false;

  for(auto aGenObj: genObjVec){
    genMuMom.SetPtEtaPhi(aGenObj.pt(),
			 aGenObj.eta(),
			 aGenObj.phi());
    fillHistosForGenMuon();
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
