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
  
  if(sysType.find("OMTF")!=std::string::npos){    
    return aL1Cand.type==L1Obj::OMTF_emu && aL1Cand.q>4;
  }
  else return false;
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

  //std::cout<<selectedCand<<std::endl;
  ///////////////////////////////////////////
  bool passPtCut = selectedCand.ptValue()>=ptCut;
  
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, genMuMom.Pt(), passPtCut);

  ///Fill histos for eff vs eta/phi only for events at the plateau.
  if(selType.size()==0 && genMuMom.Pt()<(ptCut + 20)) return; 
  tmpName = hName+"EtaVx"+std::to_string(ptCut); 
  myHistos_->fill2DHistogram(tmpName, genMuMom.Eta(), passPtCut);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, genMuMom.Phi(), passPtCut);

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
    
  genMuMom.SetPtEtaPhi(genObjVec[0].pt(),
		       genObjVec[0].eta(),
		       genObjVec[0].phi());

  bool isOMTFAcceptance = fabs(genMuMom.Eta())>0.83 && fabs(genMuMom.Eta())<1.24;
  if(!isOMTFAcceptance) return true;

  std::string sysType = "OMTF";
  std::string selType = "";
  for(int iCut=0;iCut<22;++iCut){
    fillTurnOnCurve(iCut, sysType, selType);
  }

  int iCut = 19;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = OMTFHistograms::ptBins[19];
    if(iType==0) pass = genMuMom.Pt()>24;
    if(iType==1) pass = genMuMom.Pt()>ptCut && genMuMom.Pt()<(ptCut+5);
    if(iType==2) pass = genMuMom.Pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(iCut, sysType, selType);
  }
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
