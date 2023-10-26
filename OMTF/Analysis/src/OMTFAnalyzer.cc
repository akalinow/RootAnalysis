#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "OMTFAnalyzer.h"
#include "EventProxyOMTF.h"

#include "TF1.h"
                                      
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double OMTFAnalyzer::calibratedPt(const std::string & sysType, const double & ptRaw){

  return ptRaw;
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::isInEtaAcceptance(const double & eta){

  return fabs(eta)>0.83 && fabs(eta)<1.24;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::OMTFAnalyzer(const std::string & aName):Analyzer(aName){
  selectionFlavours_.push_back(aName);
}
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
void OMTFAnalyzer::finalize(){ myHistos_->finalizeHistograms();}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType,
			       const std::string & selType){
			       
  bool qualitySelection = aL1Cand.q>=12 && aL1Cand.bx==0;        
  if(sysType.find("Vx")!=std::string::npos) qualitySelection = true;
  
  return qualitySelection;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillTurnOnCurve(const int & iPtCut,
				  const std::string & sysType,
				  const std::string & selType){

  //int important for histo name construction
  int ptCut = OMTFHistograms::OMTFHistograms::ptBins.at(iPtCut);

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2DGmt"+selType;
  if(sysType=="OMTF") {   
    hName = "h2DOMTF"+selType;
  }
  if(sysType=="OTTF") {   
    hName = "h2DBMTF"+selType;
  }

  ///Find best matching L1 candidate
  float deltaR = 0.4, tmpR = 999;
  L1Obj selectedCand;

  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    double deltaEta = std::abs(myGenObj.eta()-aCand.etaValue());    
    tmpR = deltaEta;
    if(tmpR<deltaR){
      deltaR = tmpR;
      selectedCand = aCand;
    }    
  }

  float val = calibratedPt(sysType, selectedCand.ptValue());
  bool passPtCut = int(2*val+1)>=int(2*ptCut+1);

  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), val);

  //Generic eff vs eta/phi calculated for muons on plateau
  if(!selType.size() && myGenObj.pt()<ptCut+20) return;
  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.eta(), passPtCut);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.phi(), passPtCut);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillRateHisto(const std::string & sysType,
				 const std::string & selType){

  if(name()=="NU_RATEAnalyzer" && myGenObj.pt()>0.0) return;

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;
  }

  float val = calibratedPt(sysType, selectedCand.ptValue());
  bool pass = val>=21;

  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),val);
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.pt()+(!pass)*(-100));
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillHistosForGenMuon(){
  
  if(!isInEtaAcceptance(myGenObj.eta())) return;

  std::string selType = "";
  for(unsigned int iCut=0;iCut<OMTFHistograms::ptBins.size();++iCut){
      fillTurnOnCurve(iCut, "OMTF", selType);
      fillTurnOnCurve(iCut, "BMTF", selType);
      fillTurnOnCurve(iCut, "EMTF", selType);
  }
  
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    int iPtCut = 14;
    float ptCut = OMTFHistograms::ptBins.at(iPtCut);
    
    if(iType==0) pass = myGenObj.pt()>ptCut + 20;
    else if(iType==1) pass = myGenObj.pt()>ptCut && myGenObj.pt()<(ptCut+5);
    else if(iType==2) pass = myGenObj.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(iPtCut, "OMTF", selType);
    fillTurnOnCurve(iPtCut, "BMTF", selType);
    fillTurnOnCurve(iPtCut, "EMTF", selType);
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
  myGenObj = GenObj();

  const std::vector<GenObj> genObjVec = myGenObjColl->data();  
  if(genObjVec.empty()) return false;

  for(auto aGenObj: genObjVec){
    if(std::abs(aGenObj.pdgId())!=13) continue;
    if(std::abs(aGenObj.status())!=1) continue;
    myGenObj = aGenObj;
    fillHistosForGenMuon();
  }
  
  std::vector<int> ptCuts = {0, 10, 13, 15, 16, 18, 19, 20, 21, 22, 23};
  for(int iQuality=0;iQuality<4;++iQuality){
    std::string selType = std::string(TString::Format("quality%d",iQuality));
    fillRateHisto("OMTF","Tot_"+selType);
    fillRateHisto("BMTF","Tot_"+selType);
    fillRateHisto("EMTF","Tot_"+selType);

    fillRateHisto("Vx","VsEta_"+selType);
    fillRateHisto("OMTF","VsEta_"+selType);
    fillRateHisto("BMTF","VsEta_"+selType);
    fillRateHisto("EMTF","VsEta_"+selType);
  }

  fillRateHisto("Vx","Tot");
  fillRateHisto("OMTF","Tot");
  fillRateHisto("BMTF","Tot");
  fillRateHisto("EMTF","Tot");

  fillRateHisto("Vx","VsPt");
  fillRateHisto("OMTF","VsPt");
  fillRateHisto("BMTF","VsPt");
  fillRateHisto("EMTF","VsPt");

  fillRateHisto("Vx","VsEta");
  fillRateHisto("OMTF","VsEta");
  fillRateHisto("BMTF","VsEta");
  fillRateHisto("EMTF","VsEta");
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
