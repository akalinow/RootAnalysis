#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>
#include <cmath> 
#include <iostream>
#include "TCanvas.h"
#include "TreeAnalyzer.h"
#include "GMTAnalyzer.h"
#include "EventProxyOMTF.h"
#include "TLorentzVector.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include "TF1.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TMath.h"
#include "TROOT.h"

using namespace TMath;
using namespace std;
std::vector<double> ptRanges = {0.,   1.,   2.,   3.,   4.,   
                                  5.,   6.,   7.,   8.,   9., 
                                  10.,  11.,  12.,  13.,  14.,
                                  15.,  16.,  17.,  18.,  19.,
                                  20.,  21.,  22.,  23.,  24.,
                                  25.,  26.,  28.,  30.,  32.,  34.,
                                  36.,  38.,  40.,  50.,  60.,
                                  70.,  80.,  90.,  100., 200., 99999};
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
GMTAnalyzer::GMTAnalyzer(const std::string & aName):Analyzer(aName){

  selectionFlavours_.push_back(aName);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
GMTAnalyzer::~GMTAnalyzer(){
  
  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::initialize(TDirectory* aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;
  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new GMTHistograms(aDir, selectionFlavours_);

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
Analyzer* GMTAnalyzer::clone() const{

  GMTAnalyzer* clone = new GMTAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::finalize(){

  myHistos_->finalizeHistograms();

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType,
			       const std::string & selType){
  
   if(sysType.find("uGMT")!=std::string::npos){
     return aL1Cand.type==L1Obj::uGMT && aL1Cand.q>=12 && aL1Cand.bx==0 ;
   }
   else if(sysType.find("OMTF")!=std::string::npos){
    return aL1Cand.type==L1Obj::OMTF && aL1Cand.q>=12 && aL1Cand.bx==0;
   }
   return false;
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillTurnOnCurveGen(const GenObj & aGenObj,
                                  const int & iPtCut,
				                          const std::string & sysType,
				                          const std::string & selType){

  //int is important for histo name construction
  int ptCut = GMTHistograms::ptBins[iPtCut];

  /*const std::vector<L1PhaseIIObj> & myL1Coll = myL1PhaseIIObjColl->getL1PhaseIIObjs();*/
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2DGmt"+selType;
  if(sysType=="OMTF") {   
    hName = "h2DOMTF"+selType;
  }
  if(sysType=="uGMT_emu") {   
    hName = "h2DuGMT_emu"+selType;
  }
  if(sysType=="uGMT") {   
    hName = "h2DuGMT"+selType;
  }
  if(sysType=="EMTF") {   
    hName = "h2DEMTF"+selType;
  }

  ///Find the best matching L1 candidate
  float deltaEta = 0.4;
  L1Obj selectedCand; /*L1PhaseIIObj selectedCand;*/
  
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    /*std::cout<<aCand<<std::endl;*/
    double delta = std::abs(aGenObj.eta()-aCand.etaValue());    
    if(delta<deltaEta){
      deltaEta = delta;
      selectedCand = aCand;      
    }    
  }

  bool passPtCut = selectedCand.ptValue()>=ptCut && selectedCand.ptValue()>0;

  std::string tmpName = hName+"Pt_Gen"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aGenObj.pt(), passPtCut);

  tmpName = hName+"HighPt_Gen"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aGenObj.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, aGenObj.pt(), selectedCand.ptValue());
  

}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillTurnOnCurve(const MuonObj & aMuonCand,
                                  const int & iPtCut,
				                          const std::string & sysType,
				                          const std::string & selType){
  
 //int is important for histo name construction
  int ptCut = GMTHistograms::ptBins[iPtCut];
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2DGmt"+selType;
  if(sysType=="OMTF") {   
    hName = "h2DOMTF"+selType;
  }
  if(sysType=="uGMT") {   
    hName = "h2DuGMT"+selType;
  }
  if(sysType=="EMTF") {   
    hName = "h2DEMTF"+selType;
  }

  ///Find the best matching L1 candidate
  double deltaR = 0.15;
  L1Obj selectedCand;
 
  for(auto aCand: myL1Coll){  
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    
    double phiValue = aCand.phiValue();
    if(phiValue>M_PI) phiValue-=2*M_PI;
    
    double dEta = std::abs(aMuonCand.l1eta()-aCand.etaValue());
    double dPhi = std::abs(aMuonCand.l1phi()-phiValue);
    if(dPhi>2*M_PI) dPhi=-2*M_PI;
    double delta = sqrt(dEta*dEta + dPhi*dPhi);
    if(delta<deltaR && selectedCand.ptValue()<aCand.ptValue()){
      deltaR = delta;
      selectedCand = aCand;      
    }    
  }
  bool passPtCut = selectedCand.ptValue()>=ptCut && selectedCand.ptValue()>0;
   
  std::string tmpName = hName+"Pt_Reco"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), passPtCut);

  tmpName = hName+"HighPt_Reco"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtOMTF";
  myHistos_->fill2DHistogram(tmpName, aMuonCand.pt(), selectedCand.ptValue());
  
  //Generic eff vs selected variable calculated for muons on plateau
  if(!selType.size() && aMuonCand.pt()<ptCut+20) return;
  tmpName = hName+"EtauGMT"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.eta(), passPtCut);

  tmpName = hName+"PhiuGMT"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aMuonCand.phi(), passPtCut);
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillRateHisto(const MuonObj & aRecoMuon,
                                const std::string & sysType,
				                        const std::string & selType){

  //Generator level information is not available for the neutrino sample
  if(name()=="NU_RATEAnalyzer" && aRecoMuon.pt()>0.0) return;

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;
  }

  bool pass = selectedCand.ptValue()>=20;
  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),selectedCand.ptValue());
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),pass*aRecoMuon.eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,aRecoMuon.pt(),pass*aRecoMuon.pt()+(!pass)*(-100));
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillHistosForGenMuon(const GenObj & aGenObj){   

  bool isOMTFAcceptance = fabs(aGenObj.eta())>0.83 && fabs(aGenObj.eta())<1.24;
  if(!isOMTFAcceptance) return;

  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurveGen(aGenObj, iCut, "OMTF", selType);
      
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    
    if(iType==0) pass = aGenObj.pt()>ptCut + 20;
    else if(iType==1) pass = aGenObj.pt()>ptCut && aGenObj.pt()<(ptCut+5);
    else if(iType==2) pass = aGenObj.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurveGen(aGenObj, iCut, "OMTF", selType);
    
  }
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillHistosForRecoMuon(const MuonObj & aRecoMuon){     
  bool isOMTFAcceptance = fabs(aRecoMuon.eta())>0.83 && fabs(aRecoMuon.eta())<1.24;
  if(!isOMTFAcceptance) return;
  
  myHistos_->fill1DHistogram("h1DPtProbe", aRecoMuon.pt());
  myHistos_->fill1DHistogram("h1DAbsEtaProbe", std::abs(aRecoMuon.eta()));
  
  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(aRecoMuon, iCut, "OMTF", selType);
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
    
    if(iType==0) pass = aRecoMuon.pt()>ptCut + 20;
    else if(iType==1) pass = aRecoMuon.pt()>ptCut && aRecoMuon.pt()<(ptCut+5);
    else if(iType==2) pass = aRecoMuon.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(aRecoMuon, iCut, "OMTF", selType);
  }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::analyze(const EventProxyBase& iEvent){
   
   clear();
   
   const EventProxyOMTF & myProxy = static_cast<const EventProxyOMTF&>(iEvent);

  myEventId = myProxy.getEventId();
  myMuonObjColl = myProxy.getRecoMuonObjColl();
  myL1ObjColl = myProxy.getL1ObjColl();
  myL1PhaseIIObjColl = myProxy.getL1PhaseIIObjColl();
  myGenObjColl = myProxy.getGenObjColl();

  //////////////////////////////////////////////////////////////////////
   /////////For  muonColl as Reference with L1////////////////////////////////
   /////////////////////////////////////////////////////////////////////
  const std::vector<MuonObj> & myMuonColl = myMuonObjColl->getMuonObjs();
  if(myMuonColl.size()<2) return false;
  
  MuonObj aTag =  myMuonColl.at(0);
  bool tagPass = aTag.pt()>20 && aTag.matchedisohlt();
  if(!tagPass) return true;
  
  double m_muon = 0.10566;
  ROOT::Math::PtEtaPhiMVector p4Tag(aTag.pt(), aTag.eta(), aTag.phi(), m_muon);
  myHistos_->fill1DHistogram("h1DPtTag", aTag.pt());
  myHistos_->fill1DHistogram("h1DAbsEtaTag", std::abs(aTag.eta()));
  
  MuonObj aProbe;
  double m_Z = 90;
  double deltaM_Z = 20;
  double tmpDelta = 2*deltaM_Z;
  for (auto aMuonCand: myMuonColl){   
      ROOT::Math::PtEtaPhiMVector p4MuonCand(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), m_muon);
      tmpDelta = std::abs((p4MuonCand+p4Tag).mass()-m_Z);
      if(aMuonCand.tightID() && tmpDelta<deltaM_Z){
      deltaM_Z = tmpDelta;
      aProbe = aMuonCand;   
    }
  }
  if(aProbe.pt()<1) return true;
  
  ROOT::Math::PtEtaPhiMVector p4Probe(aProbe.pt(), aProbe.eta(), aProbe.phi(), m_muon);
  myHistos_->fill1DHistogram("h1DDiMuonMassTagProbe",(p4Probe+p4Tag).mass());   
  fillHistosForRecoMuon(aProbe);
  
  ///Add PU muons to analysis
  double deltaRCut = 0.6;
  for (auto aMuonCand: myMuonColl){   
      ROOT::Math::PtEtaPhiMVector p4MuonCand(aMuonCand.pt(), aMuonCand.eta(), aMuonCand.phi(), m_muon);
      if(aMuonCand.tightID() && 
         ROOT::Math::VectorUtil::DeltaR(p4MuonCand, p4Tag)>deltaRCut && 
         ROOT::Math::VectorUtil::DeltaR(p4MuonCand, p4Probe)>deltaRCut){
        fillHistosForRecoMuon(aMuonCand);
      }
    }


  //////////////////////////////////////////////////////////////////////
   /////////For Gen muon as Reference with L1////////////////////////////////
   /////////////////////////////////////////////////////////////////////

  const std::vector<GenObj> genObjVec = myGenObjColl->data();  
  if(genObjVec.empty()) return false;

  for(auto aGenObj: genObjVec){
    if(std::abs(aGenObj.pdgId())!=13) continue;
    if(std::abs(aGenObj.status())!=1) continue;   
    fillHistosForGenMuon(aGenObj); ////////////ONLY this
  }  
    
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
