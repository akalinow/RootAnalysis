#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "GMTAnalyzer.h"
#include "EventProxyOMTF.h"

#include "TF1.h"

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
  
  bool lowPtVeto = false;

   if(sysType.find("uGMT_emu")!=std::string::npos){
     return aL1Cand.type==L1Obj::uGMT_emu && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("Vx")!=std::string::npos){
     return true;
   }   
   return false;
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
void GMTAnalyzer::fillTurnOnCurve(const RecoMuonObj & aRecoMuon,
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
  if(sysType=="uGMT_emu") {   
    hName = "h2DuGMT_emu"+selType;
  }
  if(sysType=="EMTF") {   
    hName = "h2DEMTF"+selType;
  }

  ///Find the best matching L1 candidate
  float deltaEta = 0.4;
  L1Obj selectedCand;
  
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    double delta = std::abs(aRecoMuon.eta()-aCand.etaValue());    
    if(delta<deltaEta){
      deltaEta = delta;
      selectedCand = aCand;      
    }    
  }

  bool passPtCut = selectedCand.ptValue()>=ptCut && selectedCand.ptValue()>0;

  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.pt(), selectedCand.ptValue());
  
  //Generic eff vs selected variable calculated for muons on plateau
  if(!selType.size() && aRecoMuon.pt()<ptCut+20) return;
  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.eta(), passPtCut);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.phi(), passPtCut);
/*
  tmpName = hName+"BetaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.beta(), passPtCut);

  tmpName = hName+"vxVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.vx(), passPtCut);

  tmpName = hName+"vzVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.vz(), passPtCut);

  tmpName = hName+"vyVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, aRecoMuon.vy(), passPtCut);

  float vertex_distance = sqrt(pow(aRecoMuon.vx(), 2) + pow(aRecoMuon.vy(), 2) + pow(aRecoMuon.vz(), 2));

  tmpName = hName+"vert_dist"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, vertex_distance, passPtCut);
*/
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////

void GMTAnalyzer::fillRateHisto(const RecoMuonObj & aRecoMuon,
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
void GMTAnalyzer::fillHistosForRecoMuon(const RecoMuonObj & aRecoMuon){   

  bool isGMTAcceptance = fabs(aRecoMuon.eta())<2.4;
  if(!isGMTAcceptance) return;

  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(aRecoMuon, iCut, "OMTF", selType);
      fillTurnOnCurve(aRecoMuon, iCut, "uGMT_emu", selType);
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
    fillTurnOnCurve(aRecoMuon, iCut, "uGMT_emu", selType);
  }
}
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////
bool GMTAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  const EventProxyOMTF & myProxy = static_cast<const EventProxyOMTF&>(iEvent);

  myEventId = myProxy.getEventId();
  mymuonColl = myProxy.getRecoMuonObjColl();
  myL1ObjColl = myProxy.getL1ObjColl();
  myL1PhaseIIObjColl = myProxy.getL1PhaseIIObjColl();

  const std::vector<RecoMuonObj> RecoMuonVec = mymuonColl->data();  
  if(RecoMuonVec.empty()) return false;

  for(auto aRecoMuonObj: RecoMuonVec){
   // if(std::abs(aGenObj.pdgId())!=13) continue;
   // if(std::abs(aGenObj.status())!=1) continue;    
    fillHistosForRecoMuon(aRecoMuonObj); 
  
    fillRateHisto(aRecoMuonObj, "Vx","Tot");
    fillRateHisto(aRecoMuonObj, "uGMT_emu","Tot");
  
    fillRateHisto(aRecoMuonObj, "Vx","VsPt");
    fillRateHisto(aRecoMuonObj, "uGMT_emu","VsPt");
   
    fillRateHisto(aRecoMuonObj, "Vx","VsEta");
    fillRateHisto(aRecoMuonObj, "uGMT_emu","VsEta");
  }
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
