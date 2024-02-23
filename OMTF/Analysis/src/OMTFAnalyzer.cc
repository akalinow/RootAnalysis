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
double OMTFAnalyzer::calibratedPt(const std::string & sysType, const L1Obj & aCand){

  double value = aCand.ptValue();
  if(sysType=="OMTFDispU") value = aCand.ptUnconstrainedValue();  
  else if(sysType=="OMTFDisp") value = std::max(aCand.ptUnconstrainedValue(), aCand.ptValue());  
  else if(sysType=="NN") value = aCand.ptValue();
  
  return value;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::isPtGeq(const double & pt1, const double & pt2){

  ///Go back from floating point to integer pt scale used by GMT
  return int(pt1*2+1)>=int(pt2*2+1);
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::isInEtaAcceptance(const GenObj & aGenObj){

  double eta = std::abs(aGenObj.eta()); 
  double rVtx = sqrt(aGenObj.vx()*aGenObj.vx() + aGenObj.vy()*aGenObj.vy());
  double zVtx = std::abs(aGenObj.vz());
  double eta_min = 0.83;
  double eta_max = 1.24;

  bool decision = false;
  if(true) { //use eta extrapolation for all muons
    double rMB1 = 420.0; //cm
    double tan_theta_min = tan(2*atan(exp(-eta_max)));
    double tan_theta_max = tan(2*atan(exp(-eta_min)));
    double zMB1_min = rMB1/tan_theta_max;
    double zMB1_max = rMB1/tan_theta_min;
    double tan_theta = tan(2*atan(exp(-eta)));
    double deltaR = rMB1 - rVtx;
    double deltaZ = deltaR/tan_theta;
    double zExtrapol = std::abs(zVtx + deltaZ);
    decision = zExtrapol<zMB1_max && zExtrapol>zMB1_min;
  }
  else decision = fabs(eta)>eta_min && fabs(eta)<eta_max;

  return decision;
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
  
  if(sysType=="OMTF") qualitySelection &= (aL1Cand.type==L1Obj::OMTF_emu);
  else if(sysType=="NN")  qualitySelection &= aL1Cand.type==L1Obj::EMTF;
  else if(sysType=="LUT") qualitySelection = aL1Cand.q>=8 && aL1Cand.bx==0 && aL1Cand.type==L1Obj::BMTF;
  else if(sysType=="GMT") qualitySelection &= aL1Cand.type==L1Obj::uGMT_emu;
  else if(sysType=="GMTPhase2") qualitySelection &= aL1Cand.type==L1Obj::uGMTPhase2_emu;
  else if(sysType.find("Vx")!=std::string::npos) qualitySelection = true;

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
  std::string hName = "h2D"+sysType+selType;
    
  ///Find best matching L1 candidate
  float deltaR = 0.4, tmpR = 999;
  L1Obj selectedCand;

  for(auto & aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);
    if(!pass) continue;    
    double deltaEta = std::abs(myGenObj.eta()-aCand.etaValue());  
    deltaEta =   -myGenObj.eta()*aCand.etaValue();//takes value of +1 when eta signs do not match
    tmpR = deltaEta;
    if(tmpR<deltaR){
      deltaR = tmpR;
      selectedCand = aCand;
    }    
  }

  float candPt = calibratedPt(sysType, selectedCand);
  bool passPtCut = isPtGeq(candPt, ptCut);
  
  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  if(candPt>0){ 
    tmpName = hName+"PtRecVsPtGen";
    myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), candPt);
    tmpName = hName+"dxyVsPhiB";
    myHistos_->fill2DHistogram(tmpName, myGenObj.dxy(), selectedCand.d0Value());
    tmpName = hName+"dxyVsPhiBRefLayer0";
    if(selectedCand.refLayer==0) myHistos_->fill2DHistogram(tmpName, myGenObj.dxy(), selectedCand.d0Value());
    tmpName = hName+"dxyVsPhiBRefLayer2";
    if(selectedCand.refLayer==2) myHistos_->fill2DHistogram(tmpName, myGenObj.dxy(), selectedCand.d0Value());

  }

  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, std::abs(myGenObj.eta()), passPtCut);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.phi(), passPtCut);

  tmpName = hName+"dxy"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, std::abs(myGenObj.dxy()), passPtCut);

  tmpName = hName+"dz"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.dz(), passPtCut);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillRateHisto(const std::string & sysType,
				 const std::string & selType){

  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
  std::string hName = "h2D"+sysType+"Rate"+selType;

  L1Obj selectedCand;
  for(auto & aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType, selType);    
    if(pass && selectedCand.ptValue()<aCand.ptValue()) selectedCand = aCand;
  }

  float candPt = calibratedPt(sysType, selectedCand);
  int iPtCut = OMTFHistograms::iPtCuts.at(3);
  float ptCut = OMTFHistograms::ptBins.at(iPtCut);
  bool passPtCut = isPtGeq(candPt, ptCut);

  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(), candPt);
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(), passPtCut*std::abs(myGenObj.eta())+(!passPtCut)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),  passPtCut*myGenObj.pt()+(!passPtCut)*(-100));
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillRateHistos(){

  if(name()=="NU_RATEAnalyzer" && myGenObj.pt()>0.0) return;

  for(auto & anAlgo : OMTFHistograms::algos){    
    fillRateHisto(anAlgo,"Tot");
    fillRateHisto(anAlgo,"VsPt");
    fillRateHisto(anAlgo,"VsEta");
  }
  fillRateHisto("Vx","Tot");
  fillRateHisto("Vx","VsPt");
  fillRateHisto("Vx","VsEta");
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillHistosForGenMuon(){
  
  myHistos_->fill1DHistogram("h1DGenEtaAll", myGenObj.eta());
  myHistos_->fill1DHistogram("h1DGenDxyAll", myGenObj.dxy());
  myHistos_->fill1DHistogram("h1DGenDzAll", myGenObj.dz());

  if(!isInEtaAcceptance(myGenObj)) return;

  myHistos_->fill1DHistogram("h1DGenPt", myGenObj.pt());
  myHistos_->fill1DHistogram("h1DGenEta", myGenObj.eta());
  myHistos_->fill1DHistogram("h1DGenDxy", myGenObj.dxy());
  myHistos_->fill1DHistogram("h1DGenDz", myGenObj.dz());

  ///Generic turn on curves
  std::string selType = "";
  for(unsigned int iCut=0;iCut<OMTFHistograms::ptBins.size();++iCut){
    for(auto & anAlgo : OMTFHistograms::algos){    
      fillTurnOnCurve(iCut, anAlgo, selType);      
    }
  }
    
  ///Efficiency vs eta and phi for three selected point on
  ///a turn on curve
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    int iPtCut = OMTFHistograms::iPtCuts.at(3);
    float ptCut = OMTFHistograms::ptBins.at(iPtCut);
    
    if(iType==0) pass = myGenObj.pt()>ptCut + 20;
    else if(iType==1) pass = myGenObj.pt()>ptCut && myGenObj.pt()<(ptCut+5);
    else if(iType==2) pass = myGenObj.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    for(auto & anAlgo : OMTFHistograms::algos){    
      fillTurnOnCurve(iPtCut, anAlgo, selType);      
    }   
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
  if(genObjVec.empty() && name().find("NU_RATE")==std::string::npos) return false;
  
  // Filter out non-muon particles and particles with non-final state
for (auto & aGenObj : genObjVec) {
    if (std::abs(aGenObj.pdgId()) == 13 && std::abs(aGenObj.status()) == 1) {
        myGenObj = aGenObj;
        fillHistosForGenMuon();
    }
  }

  fillRateHistos();

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
