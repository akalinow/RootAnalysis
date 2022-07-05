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

std::vector<double> ptRanges = {0.,   1.,   2.,   3.,   4.,   
                                  5.,   6.,   7.,   8.,   9., 
                                  10.,  11.,  12.,  13.,  14.,
                                  15.,  16.,  17.,  18.,  19.,
                                  20.,  21.,  22.,  23.,  24.,
                                  25.,  26.,  28.,  30.,  32.,  34.,
                                  36.,  38.,  40.,  50.,  60.,
                                  70.,  80.,  90.,  100., 200., 99999};
/*
 std::vector<double> ptRanges ={0., 0.1,
				 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
				 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
				 50., 60., 70., 80., 90., 100., 120., 140.,
				 160. };
*/

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double OMTFAnalyzer::calibratedPt(const std::string & sysType, const double & ptRaw){

  if(sysType!="BMTF") return ptRaw;

  if(!calibrationFunc){
    //calibrationFile = new TFile("PtCalibration_v34.root");
    //hCalibationHisto = (TH1F*)calibrationFile->Get("ptCalibrationPSumInter035_ptToPtCalib");
    calibrationFunc = new TF1("f1","-1.042+1.037*x-0.0063*x*x", 0,100);
  }
  else if(calibrationFunc && ptRaw>0){    
    return 2*ptRaw - calibrationFunc->Eval(ptRaw);
    int iBin = hCalibationHisto->FindBin(ptRaw);
    return hCalibationHisto->GetBinContent(iBin);
  }
  return ptRaw;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::OMTFAnalyzer(const std::string & aName):Analyzer(aName){

  hGoldenPatterns = 0;
  hPtProfile = 0;
  selectionFlavours_.push_back(aName);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::~OMTFAnalyzer(){
  
  delete myHistos_;
  delete hGoldenPatterns;
  delete hPtProfile;
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
void OMTFAnalyzer::finalize(){


  fixQualityHistos();
  myHistos_->finalizeHistograms();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fixQualityHistos(){
  
  if(omp_get_max_threads()==1){
    std::cout<<"quality_index_map.size(): "<<quality_index_map.size()<<std::endl;

    std::ostringstream stringStr;
    TH2F *h = myHistos_->get2DHistogram("h2DEMTFRateVsQuality",true);
    if(h){
      for(int iBin=1;iBin<h->GetYaxis()->GetNbins();++iBin){
	stringStr.str("");
	stringStr<<iBin;
	h->GetYaxis()->SetBinLabel(iBin,stringStr.str().c_str());
      }

      for(auto it: quality_index_map){
	int iBinX = h->GetYaxis()->FindBin(it.second);
	if(iBinX>=h->GetYaxis()->GetNbins()) continue;
	if(iBinX==0){
	  std::cout<<" it.second: "<<it.second<<std::endl;
	}

	std::bitset<21> bits(it.first);

	stringStr.str("");
	stringStr<<it.first;
	std::string label = bits.to_string()+" "+stringStr.str();
	std::cout<<"iBinX: "<<iBinX<<" label: "<<label<<std::endl;
	h->GetYaxis()->SetBinLabel(iBinX,label.c_str());
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::passQuality(const L1Obj & aL1Cand,
			       const std::string & sysType,
			       const std::string & selType){
  
  bool lowPtVeto = false;

   if(sysType.find("BMTF")!=std::string::npos){
     return aL1Cand.type==L1Obj::BMTF && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("EMTF")!=std::string::npos){
     lowPtVeto = aL1Cand.disc/10>6;
     return aL1Cand.type==L1Obj::EMTF && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("OMTF")!=std::string::npos){    
     return aL1Cand.type==L1Obj::OMTF_emu && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;    
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
  if(sysType=="EMTF") {   
    hName = "h2DEMTF"+selType;
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

  std::bitset<18> hitsWord(selectedCand.hits);
  float val = calibratedPt(sysType, selectedCand.ptValue());
  bool passPtCut = val>=ptCut;

  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), val);

  if(false && iPtCut==0 && sysType=="EMTF" &&  myGenObj.pt()<10 && selectedCand.ptValue()>=20){
    std::cout<<myGenObj<<std::endl;
    std::cout<<selectedCand<<" hits: "<<selectedCand.hits<<std::endl;
    for(auto aCand: myL1Coll){
      bool pass = passQuality(aCand , "OMTF", selType);
      if(pass) std::cout<<aCand<<std::endl;
    }
    for(auto aHit : myHits){
      std::cout<<"iLayer: "<<aHit.iLayer<<" iPhi: "<<aHit.iPhi<<std::endl;
    };
    
    std::cout<<"-------"<<std::endl;
  }

  if(selectedCand.ptValue()>0 && myGenObj.pt()<10 && selectedCand.ptValue()>=20 && sysType=="EMTF"){
    myHistos_->fill1DHistogram("h1DLLH_Low", selectedCand.disc/10.0);    
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_HitCount", (int)hitsWord.count()*10);
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_RefLayer",selectedCand.refLayer*10);    
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_RefPhi",selectedCand.phi);    
  }
  if(selectedCand.ptValue()>0 && myGenObj.pt()>20 && myGenObj.pt()<30
     && selectedCand.ptValue()>=20 && sysType=="EMTF"){
    myHistos_->fill1DHistogram("h1DLLH_High", selectedCand.disc/10.0);   
    myHistos_->fill1DHistogram("h1DHitsPattern_High_HitCount", (int)hitsWord.count()*10);
    myHistos_->fill1DHistogram("h1DHitsPattern_High_RefLayer",selectedCand.refLayer*10);    
    myHistos_->fill1DHistogram("h1DHitsPattern_High_RefPhi",selectedCand.phi);    
  }

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
  bool pass = val>=20;

  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),val);
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.pt()+(!pass)*(-100));

  if(sysType.find("EMTF")!=std::string::npos && selType=="VsQuality"){
    unsigned long hits = selectedCand.hits;
    if(quality_index_map.find(hits)==quality_index_map.end()) quality_index_map[hits] = quality_index_map.size();
    int val = quality_index_map[hits];
    if(myGenObj.pt()<10) myHistos_->fill2DHistogram(hName, myGenObj.pt(), pass*val+(!pass)*(-10));
    hName = "h2D"+sysType+"Quality20";    
    if(myGenObj.pt()>20) {
      myHistos_->fill2DHistogram(hName, val, pass);
    }
  }
  else if(selType=="VsQuality"){
     myHistos_->fill2DHistogram(hName, myGenObj.pt(), 0);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillHistosForGenMuon(){

  bool isOMTFAcceptance = fabs(myGenObj.eta())>0.83 && fabs(myGenObj.eta())<1.24;
  if(!isOMTFAcceptance) return;

  std::string selType = "";
  for(int iCut=0;iCut<31;++iCut){
      fillTurnOnCurve(iCut, "OMTF", selType);
      fillTurnOnCurve(iCut, "BMTF", selType);
      fillTurnOnCurve(iCut, "EMTF", selType);
  }

  int iCut = 18;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = OMTFHistograms::ptBins[iCut];
    
    if(iType==0) pass = myGenObj.pt()>ptCut + 20;
    else if(iType==1) pass = myGenObj.pt()>ptCut && myGenObj.pt()<(ptCut+5);
    else if(iType==2) pass = myGenObj.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(iCut, "OMTF", selType);
    fillTurnOnCurve(iCut, "BMTF", selType);
    fillTurnOnCurve(iCut, "EMTF", selType);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillBendingHistos(const std::string & sysType){

   ///Find best matching L1 candidate
  float deltaR = 0.4, tmpR = 999;  
  L1Obj selectedCand;
  const std::vector<L1Obj> & myL1Coll = myL1ObjColl->getL1Objs();
 
  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType);
    if(!pass) continue;    
    double deltaEta = std::abs(myGenObj.eta()-aCand.etaValue());    
    tmpR = deltaEta;
    if(tmpR<deltaR){
      deltaR = tmpR;
      selectedCand = aCand;
    }    
  }
  
  int refHit = 9999;
  for(auto aHit : myHits){
    if(aHit.iLayer==0) refHit = aHit.iPhi;
  };
  if(refHit>999) return;

  bool badReco = false;
  bool goodReco = false;
  if(myGenObj.pt()<10 && calibratedPt(sysType, selectedCand.ptValue())>=20) badReco = true;
  else if(myGenObj.pt()>20 && calibratedPt(sysType, selectedCand.ptValue())>=20) goodReco = true;

  std::bitset<18> hitsWord(selectedCand.hits);

  double phiB = 9999.0;
  double deltaMB2 = 9999.0;
  double deltaRB1In = 9999.0;
  double deltaRB1Out = 9999.0;
  double deltaRB2In = 9999.0;
  double deltaRB2Out = 9999.0;
  double deltaRE23 = 9999.0;
  double x1, y1;

 for(auto aHit : myHits){
   int deltaPhi = (aHit.iPhi-refHit)*myGenObj.charge();
   int iLayer = aHit.iLayer;
   if(iLayer==1) phiB = aHit.iPhi;
   if(iLayer==3) phiB = deltaPhi;
   if(iLayer==1 || iLayer==3 || iLayer==5) deltaPhi +=refHit*myGenObj.charge();
   if(iLayer==1 || iLayer==3 || iLayer==5) deltaPhi = (aHit.iPhi-phiB)*myGenObj.charge();

   myHistos_->fill3DHistogram("h3DBending", myGenObj.pt(), deltaPhi, iLayer);   
   if(iLayer==3) deltaMB2 = deltaPhi;
   if(iLayer==10) deltaRB1In = deltaPhi;
   if(iLayer==11) deltaRB1Out = deltaPhi;
   if(iLayer==12) deltaRB2In = deltaPhi;
   if(iLayer==13) deltaRB2Out = deltaPhi;
   if(iLayer==16) deltaRE23 = deltaPhi;

   x1 = phiB;
   y1 = deltaMB2;
   if(iLayer==3 && goodReco) myHistos_->fill2DHistogram("h2DDeltaVsDelta_xy", x1, y1);
   if(iLayer==3 && badReco) myHistos_->fill2DHistogram("h2DDeltaVsDelta_xy_badReco", x1, y1);
   
  }

 myHistos_->fill2DHistogram("h2DDeltaVsDelta_MB2",phiB, deltaMB2);
 myHistos_->fill2DHistogram("h2DDeltaVsDelta_RB1In",phiB, deltaRB1In);
 myHistos_->fill2DHistogram("h2DDeltaVsDelta_RB1Out",phiB, deltaRB1Out);  
 myHistos_->fill2DHistogram("h2DDeltaVsDelta_RB2In",phiB, deltaRB2In);
 myHistos_->fill2DHistogram("h2DDeltaVsDelta_RB2Out",phiB, deltaRB2Out);
 myHistos_->fill2DHistogram("h2DDeltaVsDelta_RE23",phiB, deltaRE23);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<double, double> OMTFAnalyzer::getPtProfile(){

  #pragma omp critical(GOLDENPATTERNS_INITIALIZATION)
  {
    if(!hGoldenPatterns){
      
      TFile f("GPs.root");
      hGoldenPatterns = (TH3F*)f.Get("EMULAnalyzer/h3DBending")->Clone("hGoldenPatterns");
      hPtProfile = (TH1D*)hGoldenPatterns->ProjectionX();
      
      hGoldenPatterns->SetDirectory(0);
      hPtProfile->Reset();
      hPtProfile->SetDirectory(0);
    }
  }

  std::vector<int> * myHits = 0;//TEST
  if(!myHits || !myHits->size()) return std::pair<double, double>(0,0);

  double meanHit = 0;
  int layerCount = 0;
  for(int iLayer=0;iLayer<18;++iLayer){
    int phi = myHits->at(iLayer);
    if(phi==5400 || iLayer==1 || iLayer==3 || iLayer==5) continue; 
    meanHit += phi;
    layerCount++;
  }
  meanHit/=layerCount;
  meanHit = (int)meanHit;
  if(layerCount<4) return std::pair<double, double>(0,0);
  
  int iBinY, iBinZ;
  double score = 0.0;  
  hPtProfile->Reset();  
  for(int iBinX=1;iBinX<=hGoldenPatterns->GetNbinsX();++iBinX){
    score = 0.0;
    layerCount = 0;
    for(int iLayer=0;iLayer<18;++iLayer){
      int phi = myHits->at(iLayer);
      if(phi==5400 || iLayer==1 || iLayer==3 || iLayer==5) continue;
      int deltaPhi = phi-meanHit;
      iBinZ = hGoldenPatterns->GetZaxis()->FindBin(iLayer);
      iBinY = hGoldenPatterns->GetYaxis()->FindBin(deltaPhi);            
      score += hGoldenPatterns->GetBinContent(iBinX, iBinY, iBinZ);      
    }
    hPtProfile->SetBinContent(iBinX, score);
  }

  double maximumBinCenter = hPtProfile->GetBinCenter(hPtProfile->GetMaximumBin());
  double maximum = hPtProfile->GetMaximum();
  double sigma = hPtProfile->GetStdDev();
  //double skewness = hPtProfile->GetSkewness();
  maximum = sigma;
  return std::pair<double, double>(maximum, maximumBinCenter);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  const EventProxyOMTF & myProxy = static_cast<const EventProxyOMTF&>(iEvent);

  myEventId = myProxy.getEventId();
  myGenObjColl = myProxy.getGenObjColl();
  myL1ObjColl = myProxy.getL1ObjColl();
  myHits = myProxy.getHits();
  myGenObj = GenObj();

  const std::vector<GenObj> genObjVec = myGenObjColl->data();  
  if(genObjVec.empty()) return false;

  for(auto aGenObj: genObjVec){
    if(std::abs(aGenObj.pdgId())!=13) continue;
    if(std::abs(aGenObj.status())!=1) continue;
    myGenObj = aGenObj;
    fillHistosForGenMuon();
    //fillBendingHistos("OMTF");
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

    for(auto iCut: ptCuts){
      bool isOMTFAcceptance = fabs(myGenObj.eta())>0.83 && fabs(myGenObj.eta())<1.24;
      if(!isOMTFAcceptance) continue;
      fillTurnOnCurve(iCut, "OMTF", selType);
      fillTurnOnCurve(iCut, "BMTF", selType);
      fillTurnOnCurve(iCut, "EMTF", selType);
    } 
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

  fillRateHisto("Vx","VsQuality");
  fillRateHisto("OMTF","VsQuality");
  fillRateHisto("BMTF","VsQuality");
  fillRateHisto("EMTF","VsQuality");
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
