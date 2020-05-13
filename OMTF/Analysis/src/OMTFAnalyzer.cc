#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "OMTFAnalyzer.h"
#include "EventProxyOMTF.h"

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
double rescaledPt(const double & ptRaw){

  return ptRaw;

  if(ptRaw>40) return ptRaw+30;
  if(ptRaw>60) return ptRaw+40;
  if(ptRaw<15) return 1.3*ptRaw;
  double a = -0.0042;
  double b = 0.90;
  double c = 0.15;

  //double a = -0.0048;
  //double b = 0.92;
  //double c = -0.36;
  double delta = std::pow(b,2)-4*a*(c-ptRaw);
  double x = -b+sqrt(delta);
  x/=(2*a);
  return x;

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
    TH2F *h = myHistos_->get2DHistogram("h2DOMTFRateVsQuality",true);
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

  std::map<int,bool> hitsMaskRefLayer;
  std::map<int,bool> hitsMaskRefLayer4hits;
 
  /*
  hitsMaskRefLayer[1027] = true;
  hitsMaskRefLayer[3075] = true;
  hitsMaskRefLayer[2051] = true;
  hitsMaskRefLayer[1664] = true;
  
  hitsMaskRefLayer[12300] = true;
  hitsMaskRefLayer[36940] = true;
  hitsMaskRefLayer[45132] = true;
  */

  /*
  hitsMaskRefLayer[1539] = true;
  hitsMaskRefLayer[515] = true;
  hitsMaskRefLayer[1664] = true;
  hitsMaskRefLayer[13312] = true;
  */
  //hitsMaskRefLayer[131968] = true;
  //hitsMaskRefLayer[65920] = true;
  

  
  //double phiB = 0.0;
  bool hasHits = false;
  for(auto aHit : myHits){
    if(aHit.iLayer==1) {
      //phiB = aHit.iPhi;
      hasHits = true;
    }
    if(aHit.iLayer==3){
      //deltaPhi = (aHit.iPhi-phiB)*aL1Cand.chargeValue();
      hasHits &= true;
    }
  }

  std::bitset<18> hitsWord(aL1Cand.hits);

  /*
  if(sysType=="BMTF" && aL1Cand.q==12 && hitsWord.test(0) && hitsWord.test(2)) {
    std::cout<<"aL1Cand.pt: "<<aL1Cand.pt<<" hasHits: "<<hasHits<<" "<<hitsWord<<" myHits.size(): "<<myHits.size()<<std::endl;
     for(auto aHit : myHits){
       std::cout<<"iLayer: "<<aHit.iLayer<<std::endl;
     }
  }
  */
   bool lowPtVeto = false;

   
   //lowPtVeto |= aL1Cand.refLayer!=0;
    
   //default OMTF
   //lowPtVeto |= (hitsMaskRefLayer[aL1Cand.hits]) && aL1Cand.pt>120;
   //lowPtVeto |= (!hitsWord.test(0) && aL1Cand.phi<20 && aL1Cand.refLayer==1);
   //lowPtVeto |= (hitsWord.count()<4 && !hitsMaskRefLayer4hits[aL1Cand.hits]);
   //lowPtVeto |= (hitsWord.count()<4);
   //lowPtVeto |= (aL1Cand.refLayer==1 && hitsWord.count()==3);
   //lowPtVeto |= (hitsWord.count()==3);
   //lowPtVeto |= !passRotatedVeto(aL1Cand);
   if(selType.find("quality0")!=std::string::npos) lowPtVeto = false;
   if(selType.find("quality1")!=std::string::npos) lowPtVeto = (aL1Cand.refLayer==1 && hitsWord.count()==3);
   if(selType.find("quality2")!=std::string::npos) lowPtVeto = (hitsWord.count()==3);
   if(selType.find("quality3")!=std::string::npos) lowPtVeto = (aL1Cand.refLayer==1 && hitsWord.count()<5);

   if(sysType.find("BMTF")!=std::string::npos){
     return aL1Cand.type==L1Obj::BMTF && aL1Cand.q>=12 && aL1Cand.bx==0 && !lowPtVeto;
   }
   else if(sysType.find("EMTF")!=std::string::npos){
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
  bool passPtCut = selectedCand.ptValue()>=ptCut;
  if(sysType=="BMTF" || sysType=="EMTF") {   
    passPtCut = rescaledPt( selectedCand.ptValue())>=ptCut && selectedCand.ptValue()>0;
    //passPtCut =  selectedCand.ptValue()>=rescaledPt(ptCut) && selectedCand.ptValue()>0;
  }

  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"HighPt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), passPtCut);

  tmpName = hName+"PtRecVsPtGen";
  myHistos_->fill2DHistogram(tmpName, myGenObj.pt(), selectedCand.ptValue());

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

  if(selectedCand.ptValue()>0 && myGenObj.pt()<10 && selectedCand.ptValue()>=20 && sysType=="BMTF"){
    myHistos_->fill1DHistogram("h1DLLH_Low", 10*selectedCand.disc/hitsWord.count());    
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_HitCount", (int)hitsWord.count()*10);
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_RefLayer",selectedCand.refLayer*10);    
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_RefPhi",selectedCand.phi);    
  }
  if(selectedCand.ptValue()>0 && myGenObj.pt()>20 && selectedCand.ptValue()>=20 && sysType=="BMTF"){
    myHistos_->fill1DHistogram("h1DLLH_High", 10*selectedCand.disc/hitsWord.count());   
    myHistos_->fill1DHistogram("h1DHitsPattern_High_HitCount", (int)hitsWord.count()*10);
    myHistos_->fill1DHistogram("h1DHitsPattern_High_RefLayer",selectedCand.refLayer*10);    
    myHistos_->fill1DHistogram("h1DHitsPattern_High_RefPhi",selectedCand.phi);    
  }
  
  ///Fill histos for eff vs eta/phi only for events at the plateau.
  if(selType.size()==0) return;
  if(myGenObj.pt()<ptCut || myGenObj.pt()<ptCut+50) return; 
  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.eta(), passPtCut);

  if(myGenObj.pt()<100){
  tmpName = hName+selType+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName, myGenObj.phi(), passPtCut);
  }
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

  float val = selectedCand.ptValue();
  if(sysType.find("BMTF")!=std::string::npos ||
     sysType.find("EMTF")!=std::string::npos) val = rescaledPt(selectedCand.ptValue());
  bool pass = val>=25;

  if(selType.find("Tot")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),val);
  if(selType.find("VsEta")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.eta()+(!pass)*99);
  if(selType.find("VsPt")!=std::string::npos) myHistos_->fill2DHistogram(hName,myGenObj.pt(),pass*myGenObj.pt()+(!pass)*(-100));

  if(sysType.find("BMTF")!=std::string::npos && selType=="VsQuality"){
    unsigned long hits = selectedCand.hits;
    if(quality_index_map.find(hits)==quality_index_map.end()) quality_index_map[hits] = quality_index_map.size();
    int val = quality_index_map[hits];
    if(myGenObj.pt()<5) myHistos_->fill2DHistogram(hName, myGenObj.pt(), pass*val+(!pass)*(-10));
    hName = "h2D"+sysType+"Quality20";    
    if(myGenObj.pt()>30) {
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

  int iCut = 19;
  bool pass = false;
  for(int iType=0;iType<=3;++iType){
    float ptCut = OMTFHistograms::ptBins[19];    
    if(iType==0) pass = myGenObj.pt()>ptCut + 20;
    if(iType==1) pass = myGenObj.pt()>ptCut && myGenObj.pt()<(ptCut+5);
    if(iType==2) pass = myGenObj.pt()<10;
    if(!pass) continue;
    
    selType = std::string(TString::Format("Type%d",iType));
    fillTurnOnCurve(iCut, "OMTF", selType);
    fillTurnOnCurve(iCut, "BMTF", selType);
    fillTurnOnCurve(iCut, "EMTF", selType);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::passRotatedVeto(const L1Obj & aL1Cand){

  int refHit = -999;
  int phiB = -999;
  int delta = -999;
  double x1, y1;
  for(auto aHit : myHits){
    if(aHit.iLayer==0) refHit = aHit.iPhi*aL1Cand.chargeValue();
    if(aHit.iLayer==1) phiB = aHit.iPhi*aL1Cand.chargeValue();
    delta = aHit.iPhi*aL1Cand.chargeValue() - refHit;
    std::tie(x1, y1) = getRotatedPair(phiB, delta, aHit.iLayer);
    if(refHit<-990) continue;
    if(aHit.iLayer==2 || aHit.iLayer==10 || aHit.iLayer==11 || aHit.iLayer==12 || aHit.iLayer==13){
      bool pass = x1>-20;
      return pass;
    }
  }
  return true;  
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
  if(myGenObj.pt()<10 && rescaledPt(selectedCand.ptValue())>=20) badReco = true;
  else if(myGenObj.pt()>20 && rescaledPt(selectedCand.ptValue())>=20) goodReco = true;

  std::bitset<18> hitsWord(selectedCand.hits);

  if(goodReco) return;//TEST

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
   //if(iLayer==3) phiB = deltaPhi;//TEST  
   if(iLayer==1 || iLayer==3 || iLayer==5) deltaPhi +=refHit*myGenObj.charge();
   if(iLayer==1 || iLayer==3 || iLayer==5) deltaPhi = (aHit.iPhi-phiB)*myGenObj.charge();

   myHistos_->fill3DHistogram("h3DBending", myGenObj.pt(), deltaPhi, iLayer);   
   if(iLayer==3) deltaMB2 = deltaPhi;
   if(iLayer==10) deltaRB1In = deltaPhi;
   if(iLayer==11) deltaRB1Out = deltaPhi;
   if(iLayer==12) deltaRB2In = deltaPhi;
   if(iLayer==13) deltaRB2Out = deltaPhi;
   if(iLayer==16) deltaRE23 = deltaPhi;
   std::tie(x1, y1) = getRotatedPair(phiB, deltaPhi, iLayer);
   myHistos_->fill3DHistogram("h3DBendingRotated", myGenObj.pt(), x1, iLayer);
   myHistos_->fill3DHistogram("h3DBendingRotated", myGenObj.pt(), y1, iLayer+4);

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
std::pair<double, double> OMTFAnalyzer::getRotatedPair(double phiB, double deltaPhi, int iLayer){

  int index = -1;
  if(iLayer==3) index = 0;
  /*
  else if(iLayer==11) index = 1;
  else if(iLayer==2) index = 2;
  else if(iLayer==12) index = 3;
  else if(iLayer==13) index = 4;
  else if(iLayer==16) index = 5;
  */
  else return  std::pair<double, double>(-999, -999);
  
  std::vector<double> p0 = {-0.478105, -0.534923, -0.0552107, 0.0916314, -0.562927, -0.641144};
  std::vector<double> p1 = {-0.0689445, 0.0577801, 0.217325, 0.189808, 0.246022, 0.297808};

  p0[0]= {-2.8734};
  p1[0] = {0.460671};
  
  double angle = atan(p1.at(index));  
  double x = phiB;
  double y = (deltaPhi-p0.at(index));
  double x1 = x*cos(angle) + y*sin(angle);
  double y1 = -x*sin(angle) + y*cos(angle);

  return std::pair<double, double>(x1, y1);
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

  const std::vector<GenObj> genObjVec = myGenObjColl->data();  
  if(genObjVec.empty()) return false;

  for(auto aGenObj: genObjVec){
    if(std::abs(aGenObj.pdgId())!=13) continue;
    if(std::abs(aGenObj.status())!=1) continue;
    myGenObj = aGenObj;
    fillHistosForGenMuon();
    fillBendingHistos("OMTF");
  }

  std::vector<int> ptCuts = {10, 15, 16, 18, 19, 20, 21, 22, 23};
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
