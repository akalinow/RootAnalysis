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

  std::map<int,bool> hitsMask; 

 hitsMask[1027] = true;
 hitsMask[3075] = true;
 hitsMask[4108] = true;
 hitsMask[8204] = true;
 hitsMask[12300] = true;
 hitsMask[14348] = true;
 hitsMask[11264] = true;
 hitsMask[6156] = true;
 hitsMask[2051] = true;
 hitsMask[10252] = true;
 hitsMask[14336] = true;
 hitsMask[16432] = true;
 hitsMask[13312] = true;
 hitsMask[15360] = true;
 hitsMask[13324] = true;
 hitsMask[30732] = true;
 hitsMask[28684] = true;
 
  /*
  iBin: 1 Quality: 000000000111000001100 3596 rate: 645.541 efficiency: 0.000212619
  iBin: 2 Quality: 000000000100000001100 2060 rate: 456.883 efficiency: 0.000558125
  iBin: 3 Quality: 000000001000000110000 4144 rate: 426.669 efficiency: 0.00438527
  iBin: 4 Quality: 000000000001100000011 771 rate: 118.728 efficiency: 0.00608622
  iBin: 5 Quality: 000000000000100000011 259 rate: 52.8807 efficiency: 0.000956785

  iBin: 2 Quality: 000000011000011110000 12528 rate: 31.6062 efficiency: 0.00792084
  iBin: 3 Quality: 000000001000011110000 4336 rate: 25.113 efficiency: 0.0061051
  iBin: 5 Quality: 000000000110100001111 3343 rate: 16.032 efficiency: 0.00186837
  iBin: 6 Quality: 000000000100000111100 2108 rate: 15.5961 efficiency: 0.000578932
  iBin: 13 Quality: 000000011100011111100 14588 rate: 7.41043 efficiency: 0.00692087
  */
 /*
  hitsMask[3596] = true;
  hitsMask[2060] = true;
  hitsMask[4144] = true;
  hitsMask[771] = true;
  hitsMask[259] = true;
  hitsMask[12528] = true;
  hitsMask[4336] = true;
  hitsMask[3343] = true;
  hitsMask[2108] = true;
  hitsMask[14588] = true;
 */
  
  if(sysType.find("OMTF")!=std::string::npos){
    //return aL1Cand.type==L1Obj::OMTF_emu &&
    //aL1Cand.q==12 && aL1Cand.bx==0;
    
    bool lowPtVeto = aL1Cand.disc<250 &&
				  (hitsMask[aL1Cand.hits] ||
				   aL1Cand.refLayer==2 ||
				   aL1Cand.refLayer==3 ||
				   aL1Cand.refLayer==4 ||
				   aL1Cand.refLayer==5);

    //bool lowPtVeto = aL1Cand.disc<250 && hitsMask[aL1Cand.hits];
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
  float deltaR = 0.4, tmpR = 999;
  L1Obj selectedCand;

  for(auto aCand: myL1Coll){
    bool pass = passQuality(aCand ,sysType);
    if(!pass) continue;    
    double deltaEta = std::abs(genMuMom.Eta()-aCand.etaValue());    
    tmpR = deltaEta;
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
    /*
    std::cout<<" genMuMom.Pt(): "<<genMuMom.Pt()
	     <<" genMuMom.Eta(): "<<genMuMom.Eta()
	     <<" selectedCand.etaValue(): "<<selectedCand.etaValue()
	     <<" selectedCand.ptValue(): "<<selectedCand.ptValue()
	     <<std::endl;
    */
    myHistos_->fill1DHistogram("h1DLLH_Low", selectedCand.disc);    
    myHistos_->fill1DHistogram("h1DHitsPattern_Low_RefLayer", selectedCand.refLayer);    
  }
  if(genMuMom.Pt()>20 && selectedCand.ptValue()>=20 && sysType=="OMTF"){
    myHistos_->fill1DHistogram("h1DLLH_High", selectedCand.disc);   
    myHistos_->fill1DHistogram("h1DHitsPattern_High_RefLayer", selectedCand.refLayer);    
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
    if(pass && selectedCand.ptValue()< aCand.ptValue()) selectedCand = aCand;
  }

  float val = selectedCand.ptValue();
  if(selType=="Tot") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),val);

  bool pass = val>=20;
  if(selType=="VsEta") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),pass*genMuMom.Eta()+(!pass)*99);
  if(selType=="VsPt") myHistos_->fill2DHistogram(hName,genMuMom.Pt(),pass*genMuMom.Pt()+(!pass)*(-100));

  unsigned long hits = selectedCand.hits;
  bool passPtCut = selectedCand.ptValue()>=20;
  if(sysType=="OMTF" && selType=="VsQuality"){
    if(quality_index_map.find(hits)==quality_index_map.end()) quality_index_map[hits] = quality_index_map.size();
    int val = quality_index_map[hits];
    myHistos_->fill2DHistogram(hName, genMuMom.Pt(), pass*val+(!pass)*(-10));
    hName = "h2D"+sysType+"Quality20";    
    if(genMuMom.Pt()>40) {
      myHistos_->fill2DHistogram(hName, val, passPtCut);
    }
  }
  else if(selType=="VsQuality"){
     myHistos_->fill2DHistogram(hName, genMuMom.Pt(), 0);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillHistosForGenMuon(){

  //bool isOMTFAcceptance = fabs(genMuMom.Eta())>0.83 && fabs(genMuMom.Eta())<1.24;
  //if(!isOMTFAcceptance) return;
  
  bool isBMTFAcceptance = fabs(genMuMom.Eta())<0.83;
  if(!isBMTFAcceptance) return;

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
    if(std::abs(aGenObj.pdgId())!=13) continue;
    if(std::abs(aGenObj.status())!=1) continue;    
    genMuMom.SetPtEtaPhi(aGenObj.pt(),
			 aGenObj.eta(),
			 aGenObj.phi());
    fillHistosForGenMuon();
  }

  fillRateHisto("Vx","Tot");
  fillRateHisto("OMTF","Tot");
  fillRateHisto("kBMTF","Tot");

  fillRateHisto("Vx","VsPt");
  fillRateHisto("OMTF","VsPt");
  fillRateHisto("kBMTF","VsPt");

  fillRateHisto("Vx","VsEta");
  fillRateHisto("OMTF","VsEta");
  fillRateHisto("kBMTF","VsEta");

  fillRateHisto("OMTF","VsQuality");
  fillRateHisto("Vx","VsQuality");
  fillRateHisto("kBMTF","VsQuality");
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
