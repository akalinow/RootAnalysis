#include <cstdlib>
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "OMTFAnalyzer.h"
#include "EventProxyOMTF.h"


int iCandOMTF = 0;
int iCandNN = 1;
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OMTFAnalyzer::OMTFAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

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

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OMTFHistograms(aDir, selectionFlavours_);

  registerCuts();

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

  if(omp_get_max_threads()==1){
    std::cout<<"tmpMap.size(): "<<tmpMap.size()<<std::endl;

    std::ostringstream stringStr;
    TH2F *h = myHistos_->get2DHistogram("h2DRateVsQualityOtf",true);
    if(h){
      for(int iBin=1;iBin<h->GetYaxis()->GetNbins();++iBin){
	stringStr.str("");
	stringStr<<iBin;
	h->GetYaxis()->SetBinLabel(iBin,stringStr.str().c_str());
      }

      for(auto it: tmpMap){
	int iBinX = h->GetYaxis()->FindFixBin(it.second);
	if(iBinX>=h->GetYaxis()->GetNbins()) continue;

	std::bitset<18> bits(it.first);

	stringStr.str("");
	stringStr<<it.first;
	std::string label = bits.to_string()+" "+stringStr.str();
	h->GetYaxis()->SetBinLabel(iBinX,label.c_str());
      }
    }
  }

  myHistos_->finalizeHistograms(0,1.0);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::registerCuts(){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
    mySelections_->push_back("HLT"+selectionFlavours_[i]);
  /*
 ///Register tree variables
 for(unsigned int i=0;i<triggerItemNames_.size();++i){
   treeVariables_[triggerItemNames_[i]] = -999.0;
 }
 */
 }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::passQuality(std::vector<L1Obj> * myL1Coll,
		 const std::string & sysType,
		 unsigned int iCand){

  if(sysType.find("Gmt")!=std::string::npos){
    return myL1Coll->size()>iCand &&
      fabs(myL1Coll->operator[](iCand).eta)>0.8 &&
      fabs(myL1Coll->operator[](iCand).eta)<1.3 &&
      true;
  }

  if(sysType.find("Otf")!=std::string::npos){
  return myL1Coll->size()>iCand &&
    myL1Coll->operator[](iCand).type==L1Obj::OMTF_emu &&
         myL1Coll->operator[](iCand).q>4 &&
    true;
  }
  else return myL1Coll->size();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillTurnOnCurve(const int & iPtCut,
				  const std::string & sysType,
				  const std::string & selType){
  /*

  int ptCut = OMTFHistograms::ptBins[iPtCut];
  int iCand = 0;

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
    iCand = iCandOMTF;
    myL1Coll = &(theEvent->l1ObjectsOtf);
    //if(myL1Coll->size()>1 && myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
    hName = "h2DOtf"+selType;
  }
  if(sysType=="NN") {
    iCand = iCandNN;
    myL1Coll = &(theEvent->l1ObjectsOtf);
    //if(myL1Coll->size()>1 && myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
    hName = "h2DNN"+selType;
  }


  ///Find best matching L1 candidate
  float deltaR = 0.5, tmpR = 999;
  for(unsigned int index=0;index<myL1Coll->size();++index){
    if(!passQuality(myL1Coll,sysType,index)) continue;
    L1Obj aCand = myL1Coll->operator[](index);
    tmpR = sqrt(0*pow(theEvent->phi-aCand.phi,2) + //phi is not propaged to vertex
		pow(theEvent->eta-aCand.eta,2));
    if(tmpR<deltaR){
      deltaR = tmpR;
      iCand = index;
    }
  }

  bool qualityCut = passQuality(myL1Coll,sysType,iCand);
  bool pass = myL1Coll->size() && myL1Coll->operator[](iCand).pt>=ptCut && qualityCut && deltaR<100.2;///TEST

  std::string tmpName = hName+"Pt"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->pt, pass);

  ///Fill histos for eff vs eta/phi only for events at the plateau.
  if(selType.size()==0 && theEvent->pt<(ptCut + 20)) return;

  tmpName = hName+"EtaHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->etaHit, pass);

  tmpName = hName+"PhiHit"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phiHit, pass);

  tmpName = hName+"EtaVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->eta, pass);

  tmpName = hName+"PhiVx"+std::to_string(ptCut);
  myHistos_->fill2DHistogram(tmpName,theEvent->phi, pass);

  tmpName = hName+"Quality"+std::to_string(ptCut);
  int q = -10;
  if(myL1Coll->size()) q = myL1Coll->operator[](iCand).hits;
  if(tmpMap.find(q)==tmpMap.end()) tmpMap[q] = tmpMap.size();
  int xPosition = tmpMap[q];
  myHistos_->fill2DHistogram(tmpName,xPosition,pass);

  */
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillRateHisto(const std::string & sysType,
				const std::string & selType){
  /*

	int iCut = 2;
	int iCand = 0;
	float ptCut = OMTFHistograms::ptBins[OMTFHistograms::ptCutsGmt[iCut]];

	std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsGmt);
	std::string hName = "h2DRate"+selType+"Gmt";


	if(sysType=="Rpc"){
			myL1Coll = &(theEvent->l1ObjectsRpc);
			hName = "h2DRate"+selType+"Rpc";
	}
	if(sysType=="Other"){
			myL1Coll = &(theEvent->l1ObjectsOther);
			hName = "h2DRate"+selType+"Other";
	}
	if(sysType=="Otf") {
	  iCand = iCandOMTF;
	  myL1Coll = &(theEvent->l1ObjectsOtf);
	  if(myL1Coll->size()>1){
	    ///Take higher quality
	    if(myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
	    ///take lower pt in case of the same quality
	    else if (myL1Coll->operator[](0).pt>myL1Coll->operator[](1).pt) iCand = 1;
	  }
	  hName = "h2DRate"+selType+"Otf";
	  ptCut = OMTFHistograms::ptBins[OMTFHistograms::ptCutsOtf[iCut]];
	}
  if(sysType=="NN") {
	  iCand = iCandNN;
	  myL1Coll = &(theEvent->l1ObjectsOtf);
	  if(myL1Coll->size()>1){
	    ///Take higher quality
	    if(myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
	    ///take lower pt in case of the same quality
	    else if (myL1Coll->operator[](0).pt>myL1Coll->operator[](1).pt) iCand = 1;
	  }
	  hName = "h2DRate"+selType+"NN";
	  ptCut = OMTFHistograms::ptBins[OMTFHistograms::ptCutsOtf[iCut]];
	}

	bool qualityCut = passQuality(myL1Coll,sysType,iCand);
	bool pass = myL1Coll->size() && qualityCut;
	float val = 0;

	if(pass) val = myL1Coll->operator[](iCand).pt;
        if(selType=="Tot") myHistos_->fill2DHistogram(hName,theEvent->pt,val);

	///Rate vs selected variable is plotted for given pt cut.
	pass = pass && (myL1Coll->operator[](iCand).pt>=ptCut);
	if(selType=="VsEta") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*theEvent->eta+(!pass)*99);
	if(selType=="VsPt") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*theEvent->pt+(!pass)*(-100));
	int q = -10;
	if(pass) q = myL1Coll->operator[](iCand).hits;

	std::bitset<18> hitsWord(q);
	if(sysType=="Otf" && selType=="VsQuality"){
	  if(tmpMap.find(q)==tmpMap.end()) tmpMap[q] = tmpMap.size();
	  int val = tmpMap[q];
	  float ptGen = theEvent->pt;
	  //if(ptGen<1) ptGen = 1.1;
	  myHistos_->fill2DHistogram(hName,ptGen,pass*val+(!pass)*(-10));
	}
	if(sysType=="Gmt" && selType=="VsQuality") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*q+(!pass)*(-10));
  */

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::fillGhostHisto(const std::string & sysType,
				 const std::string & selType){
  
  /*
  std::vector<L1Obj> * myL1Coll = &(theEvent->l1ObjectsGmt);
  std::string hName = "h2DGhostsVsProcessor"+sysType;

  if(sysType=="Otf") myL1Coll = &(theEvent->l1ObjectsOtf);
  if(sysType=="Gmt") myL1Coll = &(theEvent->l1ObjectsGmt);

  std::vector<unsigned int> myCounts(6), myCountsMinus(6), myCountsPlus(6);

  for(unsigned int iCandTmp=0;iCandTmp<myL1Coll->size();++iCandTmp){
    if(!passQuality(myL1Coll,sysType,iCandTmp)) continue;
    unsigned int iProcessor = myL1Coll->operator[](iCandTmp).phi;
    if(sysType=="Gmt") iProcessor = 2;
    unsigned int iCharge = myL1Coll->operator[](iCandTmp).charge;
    ++myCounts[iProcessor];
    if(iCharge==0) ++myCountsMinus[iProcessor];
    if(iCharge==1) ++myCountsPlus[iProcessor];
  }


  for(unsigned int iProcessor=0;iProcessor<6;++iProcessor){
    myHistos_->fill2DHistogram(hName,iProcessor,myCounts[iProcessor]);
    myHistos_->fill2DHistogram(hName+"SS",iProcessor,abs(myCountsPlus[iProcessor]-myCountsMinus[iProcessor]));
    myHistos_->fill2DHistogram(hName+"OS",iProcessor,abs(myCountsPlus[iProcessor]*(myCountsPlus[iProcessor]<2)
							 +myCountsMinus[iProcessor]*(myCountsMinus[iProcessor]<2)));
  }
  */
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OMTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOMTF & myProxy = static_cast<const EventProxyOMTF&>(iEvent);
  std::cout<<"addr: "<<*(myEvent.myGenObjColl)<<std::endl;

  std::string selType = "";
  std::string sysTypeGmt="Gmt";
  std::string sysTypeOtf="Otf";
  std::string sysTypeNN="NN";
  std::string sysTypeRpc="Rpc";
  std::string sysTypeOther="Other";


  if(theEvent->pt<0.01){

    fillRateHisto("Otf","Tot");
    fillRateHisto("Gmt","Tot");

    fillRateHisto("Gmt","VsEta");
    fillRateHisto("Gmt","VsPt");

    fillRateHisto("Otf","VsEta");
    fillRateHisto("Otf","VsPt");

    fillRateHisto("Otf","VsQuality");
    fillRateHisto("Gmt","VsQuality");
  }
  else if(fabs(theEvent->eta)>0.83 && fabs(theEvent->eta)<1.24){
    for(int iCut=0;iCut<22;++iCut){
      if(iCut>0 && iCut<14) continue;
      fillTurnOnCurve(iCut,sysTypeGmt,selType);
      fillTurnOnCurve(iCut,sysTypeOtf,selType);
      fillTurnOnCurve(iCut,sysTypeNN,selType);
      fillTurnOnCurve(iCut,sysTypeRpc,selType);
      fillTurnOnCurve(iCut,sysTypeOther,selType);
    }
    ////////////////
    int iCut = 2;
    bool pass = false;
    for(int iType=0;iType<=3;++iType){
      float ptCut = OMTFHistograms::ptBins[OMTFHistograms::ptCutsGmt[iCut]];
      //if(iType==0) pass = theEvent->pt>(ptCut + 20);
      if(iType==0) pass = theEvent->pt>24;
      if(iType==1) pass = theEvent->pt>ptCut && theEvent->pt<(ptCut+5);
      if(iType==2) pass = theEvent->pt<10;
      if(!pass) continue;

      selType = std::string(TString::Format("Type%d",iType));
      fillTurnOnCurve(OMTFHistograms::ptCutsGmt[iCut],sysTypeGmt,selType);
      fillTurnOnCurve(OMTFHistograms::ptCutsGmt[iCut],sysTypeRpc,selType);
      fillTurnOnCurve(OMTFHistograms::ptCutsGmt[iCut],sysTypeOther,selType);
      fillTurnOnCurve(OMTFHistograms::ptCutsOtf[iCut],sysTypeOtf,selType);
      fillTurnOnCurve(OMTFHistograms::ptCutsOtf[iCut],sysTypeNN,selType);
    }
    /////////////////
  }
  */
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  OMTFAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  OMTFAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OMTFAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
