#include <cstdlib> 
#include <string>
#include <omp.h>
#include <bitset>
#include <sstream>

#include <iostream>

#include "TreeAnalyzer.h"
#include "OTFAnalyzer.h"
#include "EventProxyOTF.h"


int iCandOTF = 0;
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::OTFAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::~OTFAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OTFHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* OTFAnalyzer::clone() const{

  OTFAnalyzer* clone = new OTFAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::finalize(){ 
  /*
 std::cout<<"tmpMap.size(): "<<tmpMap.size()<<std::endl;
 
 std::ostringstream stringStr;
 TH2F *h = myHistos_->get2DHistogram("h2DRateVsQualityOtf",true);
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
  */

  myHistos_->finalizeHistograms(0,1.0); 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::registerCuts(){

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
bool OTFAnalyzer::passQuality(std::vector<L1Obj> * myL1Coll,
		 const std::string & sysType, 
		 int iCand){

  if(sysType.find("Gmt")!=std::string::npos){
    return myL1Coll->size()>iCand &&
      //myL1Coll->operator[](iCand).eta>0.8 &&
      //myL1Coll->operator[](iCand).eta<1.3 &&
      true;
  }
  
  if(sysType.find("Otf")!=std::string::npos){
  return myL1Coll->size()>iCand &&
    ///Barrel (l1t::tftype::bmtf)
    /*
    myL1Coll->operator[](iCand).q!=9728 &&
    myL1Coll->operator[](iCand).q!=8960 &&
    myL1Coll->operator[](iCand).q!=3084 &&
    myL1Coll->operator[](iCand).q!=771 &&
    myL1Coll->operator[](iCand).q!=3596 &&
    myL1Coll->operator[](iCand).q!=7180 &&
    myL1Coll->operator[](iCand).q!=7168 &&
    */
    ///Endcap (l1t::tftype::emtf)
    /*
    myL1Coll->operator[](iCand).q!=12288 &&
    myL1Coll->operator[](iCand).q!=6144 &&
    myL1Coll->operator[](iCand).q!=9216 &&
    myL1Coll->operator[](iCand).q!=3072 &&
    myL1Coll->operator[](iCand).q!=3264 &&
    myL1Coll->operator[](iCand).q!=1036 &&
    myL1Coll->operator[](iCand).q!=8195  &&
    myL1Coll->operator[](iCand).q!=4099  &&
    myL1Coll->operator[](iCand).q!=10240  &&
    myL1Coll->operator[](iCand).q!=2051  &&
    myL1Coll->operator[](iCand).q!=3084  &&
    myL1Coll->operator[](iCand).q!=1216  &&
    myL1Coll->operator[](iCand).q!=2240  &&
    myL1Coll->operator[](iCand).q!=5120  &&
    myL1Coll->operator[](iCand).q!=2060  &&
    myL1Coll->operator[](iCand).q!=3276  &&
    myL1Coll->operator[](iCand).q!=2252  &&
    myL1Coll->operator[](iCand).q!=8240  &&
    myL1Coll->operator[](iCand).q!=16380 &&
    //myL1Coll->operator[](iCand).q!=16368 &&
    //myL1Coll->operator[](iCand).q!=16188 &&
    //myL1Coll->operator[](iCand).q!=14332 &&
    //myL1Coll->operator[](iCand).q!=15356 &&
    */
    ///Barrel (l1t::tftype::omtf_pos)           
    myL1Coll->operator[](iCand).q!=99840 &&
    myL1Coll->operator[](iCand).q!=34304 &&
    myL1Coll->operator[](iCand).q!=3075 &&
    myL1Coll->operator[](iCand).q!=36928 &&
    ////
    myL1Coll->operator[](iCand).q!=12300 && 
    
    ///Endcap (l1t::tftype::omtf_pos)
    myL1Coll->operator[](iCand).q!=98816 &&
    myL1Coll->operator[](iCand).q!=98944 &&

    myL1Coll->operator[](iCand).q!=33408 &&
    myL1Coll->operator[](iCand).q!=66688 && 
    myL1Coll->operator[](iCand).q!=66176 && 
    ///    
    true;
  }
  else return myL1Coll->size();
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::fillTurnOnCurve(const int & iPtCut,
				  const std::string & sysType,
				  const std::string & selType){

  int ptCut = OTFHistograms::ptBins[iPtCut];
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
    iCand = iCandOTF;
    myL1Coll = &(theEvent->l1ObjectsOtf);
    if(myL1Coll->size()>1 && myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
    hName = "h2DOtf"+selType;
  }
  
  bool qualityCut = passQuality(myL1Coll,sysType,iCand);
  bool pass = myL1Coll->size() && myL1Coll->operator[](iCand).pt>=ptCut && qualityCut;
  
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
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::fillRateHisto(const std::string & sysType,
				const std::string & selType){

	int iCut = 2;
	int iCand = 0;
	float ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsGmt[iCut]];

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
	  iCand = iCandOTF;
	  myL1Coll = &(theEvent->l1ObjectsOtf);
	  if(myL1Coll->size()>1){
	    ///Take lower quality
	    if(myL1Coll->operator[](0).q<myL1Coll->operator[](1).q) iCand = 1;
	    ///take lower pt in case of the same quality
	    else if (myL1Coll->operator[](0).pt>myL1Coll->operator[](1).pt) iCand = 1;
	  }
	  hName = "h2DRate"+selType+"Otf";
	  ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsOtf[iCut]];

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
	if(pass) q = myL1Coll->operator[](iCand).q;

	std::bitset<18> hitsWord(q);
	if(sysType=="Otf" && selType=="VsQuality"){
	  if(tmpMap.find(q)==tmpMap.end()) tmpMap[q] = tmpMap.size();
	  int val = tmpMap[q];
	  myHistos_->fill2DHistogram(hName,theEvent->pt,pass*val+(!pass)*(-10));
	}
	if(sysType=="Gmt" && selType=="VsQuality") myHistos_->fill2DHistogram(hName,theEvent->pt,pass*q+(!pass)*(-10));
	

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::fillGhostHisto(const std::string & sysType,
				 const std::string & selType){

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
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool OTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();
  /*
  std::bitset<17> bits(99840);
  std::cout<<99840<<" bits: "<<bits.to_string()<<std::endl;

  bits = std::bitset<17>(34304);
  std::cout<<34304<<" bits: "<<bits.to_string()<<std::endl;

  bits = std::bitset<17>(3075); 
  std::cout<<3075<<" bits: "<<bits.to_string()<<std::endl;

  bits = std::bitset<17>(36928); 
  std::cout<<36928<<" bits: "<<bits.to_string()<<std::endl;

 bits = std::bitset<17>(12300); 
  std::cout<<12300<<" bits: "<<bits.to_string()<<std::endl;

 bits = std::bitset<17>(98816); 
  std::cout<<98816<<" bits: "<<bits.to_string()<<std::endl;

  bits = std::bitset<17>(98944); 
  std::cout<<98944<<" bits: "<<bits.to_string()<<std::endl;

 bits = std::bitset<17>(33408); 
  std::cout<<98944<<" bits: "<<bits.to_string()<<std::endl;

 bits = std::bitset<17>(66688); 
  std::cout<<98944<<" bits: "<<bits.to_string()<<std::endl;

 bits = std::bitset<17>(66176); 
  std::cout<<98944<<" bits: "<<bits.to_string()<<std::endl;

  exit(0);
  */

  std::string selType = "";
  std::string sysTypeGmt="Gmt";
  std::string sysTypeOtf="Otf";
  std::string sysTypeRpc="Rpc";
  std::string sysTypeOther="Other";
  
  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events[omp_get_thread_num()];

  if(theEvent->eta<0.83 || theEvent->eta>1.24) return true;

  for(int iCut=0;iCut<22;++iCut){
	  if(iCut>0 && iCut<14) continue;
	  fillTurnOnCurve(iCut,sysTypeGmt,selType);
	  fillTurnOnCurve(iCut,sysTypeOtf,selType);
	  fillTurnOnCurve(iCut,sysTypeRpc,selType);
	  fillTurnOnCurve(iCut,sysTypeOther,selType);
  }


  ////////////////
  int iCut = 2;
  bool pass = false;

  for(int iType=0;iType<=3;++iType){
	  float ptCut = OTFHistograms::ptBins[OTFHistograms::ptCutsGmt[iCut]];
	  if(iType==0) pass = theEvent->pt>(ptCut + 20);
	  if(iType==1) pass = theEvent->pt>ptCut && theEvent->pt<(ptCut+5);
	  if(iType==2) pass = theEvent->pt<10;
	  if(!pass) continue;

	  selType = std::string(TString::Format("Type%d",iType));
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeGmt,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeRpc,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsGmt[iCut],sysTypeOther,selType);
	  fillTurnOnCurve(OTFHistograms::ptCutsOtf[iCut],sysTypeOtf,selType);	  
  }
  /////////////////

  fillRateHisto("Otf","Tot");
  fillRateHisto("Gmt","Tot");

  fillRateHisto("Gmt","VsEta");
  fillRateHisto("Gmt","VsPt");

  fillRateHisto("Otf","VsEta");
  fillRateHisto("Otf","VsPt");

  fillRateHisto("Otf","VsQuality");
  fillRateHisto("Gmt","VsQuality");

  //fillGhostHisto("Otf","VsProcessor");
  //fillGhostHisto("Gmt","VsProcessor");

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  OTFAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  OTFAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

