#include "EventProxyOMTF.h"

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
EventProxyOMTF::EventProxyOMTF(){}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
EventProxyOMTF::~EventProxyOMTF(){}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
EventProxyBase* EventProxyOMTF::clone() const{

  return new EventProxyOMTF();
  
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
std::vector<OMTFHit> EventProxyOMTF::getHits() const{

  std::vector<OMTFHit> theHits;
  if(!hits || !hits->size()) return theHits;

  OMTFHit aHit;
  int sign = 1, tmp;
  for(unsigned int iHit=0; iHit<hits->size();++iHit){
    tmp = hits->at(iHit);
    if(hits->at(iHit)<0) {
      sign = 1;
      tmp = std::abs(tmp);
    }
    else sign = 0;
    
    aHit.iPhi = (tmp & 0xFFFF<<10) >> 10;
    aHit.iPhi *= std::pow(-1,sign);
    aHit.iLayer = (tmp & 0x1F<<5) >> 5;
    aHit.iHit = (tmp & 0x1F<<5);
    aHit.iQuality = (hitsQuality->at(iHit) & 0xFFFF<<10) >> 10;
    theHits.push_back(aHit);
  }

  return theHits;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void EventProxyOMTF::init(std::vector<std::string> const& iFileNames){

  treeName_ = "tOmtf";
  
  EventProxyBase::init(iFileNames);
  
  myEvent = 0;
  myGenObjColl = 0;
  myL1ObjColl = 0;
  myL1PhaseIIObjColl = 0;
  hits = 0;
  hitsQuality = 0;
  
  fChain->SetBranchStatus("*",1);
  fChain->SetMakeClass(0);
  fChain->SetBranchAddress("event",&myEvent);
  fChain->SetBranchAddress("genColl",&myGenObjColl);
  fChain->SetBranchAddress("l1ObjColl",&myL1ObjColl);
  fChain->SetBranchAddress("l1PhaseIIObjColl",&myL1PhaseIIObjColl);
  fChain->SetBranchAddress("hits",&hits);
  fChain->SetBranchAddress("hitsQuality",&hitsQuality);
  
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
std::ostream& operator<< (std::ostream& aStream, const OMTFHit& aHit){

  aStream<<"iLayer: "<<aHit.iLayer<<" iPhi: "<<aHit.iPhi<<" iQuality: "<<aHit.iQuality;

  return aStream;
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
