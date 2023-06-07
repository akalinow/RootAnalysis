#include "MuonObj.h"
#include <bitset>
#include <iomanip>


std::ostream & operator<< (std::ostream &out, const MuonObj &o)
{
  out<<"RecoMuonObj: ";
  out<<" ch: "<<o.charge();
  out<<" pt,eta,phi: ("<<o.pt() <<", "<<o.eta() <<", "<<o.phi()<<")";
  out<<" Hits: RPC|DT:"<<o.nrpchits()<<"|"<<o.ndthits();
  out<<" ID medium/tight: "<<o.mediumID()<<"|"<<o.tightID();
  out<<" l1 eta,phi: "<<o.l1eta()<<", "<<o.l1phi();
  out<<" match HLT|IsoHLT: "<<o.matchedhlt()<<"|"<<o.matchedisohlt();
  
  out<<std::endl;
  return out;
}


ClassImp(MuonObj)               
