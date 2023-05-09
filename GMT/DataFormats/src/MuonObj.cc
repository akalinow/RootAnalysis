#include "MuonObj.h"
#include <bitset>
#include <iomanip>


std::ostream & operator<< (std::ostream &out, const MuonObj &o)
{
  out<<"RecoMuonObj: ";
  out<<" charge: "<<o.charge();
  out<<" pt:     "<<o.pt() <<"  eta:  "<<o.eta() <<" phi: "<<o.phi();
  out<<" RPC Hits : "<<o.nrpchits();
  out<<" DT  Hits : "<<o.ndthits();
  out<<" Medium ID : "<<o.mediumID();
  out<<" l1Eta     : "<<o.l1eta();
  out<<" Hlt  matching  : "<<o.matchedhlt();
  out<<" Hlt Iso matching  : "<<o.matchedisohlt();
  
  out<<std::endl;

  return out;
}


ClassImp(MuonObj)               
