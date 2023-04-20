#include "MuonObj.h"
#include <bitset>
#include <iomanip>


std::ostream & operator<< (std::ostream &out, const MuonObj &o)
{
  out<<"RecoMuonObj: ";
  out<<" charge: "<<o.charge();
  out<<"pt: "<<o.pt()<<" eta: "<<o.eta()<<" phi: "<<o.phi();
  out<<" RPC Hits : "<<o.nrpchits();
  out<<" DT  Hits : "<<o.ndthits();
  out<<" Medium ID : "<<o.mediumID();




  out<<std::endl;

  return out;
}


ClassImp(MuonObj)               
