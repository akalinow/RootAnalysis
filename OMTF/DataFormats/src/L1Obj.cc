#include "L1Obj.h"
//AK #include "L1Trigger/L1TMuonOverlap/interface/OmtfName.h"
#include <bitset>
#include <iomanip>


L1Obj::L1Obj() : pt(0),eta(0),phi(0),
                 disc(0), 
                 bx(0),q(0), hits(0), charge(0), refLayer(0), 
                 type(NONE), 
                 iProcessor(-1), position(0) {};

std::ostream & operator<< (std::ostream &out, const L1Obj &o)
{
  out<<"L1Obj: ";
  switch (o.type) {
    case L1Obj::RPCb     : { out <<"RPCb    "; break; }
    case L1Obj::RPCf     : { out <<"RPCf    "; break; }
    case L1Obj::DT       : { out <<"DT      "; break; }
    case L1Obj::CSC      : { out <<"CSC     "; break; }
    case L1Obj::GMT      : { out <<"GMT     "; break; }
    case L1Obj::RPCb_emu : { out <<"RPCb_emu"; break; }
    case L1Obj::RPCf_emu : { out <<"RPCf_emu"; break; }
    case L1Obj::GMT_emu  : { out <<"GMT_emu "; break; }
    case L1Obj::OMTF     : { out <<"OMTF    "; break; }
    case L1Obj::OMTF_emu : { out <<"OMTF_emu"; break; }
    case L1Obj::BMTF     : { out <<"BMTF    "; break; }
    case L1Obj::EMTF     : { out <<"EMTF    "; break; }
    case L1Obj::uGMT     : { out <<"uGMT    "; break; }
    case L1Obj::uGMT_emu : { out <<"uGMT_emu"; break; }
    case L1Obj::NONE     : { out <<"NONE    "; break; }
    default: out <<"Unknown";
  };
  out << "  pt: "; 
    if (o.chargeValue()==1) out <<"+"; 
    else if (o.chargeValue()==-1) out <<"-";
    else out <<"?";
  out <<o.pt<<", eta: "<<o.eta;
  out <<", phi: ";  
  if (o.iProcessor >= 0){ out<<std::setw(2)<<o.iProcessor<<"_";} else { out<<std::setw(5); }
  out <<o.phi <<" (V: "<<std::setprecision(4)<<o.phiValue()<<std::setprecision(6)<<")";
  out <<", q: "<<o.q<<", bx: "<<o.bx;
  if (true || o.type ==  L1Obj::OMTF || o. type== L1Obj::OMTF_emu) {//AK true
      out <<" track: "<< std::bitset<29>(o.hits) 
          <<" disc: "<< std::bitset<12>(o.disc)
	  <<" ref layer: "<<o.refLayer
          <<", ";//AK << OmtfName(o.iProcessor, o.position);
  }
  return out;
}


ClassImp(L1Obj)
