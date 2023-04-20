#include "L1PhaseIIObj.h"
#include <bitset>
#include <iomanip>


L1PhaseIIObj::L1PhaseIIObj() : pt(0),eta(0),phi(0),
                 disc(0), 
                 bx(0),q(0), hits(0), charge(0), refLayer(0), 
                 type(NONE), 
                 iProcessor(-1), position(0) {};

std::ostream & operator<< (std::ostream &out, const L1PhaseIIObj &o)
{
  out<<"L1PhaseIIObj: ";
  switch (o.type) {
    case L1PhaseIIObj::RPCb     : { out <<"RPCb    "; break; }
    case L1PhaseIIObj::RPCf     : { out <<"RPCf    "; break; }
    case L1PhaseIIObj::DT       : { out <<"DT      "; break; }
    case L1PhaseIIObj::CSC      : { out <<"CSC     "; break; }
    case L1PhaseIIObj::GMT      : { out <<"GMT     "; break; }
    case L1PhaseIIObj::RPCb_emu : { out <<"RPCb_emu"; break; }
    case L1PhaseIIObj::RPCf_emu : { out <<"RPCf_emu"; break; }
    case L1PhaseIIObj::GMT_emu  : { out <<"GMT_emu "; break; }
    case L1PhaseIIObj::OMTF     : { out <<"OMTF    "; break; }
    case L1PhaseIIObj::OMTF_emu : { out <<"OMTF_emu"; break; }
    case L1PhaseIIObj::BMTF     : { out <<"BMTF    "; break; }
    case L1PhaseIIObj::EMTF     : { out <<"EMTF    "; break; }
    case L1PhaseIIObj::uGMT     : { out <<"uGMT    "; break; }
    case L1PhaseIIObj::uGMT_emu : { out <<"uGMT_emu"; break; }
    case L1PhaseIIObj::NONE     : { out <<"NONE    "; break; }
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
  if (o.type ==  L1PhaseIIObj::OMTF || o. type== L1PhaseIIObj::OMTF_emu) {
      out <<" track: "<< std::bitset<29>(o.hits) 
          <<" disc: "<< std::bitset<12>(o.disc)
          // <<", "<< OmtfName(o.iProcessor, o.position);
          <<", ";
  }
  return out;
}


// ClassImp(L1PhaseIIObj);
