#include "L1Obj.h"
#include <bitset>
#include <iomanip>

namespace { 
  double modulo2PI (double phi) { 
    while (phi > 2*M_PI) phi -= 2*M_PI;
    while (phi < 0.) phi += 2*M_PI;
    return phi;
  }
}


L1Obj::L1Obj() : pt(0),eta(0),phi(0),
                 disc(0), 
                 bx(0),q(0), hits(0), charge(0), refLayer(0), 
                 type(NONE), 
                 iProcessor(-1), position(0) {};
                 
double L1Obj::ptValue() const { return type==uGMTPhase2_emu ?  pt : (pt-1.)/2.; }
double L1Obj::etaValue() const { return type==uGMTPhase2_emu ? eta : eta/240.*2.61; }
double L1Obj::phiValue() const {
    if (type==OMTF || type==OMTF_emu || type==EMTF) 
    return modulo2PI( ( (15.+iProcessor*60.)/360. + phi/576. ) *2*M_PI) ;  
    else if (type==BMTF) return modulo2PI( ( (-15.+iProcessor*30.)/360. + phi/576. ) *2*M_PI);
    else if (type==uGMT || type==uGMT_emu) return modulo2PI((phi/576.)*2*M_PI);
    else if (type==uGMTPhase2_emu) return modulo2PI(phi);
    else return 9999.;
  }
int L1Obj::chargeValue() const { return type==uGMTPhase2_emu ? charge : pow(-1,charge); }

double L1Obj::ptUnconstrainedValue() const { return ptUnconstrained - 1;}
double L1Obj::z0Value() const { return z0;}
double L1Obj::d0Value() const { return d0;}

               

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
    case L1Obj::uGMTPhase2_emu : { out <<"uGMTPhase2_emu"; break; }
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
  if (o.type ==  L1Obj::OMTF || o. type== L1Obj::OMTF_emu) {
      out <<" track: "<< std::bitset<29>(o.hits) 
          <<" disc: "<< std::bitset<12>(o.disc);          
  }
  return out;
}


ClassImp(L1Obj)
