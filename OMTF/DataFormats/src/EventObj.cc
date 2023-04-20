#include "EventObj.h"


std::ostream & operator<< (std::ostream &out, const EventObj &o) {
  out <<"run: "<< o.run <<" event: "<<o.id<<" lumi: "<<o.lumi; 
  return out;
}

ClassImp(EventObj)
