#include "GenObj.h"
#include <bitset>
#include <iomanip>


std::ostream & operator<< (std::ostream &out, const GenObj &o)
{
  out<<"GenObj: ";
  out<<" PDG id: "<<o.pdgId();
  out<<" charge: "<<o.charge();
  out<<"pt: "<<o.pt()<<" eta: "<<o.eta()<<" phi: "<<o.phi();
  out<<" status: "<<o.status();
  out<<" x coordinate of vertex position: "<<o.vx();
  out<<" y coordinate of vertex position: "<<o.vy();
  out<<" z coordinate of vertex position: "<<o.vz();
  out<<" beta: "<<o.beta();




  out<<std::endl;

  return out;
}


ClassImp(GenObj)
