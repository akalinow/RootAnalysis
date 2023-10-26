#include "GenObj.h"
#include <bitset>
#include <iomanip>
#include <cmath>

void GenObj::setVertexXYZ(double x, double y, double z){

  _vx = x;
  _vy = y;
  _vz = z;
}

void GenObj::setPtEtaPhiM(double pt, double eta, double phi, double m){

  _pt = pt;
  _eta = eta;
  _phi = phi;
  _mass = m;
  
  double theta = 2.0*atan(exp(_eta));
  double p = std::abs(_pt/sin(theta));
  _beta = p/sqrt(p*p + _mass*_mass);
}

std::ostream & operator<< (std::ostream &out, const GenObj &o)
{
  out<<"GenObj: ";
  out<<" PDG id: "<<o.pdgId();
  out<<" charge: "<<o.charge();
  out<<"pt: "<<o.pt()<<" eta: "<<o.eta()<<" phi: "<<o.phi();
  out<<" status: "<<o.status()<<" motherId: "<<o.motherId();
  out<<" vertex: (";
  out<<o.vx()<<", "<<o.vy()<<", "<<o.vz()<<")";
  out<<std::endl;

  return out;
}


ClassImp(GenObj)
