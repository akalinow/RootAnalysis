#ifndef GenObj_H
#define GenObj_H
#include <ostream>

#include "TObject.h"
// #include "DataFormats/Math/interface/LorentzVector.h"

class GenObj : public TObject {
public:
 GenObj(float pt=0., float eta=0., float phi=0.,float mass=0., 
	int charge=0, int pdgid=0, int st=0, double vx=0., double vy=0., double vz=0., double beta=0):
  _pt(pt),_eta(eta),_phi(phi),_mass(mass),_charge(charge),_id(pdgid),_status(st),_vx(vx),_vy(vy),_vz(vz),_beta(beta){}
  virtual ~GenObj() {}
public:
  float pt() const { return _pt;}
  float eta() const { return _eta;}
  float phi() const { return _phi;}
  float mass() const { return _mass;}
  int pdgId() const { return _id;}
  int status() const { return _status;}
  int charge() const { return _charge;}
  double vx() const { return _vx;}
  double vy() const { return _vy;}
  double vz() const { return _vz;}
  double beta() const { return _beta;}


private:  
  float _pt,_eta,_phi,_mass; 
  int _charge,_id,_status; 
  double _vx, _vy, _vz, _beta;

public:
  ClassDef(GenObj,5);
};

std::ostream & operator<< (std::ostream &out, const GenObj &o);

#endif
