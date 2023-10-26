#ifndef GenObj_H
#define GenObj_H
#include <ostream>

#include "TObject.h"

class GenObj : public TObject {
public:

GenObj(int charge=0, int pdgid=0, int st=0, int mother=0):_charge(charge),_id(pdgid),_status(st),_mid(mother){}

  virtual ~GenObj() {}
public:

  void setVertexXYZ(double x, double y, double z);
  void setPtEtaPhiM(double pt, double eta, double phi, double m);

  float pt() const { return _pt;}
  float eta() const { return _eta;}
  float phi() const { return _phi;}
  float mass() const { return _mass;}
  int pdgId() const { return _id;}
  int status() const { return _status;}
  int motherId() const { return _mid;}
  int charge() const { return _charge;}
  double vx() const { return _vx;}
  double vy() const { return _vy;}
  double vz() const { return _vz;}
  double beta() const { return _beta;}


private:  
  float _pt,_eta,_phi,_mass; 
  int _charge,_id,_status,_mid; 
  double _vx, _vy, _vz, _beta;

public:
  ClassDef(GenObj,4);
};

std::ostream & operator<< (std::ostream &out, const GenObj &o);

#endif
