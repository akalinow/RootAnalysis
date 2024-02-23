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

  double dxy() const;
  double dz() const;

private:  
  float _pt{0},_eta{-99},_phi{-99},_mass{0}; 
  int _charge{0},_id{0},_status{0},_mid{0}; 
  double _vx{-99}, _vy{-99}, _vz{-99}, _beta{0};

public:
  ClassDef(GenObj,4);
};

std::ostream & operator<< (std::ostream &out, const GenObj &o);

#endif
