#ifndef EventData_H
#define EventData_H
#include "TObject.h"
#include <iostream>

#include "L1Obj.h"

///Small class to create a TTree
///holding trigger response.

struct EventData: public TObject{

  EventData(){};
  ~EventData(){}

  void clear();
  
  ///Header
  float weight;
  ///Generated kinematics
  float pt;
  float eta;
  float phi;

  float pt1;
  float eta1;
  float phi1;

  float phiHit;
  float etaHit;
  int charge;

  inline float etaMuon(unsigned int iMuon){
    if(iMuon==0) return eta;
    else if(iMuon==1) return eta1;
    else return 99;
  };

  inline float phiMuon(unsigned int iMuon){
    if(iMuon==0) return phi;
    else if(iMuon==1) return phi1;
    else return 99;
  };
  
  inline float ptMuon(unsigned int iMuon){
    if(iMuon==0) return pt;
    else if(iMuon==1) return pt1;
    else return 0;
  };

  std::vector<L1Obj> l1ObjectsOtf;
  std::vector<L1Obj> l1ObjectsGmt;
  std::vector<L1Obj> l1ObjectsRpc;
  std::vector<L1Obj> l1ObjectsOther;

  ClassDef(EventData,2)
};
///////////////////////////////////////
#endif
