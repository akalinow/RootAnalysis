#ifndef L1Obj_H
#define L1Obj_H
#include "TObject.h"
#include <iostream>
#include <math.h>

class L1Obj : public TObject {

public:
  
  enum TYPE { NONE, RPCb, RPCf, DT, CSC, GMT, RPCb_emu, RPCf_emu, GMT_emu, OMTF, OMTF_emu, BMTF, EMTF, uGMT, uGMT_emu, uGMTPhase2_emu};

  double pt{-999}, eta{-999}, phi{-999};
  double ptUnconstrained{-999};
  double z0{-999}, d0{-999};
  int disc{-999};
  int   bx{-999}, q{-999}, hits{-999}, charge{-999}, refLayer{-999};
  TYPE  type{NONE};
  int   iProcessor{-999}, position{-999};

  L1Obj();

  bool isValid() const { return type!=NONE && pt>0;}

  double ptValue() const;   
  double etaValue() const; 
  double phiValue() const; 
  int chargeValue() const;
  
  double ptUnconstrainedValue() const; 
  double z0Value() const;
  double d0Value() const;

  ClassDef(L1Obj,6)
};


std::ostream & operator<< (std::ostream &out, const L1Obj &o);

#endif
