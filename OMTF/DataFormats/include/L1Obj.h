#ifndef L1Obj_H
#define L1Obj_H
#include "TObject.h"
#include <iostream>
#include <math.h>

class L1Obj : public TObject {

public:
  
  enum TYPE { NONE, RPCb, RPCf, DT, CSC, GMT, RPCb_emu, RPCf_emu, GMT_emu, OMTF, OMTF_emu, BMTF, EMTF, uGMT, uGMT_emu, uGMTPhase2_emu};

  double pt{0}, eta{0}, phi{0};
  double ptUnconstrained{0};
  double z0{0}, d0{0};
  int disc{0};
  int   bx{0}, q{0}, hits{0}, charge{0}, refLayer{0};
  TYPE  type{NONE};
  int   iProcessor{-1}, position{0};

  L1Obj(){};

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

bool operator< (const L1Obj &a, const L1Obj &b);
bool operator> (const L1Obj &a, const L1Obj &b);

std::ostream & operator<< (std::ostream &out, const L1Obj &o);

bool operator< (const L1Obj &a, const L1Obj &b);
bool operator> (const L1Obj &a, const L1Obj &b);

#endif
