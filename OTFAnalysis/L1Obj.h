#ifndef L1Obj_H
#define L1Obj_H
#include "TObject.h"
#include <iostream>

class L1Obj : public TObject {

 public:
  
  enum TYPE { NONE, RPCb, RPCf, DT, CSC, GMT, RPCb_emu, RPCf_emu, GMT_emu, OTF, simOMTF, hwOMTF };

  float pt, eta, phi;
  float disc;
  int   bx, q, hits, charge, refLayer;
  int   iProcessor, mtfType;
  TYPE  type;

  L1Obj();
  bool isValid() const { return q >= 0;}

  ClassDef(L1Obj,2)
};


std::ostream & operator<< (std::ostream &out, const L1Obj &o);

#endif
