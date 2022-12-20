#ifndef L1PhaseIIObj_H
#define L1PhaseIIObj_H
#include "TObject.h"
#include <iostream>
#include <math.h>


class L1PhaseIIObj : public TObject {

public:
  
  enum TYPE { NONE, RPCb, RPCf, DT, CSC, GMT, RPCb_emu, RPCf_emu, GMT_emu, OMTF, OMTF_emu, BMTF, EMTF, uGMT, uGMT_emu };

  double pt, eta, phi;
  int disc;
  int   bx, q, hits, charge, refLayer;
  TYPE  type;
  int   iProcessor, position;

  L1PhaseIIObj();

  double modulo2PI (double phi) const{ 
    while (phi > 2*M_PI) phi -= 2*M_PI;
    while (phi < 0.) phi += 2*M_PI;
    return phi;
  }

  bool isValid() const { return type!=NONE && pt>0;}

  double ptValue() const { return pt; }
  double etaValue() const { return eta; }
  double phiValue() const { return phi;}
  int chargeValue() const { return charge;}

  ClassDef(L1PhaseIIObj,2)
};


std::ostream & operator<< (std::ostream &out, const L1PhaseIIObj &o);

#endif
