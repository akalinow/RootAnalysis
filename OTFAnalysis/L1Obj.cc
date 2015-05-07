/*
 * L1Obj.cc
 *
 *  Created on: 13 lut 2014
 *      Author: akalinow
 */

#include <ostream>

#include "L1Obj.h"


std::ostream & operator<< (std::ostream &out, const L1Obj &o)
{
  out<<"L1Obj: ";
  out <<" pt: "<<o.pt<<", eta: "<<o.eta<<", phi: "<<o.phi
      <<" charge: "<<o.charge
      <<", q: "<<o.q<<", bx: "<<o.bx;
  return out;
}


