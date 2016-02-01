/*
 * L1Obj.cc
 *
 *  Created on: 13 lut 2014
 *      Author: akalinow
 */

#include <ostream>
#include <bitset>

#include "TMath.h"
#include "L1Obj.h"


std::ostream & operator<< (std::ostream &out, const L1Obj &o)
{

  std::bitset<17> bits(o.q);

  out<<"L1Obj: ";
  out <<" pt: "<<o.pt<<", eta: "<<o.eta<<", phi: "<<o.phi
      <<" charge: "<<o.charge
      <<" refLayer: "<<o.refLayer
      <<" disc: "<<o.disc
      <<", q: "<<bits.to_string()<<", bx: "<<o.bx;
  return out;
}


