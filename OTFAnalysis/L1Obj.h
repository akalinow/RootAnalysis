/*
 * L1Obj.h
 *
 *  Created on: 13 lut 2014
 *      Author: akalinow
 */

#ifndef L1OBJ_H_
#define L1OBJ_H_

#include "TObject.h"


 struct L1Obj{

 L1Obj() : pt(-1.),eta(99.),phi(99.),disc(-999), bx(0),q(-1), hits(0), charge(99), type(0) {}

   UInt_t fUniqueID;
   UInt_t fBits;
   Float_t pt;
   Float_t eta;
   Float_t phi;
   Float_t disc;
   Int_t bx;
   Int_t q;
   Int_t hits;
   Int_t charge;
   Int_t type;
   Int_t refLayer;
   
 };

std::ostream & operator<< (std::ostream &out, const L1Obj &o);

#endif /* L1OBJ_H_ */
