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

	 UInt_t fUniqueID;
	 UInt_t fBits;
     Float_t pt;
     Float_t eta;
     Float_t phi;
     Int_t bx;
     Int_t q;
     Int_t charge;
     Int_t type;

 };

#endif /* L1OBJ_H_ */
