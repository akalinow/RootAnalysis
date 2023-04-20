#ifndef MuonObj_H
#define MuonObj_H

#include <ostream>
#include "TObject.h"
#include <cmath>
class MuonObj : public TObject {
public:
 MuonObj(Float_t pt=0., Float_t eta=0., Float_t phi=0.,Int_t charge=0, UInt_t nrpchits=0., UInt_t ndthits=0., 
   UInt_t ncschits=0., UInt_t ntrackerhits=0., UInt_t nmatchedstations=0., Bool_t mediumID=false, Bool_t tightID=false):
  thePt(pt),theEta(eta),thePhi(phi),theCharge(charge),nRPCHits(nrpchits),nDTHits(ndthits),nCSCHits(ncschits),nTrackerHits(ntrackerhits),nMatchedStations(nmatchedstations),isMedium(mediumID),isTight(tightID){}
  virtual ~MuonObj() {}
public:
  Float_t pt() const { return thePt;}
  Float_t eta() const { return theEta;}
  Float_t phi() const { return thePhi;}
  Int_t charge() const { return theCharge;}
  UInt_t nrpchits() const { return nRPCHits;}
  UInt_t ndthits() const { return nDTHits;}
  UInt_t ncschits() const { return nCSCHits;}
  UInt_t ntrackerhits() const { return nTrackerHits;}
  UInt_t nmatchedstations() const { return nMatchedStations;}
  Bool_t mediumID() const { return isMedium;}
  Bool_t tightID() const { return isTight;}


private:
  Float_t thePt,theEta,thePhi;
  Int_t theCharge;
  UInt_t nRPCHits, nDTHits, nCSCHits, nTrackerHits, nMatchedStations;
  Bool_t isMedium, isTight;

public:
  ClassDef(MuonObj,6);
};

std::ostream & operator<< (std::ostream &out, const MuonObj &o);

#endif


