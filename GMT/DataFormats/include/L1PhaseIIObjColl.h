#ifndef L1PhaseIIObjColl_H
#define L1PhaseIIObjColl_H

#include "L1PhaseIIObj.h"
#include <vector>
#include <cmath>

class L1PhaseIIObjColl : public TObject {

public:
  L1PhaseIIObjColl() {}
  virtual ~L1PhaseIIObjColl(){}

  typedef L1PhaseIIObj::TYPE TYPE;
  void set(const std::vector<L1PhaseIIObj> & obj) { theL1PhaseIIObj = obj; }
  void set(const std::vector<bool> & comp) { theL1Matching = comp; }
  void set(const std::vector<double> & dr) { theDeltaR = dr; }
  void push_back(const L1PhaseIIObj & obj, bool match, double deltaR);

  const std::vector<L1PhaseIIObj> & getL1PhaseIIObjs() const { return theL1PhaseIIObj; }
  operator  const std::vector<L1PhaseIIObj> & () const {  return theL1PhaseIIObj; }

  const std::vector<bool> & getL1PhaseIIObjsMatching() const { return theL1Matching; }
  const std::vector<double> & getL1PhaseIIObjDeltaR() const { return theDeltaR; }
 
  L1PhaseIIObjColl l1RpcColl() const {
    return selectByType(L1PhaseIIObj::RPCb)+selectByType(L1PhaseIIObj::RPCf);
  }
  L1PhaseIIObjColl l1RpcCollEmu() const {
    return selectByType(L1PhaseIIObj::RPCb_emu)+selectByType(L1PhaseIIObj::RPCf_emu);
  }
  L1PhaseIIObjColl l1OthColl() const {
    return selectByType(L1PhaseIIObj::DT)+selectByType(L1PhaseIIObj::CSC);
  }

  L1PhaseIIObjColl selectByType( TYPE t1) const;
  L1PhaseIIObjColl selectByPt( double ptMin = 0., double ptMax = 161.) const;
  L1PhaseIIObjColl selectByPtMin( double ptMin = 0.) const;
  L1PhaseIIObjColl selectByEta( double etaMin = -1.61, double etaMax = 1.61) const;
  L1PhaseIIObjColl selectByBx(  int bxMin = 0, int bxMax = 0) const;
  L1PhaseIIObjColl selectByQuality( int qMin = 0, int qMax = 16) const;
  L1PhaseIIObjColl selectByMatched() const;
  L1PhaseIIObjColl selectByDeltaR( double deltaRMax) const;
  L1PhaseIIObjColl operator+(const L1PhaseIIObjColl &o) const;
  operator bool() const {return !theL1PhaseIIObj.empty(); } 
  bool operator!() const {return theL1PhaseIIObj.empty(); }

  // inline bool isMatching_DRBx_At(double deltaR, int bx, double ptMin, double eta, double phi) const {
  //   for (unsigned int i=0; i< theL1PhaseIIObj.size(); i++) if ( (bx ==  theL1PhaseIIObj[i].bx) && ( reco::deltaR(theL1PhaseIIObj[i].eta, theL1PhaseIIObj[i].phi, eta, phi) < deltaR) && (theL1PhaseIIObj[i].pt >= ptMin) ) return true;
  //   return false;
  // }
  inline bool isMatching_DRBx(double deltaR, int bx, double ptMin = 0.) const {
    for (unsigned int i=0; i< theL1PhaseIIObj.size(); i++) if ( (bx ==  theL1PhaseIIObj[i].bx) && ( theDeltaR[i] < deltaR) && (theL1PhaseIIObj[i].pt >= ptMin) ) return true;
    return false;
  }

  inline bool isMatching_PtminPtmaxBx(double ptMin, double ptMax, int bx, bool firstBXonly) const {
    bool firstBX = true; 
    bool result = false;
    for (unsigned int i=0; i< theL1PhaseIIObj.size(); i++) {
      if (theL1PhaseIIObj[i].bx < bx && theL1PhaseIIObj[i].pt >= ptMin) firstBX = false;
      if ( (theL1PhaseIIObj[i].pt >= ptMin && theL1PhaseIIObj[i].pt < ptMax) && (bx == theL1PhaseIIObj[i].bx) ) result = true;
    }
    return firstBXonly ? (result && firstBX) : result;
  }

//tmp
  std::vector<L1PhaseIIObj> getL1PhaseIIObjsMatched(double ptMin = 0) const;
  std::vector<L1PhaseIIObj> getL1PhaseIIObjsSelected(
		  bool requireMatched = true, bool requireNonMatched = false, 
		  double ptMin = 0., double ptMax = 161.,
		  int bxMin = 0, int bxMax = 0,                  
		  double etaMin = -2.5, double etaMax = 2.5,
		  double phiMin = 0., double phiMax = 7.,
		  int qMin = 0, int qMax = 99) const;
  static  std::vector<L1PhaseIIObj> typeSelector(const  std::vector<L1PhaseIIObj> & col,  TYPE t1=L1PhaseIIObj::NONE, TYPE t2=L1PhaseIIObj::NONE, TYPE t3=L1PhaseIIObj::NONE, TYPE t4=L1PhaseIIObj::NONE);
  
private:
  std::vector<L1PhaseIIObj> theL1PhaseIIObj;
  std::vector<bool> theL1Matching;
  std::vector<double> theDeltaR;

public:
  // ClassDef(L1PhaseIIObjColl,1)

friend std::ostream & operator<< (std::ostream &out, const L1PhaseIIObjColl&s);

};



#endif
