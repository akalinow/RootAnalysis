#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"

//
// class declaration
//

class DiTauData{

 public:

  ///Data common for generator and reconstruction levels.
  Int_t decModeMinus_, decModePlus_;
  TVector3 thePV_;
  TVector3 svMinus_, svPlus_;
  TVector3 nPiPlus_, nPiMinus_;

  TVector3 nPiPlusAODvx_,  nPiPlusGenvx_,  nPiPlusRefitvx_;
  TVector3 nPiMinusAODvx_,  nPiMinusGenvx_,  nPiMinusRefitvx_; 
  
  TLorentzVector p4Sum_;
  TLorentzVector met_, metNu_;
  TLorentzVector piMinus_, piPlus_;

  TLorentzVector neutralMinus_, neutralPlus_;
  TLorentzVector piZeroMinus_, piZeroPlus_;
  
  TLorentzVector tauMinus_, tauPlus_;
  TLorentzVector visTauMinus_, visTauPlus_;

  Float_t phi_, phi2_;
  Float_t rho_, phiRho_;

  Float_t         yMinus_, yPlus_, yMinus2_, yPlus2_;
  Float_t         yMinusLab_, yPlusLab_, yMinusLab2_, yPlusLab2_;
  Float_t         isoMinus_, isoPlus_, outerMinus_, outerPlus_;

  Int_t nChargedMinus_, nChargedPlus_; 
  Int_t nPiZeroMinus_,  nPiZeroPlus_;
  Int_t nNeutralMinus_, nNeutralPlus_;
    
  ///Reconstruction level only data
  TVector3 pfPV_, pt2PV_, refitPfPV_, refitPfPVNoBS_;
  float zErrorRefitPV_;

  bool isRefit_;
  int nTracksInRefit_;
  int pfPVIndex_, pt2PVindex_;

  ///Reset data members to default values.
  void clear();

  ///Pointers
  Int_t *decModeMinusPtr_, *decModePlusPtr_;
  TVector3 *thePVPtr_;
  TVector3 *svMinusPtr_, *svPlusPtr_;
  TVector3 *nPiPlusPtr_, *nPiMinusPtr_;
  
  TLorentzVector *p4SumPtr_;
  TLorentzVector *metPtr_, *metNuPtr_;
  TLorentzVector *piMinusPtr_, *piPlusPtr_;

  TLorentzVector *neutralMinusPtr_, *neutralPlusPtr_;
  TLorentzVector *piZeroMinusPtr_, *piZeroPlusPtr_;
  
  TLorentzVector *tauMinusPtr_, *tauPlusPtr_;
  TLorentzVector *visTauMinusPtr_, *visTauPlusPtr_;

  Float_t *phiPtr_, *phi2Ptr_;
  Float_t *rhoPtr_, *phiRhoPtr_;

  Float_t *yMinusPtr_, *yPlusPtr_, *yMinus2Ptr_, *yPlus2Ptr_;
  Float_t *yMinusLabPtr_, *yPlusLabPtr_, *yMinusLab2Ptr_, *yPlusLab2Ptr_;
  Float_t *isoMinusPtr_, *isoPlusPtr_, *outerMinusPtr_, *outerPlusPtr_;

  Int_t *nChargedMinusPtr_, *nChargedPlusPtr_; 
  Int_t *nPiZeroMinusPtr_,  *nPiZeroPlusPtr_;
  Int_t *nNeutralMinusPtr_, *nNeutralPlusPtr_;
  
  void setPointers();
};


class HTTEvent {

 public:

  HTTEvent();

  ~HTTEvent();

  ///Reset data members to default values.
  void clear();
  
  float run_, lumi_, event_;
  Int_t bosonId_;
  Int_t *bosonIdPtr_;

  DiTauData genEvent_, recoEvent_, toyEvent_;
    
};

#endif

