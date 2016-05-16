#ifdef CMSSW
#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::HTTEvent(){

  clear();

  bosonIdPtr_ = &bosonId_;
  
  genEvent_.setPointers();
  recoEvent_.setPointers();
  toyEvent_.setPointers();
  
}
//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::~HTTEvent(){;}
//////////////////////////////////////////////
//////////////////////////////////////////////
void HTTEvent::clear(){
  
  run_ = -999;
  lumi_ = -999;
  event_ = -999;
  bosonId_ = 0;

  genEvent_.clear();
  recoEvent_.clear();
  toyEvent_.clear();

}
//////////////////////////////////////////////
//////////////////////////////////////////////
void DiTauData::setPointers(){

  decModeMinusPtr_ = &decModeMinus_;
  decModePlusPtr_ = &decModePlus_;
  
  thePVPtr_ = &thePV_;
  svMinusPtr_ = &svMinus_;
  svPlusPtr_ = &svPlus_;

  nPiPlusPtr_ = &nPiPlus_;
  nPiMinusPtr_ = &nPiMinus_;
  
  p4SumPtr_ = &p4Sum_;
  
  metPtr_ = &met_;
  metNuPtr_ = &metNu_;
  piMinusPtr_ = &piMinus_;
  piPlusPtr_ = &piPlus_;

  neutralMinusPtr_ = &neutralMinus_;
  neutralPlusPtr_ = &neutralPlus_;
  piZeroMinusPtr_ = &piZeroMinus_;
  piZeroPlusPtr_ = &piZeroPlus_;
  
  tauMinusPtr_ = &tauMinus_;
  tauPlusPtr_ = &tauPlus_;
  visTauMinusPtr_ = &visTauMinus_;
  visTauPlusPtr_ = &visTauPlus_;

  phiPtr_ = &phi_;
  phi2Ptr_ = &phi2_;

  rhoPtr_ = &rho_;
  phiRhoPtr_ = &phiRho_;

  yMinusPtr_ = &yMinus_;
  yPlusPtr_ = &yPlus_;
  yMinus2Ptr_ = &yMinus2_;
  yPlus2Ptr_ = &yPlus2_;
 
  yMinusLabPtr_ = &yMinusLab_;
  yPlusLabPtr_ = &yPlusLab_;
  yMinusLab2Ptr_ = &yMinusLab2_;
  yPlusLab2Ptr_ = &yPlusLab2_;

  isoMinusPtr_ = &isoMinus_;
  isoPlusPtr_ = &isoPlus_;
  outerMinusPtr_ = &outerMinus_;
  outerPlusPtr_ = &outerPlus_;

  nChargedMinusPtr_ = &nChargedMinus_;
  nChargedPlusPtr_ = &nChargedPlus_;
  
  nPiZeroMinusPtr_ = &nPiZeroMinus_;
  nPiZeroPlusPtr_ = &nPiZeroPlus_;

  nNeutralMinusPtr_ = &nNeutralMinus_;
  nNeutralPlusPtr_ = &nNeutralPlus_;

}
//////////////////////////////////////////////
//////////////////////////////////////////////
void DiTauData::clear(){

#ifdef CMSSW
  decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decModeMinus_  = 99;
  decModePlus_  = 99;
#endif

  thePV_ = TVector3();
  svMinus_ = TVector3();
  svPlus_ = TVector3();
  nPiPlus_ = TVector3();
  nPiMinus_ = TVector3();

  p4Sum_ = TLorentzVector();
  piMinus_ = TLorentzVector();
  piPlus_ = TLorentzVector();
  tauMinus_ = TLorentzVector();
  tauPlus_ = TLorentzVector();
  visTauMinus_ = TLorentzVector();
  visTauPlus_ = TLorentzVector();

  pfPV_ = TVector3();
  pt2PV_ = TVector3();
  refitPfPV_ = TVector3();
  refitPfPVNoBS_ = TVector3();

  isRefit_ = false;
  nTracksInRefit_ = 0;
  
  pfPVIndex_  = -1;
  pt2PVindex_ = -1;
   
}
//////////////////////////////////////////////
//////////////////////////////////////////////
