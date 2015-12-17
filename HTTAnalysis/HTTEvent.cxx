#ifdef PROJECT_NAME
#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"
#include "WarsawAnalysis/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

//////////////////////////////////////////////
//////////////////////////////////////////////
HTTEvent::HTTEvent(){

  clear();
  
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

}
//////////////////////////////////////////////
//////////////////////////////////////////////
void DiTauData::clear(){

#ifdef PROJECT_NAME
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
Wtau::Wtau(){
  clear();
}
Wtau::~Wtau(){;}

void Wtau::clear(){

  pt_ = -999;
}

Wmu::Wmu(){
  clear();
}
Wmu::~Wmu(){;}

void Wmu::clear(){

  pt_ = -999;
}

Welectron::Welectron(){
  clear();
}
Welectron::~Welectron(){;}

void Welectron::clear(){

  pt_ = -999;
}

Wevent::Wevent(){run_ = 0; lumi_ = 0; event_ = 0;}
Wevent::~Wevent(){;}

Wpair::Wpair(){;}
Wpair::~Wpair(){;}

Wtriggers::Wtriggers(){;}
Wtriggers::~Wtriggers(){;}
