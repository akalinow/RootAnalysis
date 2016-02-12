#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif


Wtau::Wtau(){
  clear();
}

Wtau::~Wtau(){;}

void Wtau::clear(){

  pt_ = -999;

  leadingTk_ = TLorentzVector();

  sv_ = TVector3();
  nPCA_ = TVector3();
  nPCAAODvx_ = TVector3();
  nPCAGenvx_ = TVector3();
  nPCARefitvx_ = TVector3();

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

void Wevent::clear(){

#ifdef PROJECT_NAME
  decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decModeMinus_  = 99;
  decModePlus_  = 99;
#endif

  thePV_ = TVector3();
  pfPV_ = TVector3();
  refitPfPV_ = TVector3();
  refitPfPVNoBS_ = TVector3();

  isRefit_ = false;
  nTracksInRefit_ = 0;

  run_ = 0;
  lumi_ = 0;
  event_ = 0;
  nup_ = -999, npu_ = -999, npv_ = -999; 
  paircount_ = -1;
  genevtweight_ = 1;
  sample_ = -1;
  bosonId_ = -1;
#ifdef PROJECT_NAME
  decayModeBoson_ = WawGenInfoHelper::bosonDecayModes::kUndefined;
#else
  decayModeBoson_ = -1;
#endif

}

Wpair::Wpair(){;}
Wpair::~Wpair(){;}
void Wpair::clear(){}
