#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTEvent::clear(){

  runId = 0;
  eventId = 0;

  nPU = 0;
  nPV = 0;

  mcWeight = 1.0;
  lheHt = 1.0;
  lheNOutPartons = 0;
  aMCatNLOweight = 1.0;
  
  sampleType = DUMMY;

#ifdef PROJECT_NAME
  decayModeMinus = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decayModePlus = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decayModeMinus  = 99;
  decayModePlus  = 99;
#endif


#ifdef PROJECT_NAME
  decayModeBoson = WawGenInfoHelper::bosonDecayModes::kUndefined;
#else
  decayModeBoson = -1;
#endif

  genPV =  TVector3();
  AODPV =  TVector3();
  refittedPV =  TVector3();

  isRefit = false;

  nTracksInRefit = 0;
  
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTParticle::clear(){

  p4 = TLorentzVector();

  chargedP4 = TLorentzVector();
  neutralP4 = TLorentzVector();

  p4ScaleUp = TLorentzVector();
  p4ScaleDown = TLorentzVector();

  pca = TVector3();
  pcaRefitPV = TVector3();
  pcaGenPV = TVector3();

  properties.clear();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTPair::clear(){

  p4 = TLorentzVector();

  p4SVFit = TLorentzVector();  
  p4SVFitScaleUp = TLorentzVector();  
  p4SVFitScaleDown = TLorentzVector();  

  met =  TVector2();
  metSVfit = TVector2();

  leg1 = HTTParticle();
  leg2 = HTTParticle();
  
}
////////////////////////////////////////////////
////////////////////////////////////////////////
