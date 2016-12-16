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
  ptReWeight = 1.0;
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

  pca = TVector3();
  pcaRefitPV = TVector3();
  pcaGenPV = TVector3();

  properties.clear();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TLorentzVector HTTParticle::getSystScaleP4(sysEffects::sysEffectsEnum type) const{

  if(type==sysEffects::NOMINAL || type==sysEffects::NOMINAL_SVFIT) return p4;

  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==5){
    ///True taus
    if(type!=sysEffects::TESUp && type!=sysEffects::TESDown) return p4;
    float TES = 0.03;
    if(type==sysEffects::TESDown) TES*=-1;
    return getShiftedP4(1+TES);
  }
  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==3){
    ///Fake e->tau
    if(type!=sysEffects::E2TUp && type!=sysEffects::E2TDown) return p4;
    float EES = 0.03;
    if(type==sysEffects::E2TDown) EES*=-1;
    return getShiftedP4(1+EES);
  }
  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==4){
    ///Fake mu->tau
    if(type!=sysEffects::M2TUp && type!=sysEffects::M2TDown) return p4;
    float MES = 0.03;
    if(type==sysEffects::M2TDown) MES*=-1;
    return getShiftedP4(1+MES);
  }
  if(abs(getPDGid())==98){
    if(type!=sysEffects::JESUp && type!=sysEffects::JESDown) return p4;
    float JES = getProperty(PropertyEnum::jecUnc);
    if(type==sysEffects::JESDown) JES*=-1;
    return getShiftedP4(1+JES); 
  }
  return p4;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TLorentzVector HTTParticle::getShiftedP4(float scale) const{

  TVector3 momentum = p4.Vect();
  float pt = momentum.Perp();
  float energy =  p4.E();
  pt*=scale;
  float shiftedMomentum = pt/sin(momentum.Theta());
  momentum = shiftedMomentum*momentum.Unit();
  energy = sqrt(momentum.Mag() + p4.M2());
  TLorentzVector scaled(momentum,energy);
  return scaled;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTPair::clear(){

  p4Vector.clear();
  svMetVector.clear();

  p4Vector.resize((unsigned int)sysEffects::DUMMY);
  svMetVector.resize((unsigned int)sysEffects::DUMMY);

  metMatrix.clear();

  mtLeg1= -999;
  mtLeg2 = -999;

  leg1 = HTTParticle();
  leg2 = HTTParticle();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
TVector2 HTTPair::getSystScaleMET(sysEffects::sysEffectsEnum type) const{

  if(type==sysEffects::NOMINAL) return met;

  TLorentzVector metShiftedP4(met.X(), met.Y(), 0, met.Mod());

  metShiftedP4+=leg1.getP4(sysEffects::NOMINAL);
  metShiftedP4+=leg2.getP4(sysEffects::NOMINAL);

  metShiftedP4-=leg1.getP4(type);
  metShiftedP4-=leg2.getP4(type);

  TVector2 metShifted(metShiftedP4.X(), metShiftedP4.Y());
  return metShifted;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
float HTTPair::getSystScaleMT(const HTTParticle &aParticle,
			      sysEffects::sysEffectsEnum type) const{

  TVector2 metScaled = getSystScaleMET(type);
  TLorentzVector metP4(metScaled.X(), metScaled.Y(),0, metScaled.Mod());
  TLorentzVector legP4 = aParticle.getP4(type);
  legP4.SetZ(0);
  legP4.SetE(legP4.Perp());

  float mT = (metP4+legP4).M();
  return mT;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
