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

  genPV*=0;
  AODPV*=0;
  refittedPV*=0;

  isRefit = false;

  nTracksInRefit = 0;

  selectionWord.ResetAllBits();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTParticle::clear(){

  p4*=0;
  chargedP4*=0;
  neutralP4*=0;

  pca*=0;
  pcaRefitPV*=0;
  pcaGenPV*=0;

  properties.clear();

  lastSystEffect = HTTAnalysis::NOMINAL;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTParticle::getSystScaleP4(HTTAnalysis::sysEffects type) const{

  if(type==HTTAnalysis::NOMINAL) {
    lastSystEffect = type;
    return p4;
  }
  else if(lastSystEffect==type) return p4Cache;

  lastSystEffect = type;

  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==5){
    ///True taus
    if(type!=HTTAnalysis::TESUp && type!=HTTAnalysis::TESDown) return p4;
    float TES = 0.03;
    if(type==HTTAnalysis::TESDown) TES*=-1;
    return getShiftedP4(1+TES);
  }
  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==3){
    ///Fake e->tau
    if(type!=HTTAnalysis::E2TUp && type!=HTTAnalysis::E2TDown) return p4;
    float EES = 0.03;
    if(type==HTTAnalysis::E2TDown) EES*=-1;
    return getShiftedP4(1+EES);
  }
  if(abs(getPDGid())==15 && getProperty(PropertyEnum::mc_match)==4){
    ///Fake mu->tau
    if(type!=HTTAnalysis::M2TUp && type!=HTTAnalysis::M2TDown) return p4;
    float MES = 0.03;
    if(type==HTTAnalysis::M2TDown) MES*=-1;
    return getShiftedP4(1+MES);
  }
  if(abs(getPDGid())==98){
    if(type!=HTTAnalysis::JESUp && type!=HTTAnalysis::JESDown) return p4;
    float JES = getProperty(PropertyEnum::jecUnc);
    if(type==HTTAnalysis::JESDown) JES*=-1;
    return getShiftedP4(1+JES);
  }

  p4Cache = p4;
  return p4;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTParticle::getShiftedP4(float scale) const{

  double pt = p4.Perp();
  double energy =  p4.E();
  pt*=scale;
  double shiftedMomentum = pt/sin(p4.Theta());
  energy = sqrt(p4.M2() + pow(shiftedMomentum,2));
  p4Cache = p4;
  p4Cache.SetRho(shiftedMomentum);
  p4Cache.SetE(energy);
  return p4Cache;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTPair::clear(){

  if(!p4Vector.size()) p4Vector.resize(HTTAnalysis::DUMMY_SYS);
  if(!svMetVector.size()) svMetVector.resize(HTTAnalysis::DUMMY_SYS);

  for(auto &it:p4Vector) it*=0;
  for(auto &it:svMetVector) it*=0;

  metMatrix.clear();

  mtLeg1= -999;
  mtLeg2 = -999;

  leg1.clear();
  leg2.clear();

  lastSystEffect = HTTAnalysis::NOMINAL;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTPair::getP4(HTTAnalysis::sysEffects type) const {

  //std::cout<<"type: "<<(unsigned int)type<<" size: "<<p4Vector.size()<<std::endl;

  if(p4Vector.size()>(unsigned int)type) return p4Vector[(unsigned int)type];
  return p4Vector[(unsigned int)HTTAnalysis::NOMINAL];
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TVector2 & HTTPair::getSystScaleMET(HTTAnalysis::sysEffects type) const{

  if(type==HTTAnalysis::NOMINAL ||
  (unsigned int)type>(unsigned int)HTTAnalysis::DUMMY_SYS) {
    lastSystEffect = type;
    return met;
  }
  else if(lastSystEffect==type) return metCache;

  double metX = met.X();
  metX+=leg1.getP4(HTTAnalysis::NOMINAL).X();
  metX+=leg2.getP4(HTTAnalysis::NOMINAL).X();
  metX-=leg1.getP4(type).X();
  metX-=leg2.getP4(type).X();

  double metY = met.Y();
  metY+=leg1.getP4(HTTAnalysis::NOMINAL).Y();
  metY+=leg2.getP4(HTTAnalysis::NOMINAL).Y();
  metY-=leg1.getP4(type).Y();
  metY-=leg2.getP4(type).Y();

  metCache.SetX(met.X());
  metCache.SetY(met.Y());
  return metCache;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
float HTTPair::getSystScaleMT(const HTTParticle &aParticle,
			      HTTAnalysis::sysEffects type) const{

  const TVector2 & metScaled = getSystScaleMET(type);
  const TLorentzVector & legP4 = aParticle.getP4(type);
  float sumP2 = pow(metScaled.X() + legP4.X(),2) +
                pow(metScaled.Y() + legP4.Y(),2);
  float sumE2 = pow(metScaled.Mod() + legP4.Perp(),2);
  float mT = sqrt(sumE2 - sumP2);
  return mT;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
