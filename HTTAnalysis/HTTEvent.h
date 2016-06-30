#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include <vector>
#include <bitset> 

#include "PropertyEnum.h"
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTEvent{

 public:

  enum sampleTypeEnum {DUMMY = -1, DATA = 0, DY = 1, WJets=2, TTbar=3, h=3, H=4, A=5, QCD=6};

  HTTEvent(){ clear();}
    
  ~HTTEvent(){}

  ///Data member setters.
  void setRun(unsigned int x){runId = x;}
  
  void setEvent(unsigned long int x){eventId = x;}
  
  void setNPU(unsigned int x){nPU = x;}
  
  void setNPV(unsigned int x){nPV = x;}
  
  void setMCatNLOWeight(float x){aMCatNLOweight = x;}

  void setMCWeight(float x){mcWeight = x;}

  void setLHE_Ht(float x){lheHt = x;}

  void setLHEnOutPartons(int x){lheNOutPartons = x;}

  void setSampleType(sampleTypeEnum x){sampleType = x;}
  
  void setDecayModeMinus(int x){decayModeMinus = x;}
  
  void setDecayModePlus(int x){decayModePlus = x;}
  
  void setDecayModeBoson(int x){decayModeBoson = x;}

  void setGenPV(const TVector3 & aPV) {genPV = aPV;}

  void setAODPV(const TVector3 & aPV) {AODPV = aPV;}

  void setRefittedPV(const TVector3 & aPV) {refittedPV = aPV;}

  void setIsRefit(bool aBit){isRefit = aBit;};

  void setNTracksInRefit(const int & nTracks) {nTracksInRefit = nTracks;};
  ////////////////////////

  ///Reset class data members
  void clear();

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}
    
  unsigned int getNPU() const {return nPU;}
  
  unsigned int getNPV() const {return nPV;}
  
  float getMCatNLOWeight() const {return aMCatNLOweight;}

  float getMCWeight() const {return mcWeight;}

  float getLHE_Ht() const {return lheHt;}

  int getLHEnOutPartons() const {return lheNOutPartons;}
  
  sampleTypeEnum getSampleType() const {return sampleType;}
  
  int getDecayModeMinus() const {return decayModeMinus;}
  
  int getDecayModePlus() const {return decayModePlus;}
  
  int getDecayModeBoson() const {return decayModeBoson;}
  
  const TVector3 & getGenPV() const {return genPV;}

  const TVector3 & getAODPV() const {return AODPV;}

  const TVector3 & getRefittedPV() const {return refittedPV;}

  bool getIsRefit() const {return isRefit;};

  int getNTracksInRefit() const {return nTracksInRefit;}

 private:

  ///Event run and number
  unsigned int runId;
  unsigned long int eventId;

  //Generator event weight
  float mcWeight;

  ///Ht value from LHE record.
  float lheHt;

  ///Number of outgoing partons from LHE record
  int   lheNOutPartons;

  ///MCatNLO weight
  float aMCatNLOweight;

  ///Number of PU vertices from MC
  unsigned int nPU;

  //Number of reocnstructed PV
  unsigned int nPV;

  ///Type of the physics process or DATA
  sampleTypeEnum sampleType;

  ///Boson (H, Z, W) decay mode
  int decayModeBoson;
    
  ///Tau decay modes
  int decayModeMinus, decayModePlus;

  ///Primary Vertices recontructed with different methods  
  //Generated PV position 
  TVector3 genPV;

  //PV stored in miniAOD
  TVector3 AODPV;

  ///PV recontructed from PF candidates, refitted
  TVector3 refittedPV;

  ///Flag marking if refit was successfull
  bool isRefit;

  ///Number of tracks used in the refit
  int nTracksInRefit;

};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTParticle{

  public:

  HTTParticle(){ clear();}
  
  ~HTTParticle(){}

  void clear();

  ///Data member setters.
  void setP4(const TLorentzVector &aP4) { p4 = aP4;}

  void setP4ScaleUp(const TLorentzVector &aP4) { p4ScaleUp = aP4;}

  void setP4ScaleDown(const TLorentzVector &aP4) { p4ScaleDown = aP4;}

  void setProperties(const std::vector<float> & aProperties) { properties = aProperties;}

  ///Data member getters.
  TLorentzVector getP4() const {return p4;}

  TLorentzVector getP4ScaleUp() const {return p4ScaleUp;}

  TLorentzVector getP4ScaleDown() const {return p4ScaleDown;}

  int getPDGid() const {return getProperty(PropertyEnum::PDGId);}

  int getCharge() const {return getProperty(PropertyEnum::charge);}

  float getProperty(PropertyEnum index) const {return (unsigned int)index<properties.size()?  properties[(unsigned int)index]: -999;}

 private:

  ///Nominal particle four-momentum
  TLorentzVector p4;

  ///Particle four-momentum after scaling up/down.
  ///Variated scale depends on the particle type (mu, tau, e)
  TLorentzVector p4ScaleUp, p4ScaleDown;

  ///Vector of vaious particle properties.
  ///Index generated automatically during conversion from
  ///LLR ntuple format
  std::vector<float> properties;

};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTPair{

 public:

  HTTPair(){ clear();}
  
  ~HTTPair(){}

  void clear();

  ///Data member setters.
  void setP4(const TLorentzVector &aP4) {p4 = aP4;}

  void setP4SVFit(const TLorentzVector &aP4) {p4SVFit = aP4;}

  void setP4SVFitScaleUp(const TLorentzVector &aP4) {p4SVFitScaleUp = aP4;}

  void setP4SVFitScaleDown(const TLorentzVector &aP4) {p4SVFitScaleDown = aP4;}

  void setMET(const TVector2 &aVector) {met = aVector;}

  void setMTLeg1(const float & aMT) {mtLeg1 = aMT;}

  void setMTLeg2(const float & aMT) {mtLeg2 = aMT;}

  void setMETSVfit(const TVector2 &aVector) {metSVfit = aVector;}
  
  void setLeg1(const HTTParticle &aParticle){leg1 = aParticle;}

  void setLeg2(const HTTParticle &aParticle){leg2 = aParticle;}


  ///Data member getters.
  TLorentzVector getP4() const {return p4;}

  TLorentzVector getP4SVFit() const {return p4SVFit;}

  TVector2 getMET() const {return met;}

  TVector2 getMETSVfit() const {return metSVfit;}

  float getMTLeg1() const {return mtLeg1;}

  float getMTLeg2() const {return mtLeg2;}

  TLorentzVector getP4SVFitScaleUp() {return p4SVFitScaleUp;}

  TLorentzVector getP4SVFitScaleDown() {return p4SVFitScaleDown;}

  HTTParticle getLeg1() const {return leg1;}

  HTTParticle getLeg2() const {return leg2;}

  HTTParticle getMuon() const {return abs(leg1.getPDGid())==13 ? leg1 : leg2; }

  HTTParticle getTau() const {return abs(leg1.getPDGid())==15 ? leg1 : leg2; }

  float getMTMuon() const {return abs(leg1.getPDGid())==13 ? getMTLeg1() : getMTLeg2(); }

 private:

  ///Nominal pair p4 (sum of legs p4)
  TLorentzVector p4;

  ///P4 calculated with SV fit. Nominal, and with scale up/down
  TLorentzVector p4SVFit, p4SVFitScaleUp, p4SVFitScaleDown;

  ///MET vectors
  TVector2 met, metSVfit;

  ///MT calculated for (leg1,MET)
  float mtLeg1;

  ///MT calculated for (leg2,MET)
  float mtLeg2;
    
  ///Pair legs
  HTTParticle leg1, leg2;

};

#endif

