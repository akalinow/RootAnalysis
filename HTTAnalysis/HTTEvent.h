#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TBits.h"
#include <map>
#include <vector>
#include <bitset> 

#include "PropertyEnum.h"
#include "TriggerEnum.h"
#include "SelectionBitsEnum.h"
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTEvent{

 public:

  enum sampleTypeEnum {DUMMY = -1, DATA = 0, DY = 1, DYLowM = 2, WJets=3, TTbar=4, h=5, H=6, A=7};

  ///Copy from LLRHiggsTauTau/NtupleProducer/plugins/HTauTauNtuplizer.cc
  static const int ntauIds = 30;
  TString tauIDStrings[ntauIds] = {
    "byLoosePileupWeightedIsolation3Hits",
    "byMediumPileupWeightedIsolation3Hits",
    "byTightPileupWeightedIsolation3Hits",
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    "againstMuonLoose3",
    "againstMuonTight3",
    "againstElectronVLooseMVA6",
    "againstElectronLooseMVA6",
    "againstElectronMediumMVA6",
    "againstElectronTightMVA6",
    "againstElectronVTightMVA6",
    "byVLooseIsolationMVArun2v1DBoldDMwLT",
    "byLooseIsolationMVArun2v1DBoldDMwLT",
    "byMediumIsolationMVArun2v1DBoldDMwLT",
    "byTightIsolationMVArun2v1DBoldDMwLT",
    "byVTightIsolationMVArun2v1DBoldDMwLT",
    "byVLooseIsolationMVArun2v1DBnewDMwLT",
    "byLooseIsolationMVArun2v1DBnewDMwLT",
    "byMediumIsolationMVArun2v1DBnewDMwLT",
    "byTightIsolationMVArun2v1DBnewDMwLT",
    "byVTightIsolationMVArun2v1DBnewDMwLT",
    "byLooseCombinedIsolationDeltaBetaCorr3HitsdR03",
    "byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",
    "byTightCombinedIsolationDeltaBetaCorr3HitsdR03",
    "byLooseIsolationMVArun2v1DBdR03oldDMwLT",
    "byMediumIsolationMVArun2v1DBdR03oldDMwLT",
    "byTightIsolationMVArun2v1DBdR03oldDMwLT",
    "byVTightIsolationMVArun2v1DBdR03oldDMwLT"
  };

  HTTEvent(){ clear();}
    
  ~HTTEvent(){}

  ///Data member setters.
  void setRun(unsigned int x){runId = x;}
  
  void setEvent(unsigned long int x){eventId = x;}
  
  void setNPU(float x){nPU = x;}
  
  void setNPV(unsigned int x){nPV = x;}
  
  void setMCatNLOWeight(float x){aMCatNLOweight = x;}

  void setMCWeight(float x){mcWeight = x;}

  void setPtReWeight(float x){ptReWeight = x;}

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

  void setSelectionBit(SelectionBitsEnum iBit, bool value = true) {selectionWord.SetBitNumber((int)iBit, value);}

  void setMET(const TVector2 &aVector) {met = aVector;}
  ////////////////////////

  ///Reset class data members
  void clear();

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}
    
  float getNPU() const {return nPU;}
  
  unsigned int getNPV() const {return nPV;}
  
  float getMCatNLOWeight() const {return aMCatNLOweight;}

  float setPtReWeight() const {return ptReWeight;}

  float getMCWeight() const {return mcWeight;}

  float getLHE_Ht() const {return lheHt;}

  int getLHEnOutPartons() const {return lheNOutPartons;}
  
  sampleTypeEnum getSampleType() const {return sampleType;}
  
  int getDecayModeMinus() const {return decayModeMinus;}
  
  int getDecayModePlus() const {return decayModePlus;}
  
  int getDecayModeBoson() const {return decayModeBoson;}

  TVector2 getMET() const {return met;}
  
  const TVector3 & getGenPV() const {return genPV;}

  const TVector3 & getAODPV() const {return AODPV;}

  const TVector3 & getRefittedPV() const {return refittedPV;}

  bool getIsRefit() const {return isRefit;};

  int getNTracksInRefit() const {return nTracksInRefit;}

  bool checkSelectionBit(SelectionBitsEnum iBit) const {return selectionWord.TestBitNumber((unsigned int)iBit);}

 private:

  ///Event run and number
  unsigned int runId;
  unsigned long int eventId;

  //Generator event weight
  float mcWeight;

  ///Weight used to modify the pt shape.
  float ptReWeight;
  
  ///Ht value from LHE record.
  float lheHt;

  ///Number of outgoing partons from LHE record
  int   lheNOutPartons;

  ///MCatNLO weight
  float aMCatNLOweight;

  ///Number of true PU vertices from MC
  float nPU;

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

  ///Bit word coding event selection result
  TBits selectionWord;
  
  //MET vector
  TVector2 met;

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

  void setChargedP4(const TLorentzVector &aP4) { chargedP4 = aP4;}
  
  void setNeutralP4(const TLorentzVector &aP4) { neutralP4 = aP4;}

  void setP4ScaleUp(const TLorentzVector &aP4) { p4ScaleUp = aP4;}

  void setP4ScaleDown(const TLorentzVector &aP4) { p4ScaleDown = aP4;}

  void setPCA(const TVector3 &aV3) {pca = aV3;}

  void setPCARefitPV(const TVector3 &aV3) {pcaRefitPV = aV3;}
  
  void setPCAGenPV(const TVector3 &aV3) {pcaGenPV = aV3;}

  void setProperties(const std::vector<float> & aProperties) { properties = aProperties;}

  ///Data member getters.
  TLorentzVector getP4() const {return p4;}

  TLorentzVector getChargedP4() const {return chargedP4;}

  TLorentzVector getNeutralP4() const {return neutralP4;}

  TLorentzVector getP4ScaleUp() const {return p4ScaleUp;}

  TLorentzVector getP4ScaleDown() const {return p4ScaleDown;}

  TVector3 getPCA() const {return pca;}

  TVector3 getPCARefitPV() const {return pcaRefitPV;}

  TVector3 getPCAGenPV() const {return pcaGenPV;}

  int getPDGid() const {return getProperty(PropertyEnum::PDGId);}

  int getCharge() const {return getProperty(PropertyEnum::charge);}

  float getProperty(PropertyEnum index) const {return (unsigned int)index<properties.size()?  properties[(unsigned int)index]: -999;}

  bool hasTriggerMatch(TriggerEnum index) const {return (unsigned int)getProperty(PropertyEnum::isGoodTriggerType)& (1<<(unsigned int)index) &&
                                                        (unsigned int)getProperty(PropertyEnum::FilterFired)& (1<<(unsigned int)index);}
 private:

  ///Nominal particle four-momentum
  TLorentzVector p4;

  ///Charged and neuttral components four-momentum
  TLorentzVector chargedP4, neutralP4;

  ///Particle four-momentum after scaling up/down.
  ///Variated scale depends on the particle type (mu, tau, e)
  TLorentzVector p4ScaleUp, p4ScaleDown;

  ///Vectors from primary vertex to point of closest approach (PCA)
  ///calculated with respect to AOD vertex, refitted and generated vertex.
  TVector3 pca, pcaRefitPV, pcaGenPV;

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

  void setMuonTriggerSF(float aSF){muonTriggerSF = aSF;}
  
  void setMETMatrix(float m00, float m01, float m10, float m11) {metMatrix.push_back(m00); metMatrix.push_back(m01); metMatrix.push_back(m10); metMatrix.push_back(m11);}
  
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

  float getMuonTriggerSF() const {return muonTriggerSF;}
  
  std::vector<float> getMETMatrix() {return metMatrix;}

 private:

  ///Nominal pair p4 (sum of legs p4)
  TLorentzVector p4;

  ///P4 calculated with SV fit. Nominal, and with scale up/down
  TLorentzVector p4SVFit, p4SVFitScaleUp, p4SVFitScaleDown;

  ///MET vectors
  TVector2 met, metSVfit;
  
  //MVAMET covariance matrix in order 00,01,10,11
  std::vector<float> metMatrix;

  ///MT calculated for (leg1,MET)
  float mtLeg1;

  ///MT calculated for (leg2,MET)
  float mtLeg2;
    
  ///Pair legs
  HTTParticle leg1, leg2;

  ///Lepton selection scale factor
  float muonTriggerSF;
  ///Tau selection scale factor

};

#endif

