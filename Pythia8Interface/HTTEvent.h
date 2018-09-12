#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TBits.h"
#include <map>
#include <vector>
#include <bitset>
#include <iostream>

#include "PropertyEnum.h"
#include "TriggerEnum.h"
#include "SelectionBitsEnum.h"
#include "AnalysisEnums.h"

///////////////////////////////////////////////////
///////////////////////////////////////////////////
struct SVFitEvent{

  double tree_l1_type, tree_l1_decayMode, tree_l1_px, tree_l1_py, tree_l1_pz, tree_l1_pe, tree_l1_mass;
  double tree_l2_type, tree_l2_decayMode, tree_l2_px, tree_l2_py, tree_l2_pz, tree_l2_pe, tree_l2_mass;
  double tree_metX, tree_metY, tree_metCov00, tree_metCov11, tree_metCov10, tree_metCov01;
  double tree_mcMass, tree_cubaMass;
  double tree_cpuTime_MC, tree_cpuTime_Cuba;

};
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

  HTTEvent(){clear();}

  ~HTTEvent(){}

  ///Data member setters.
  void setRun(unsigned int x){runId = x;}

  void setEvent(unsigned long int x){eventId = x;}

  void setNPU(float x){nPU = x;}

  void setNPV(unsigned int x){nPV = x;}

  void setMCatNLOWeight(float x){aMCatNLOweight = x;}

  void setMCWeight(float x){mcWeight = x;}

  void setPtReWeight(float x){ptReWeight = x;}

  void setPtReWeightSUSY(float x){ptReWeightSUSY = x;}

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

  void setMETFilterDecision(unsigned int aMETFilterDecision) {metFilterDecision = aMETFilterDecision;}
  ////////////////////////

  ///Reset class data members
  void clear();

  void clearSelectionWord() {selectionWord.ResetAllBits();}

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}

  float getNPU() const {return nPU;}

  unsigned int getNPV() const {return nPV;}

  float getMCatNLOWeight() const {return aMCatNLOweight;}

  float getPtReWeight() const {return ptReWeight;}

  float getPtReWeightSUSY() const {return ptReWeightSUSY;}

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

  unsigned int getMETFilterDecision() const { return metFilterDecision;}

 private:

  ///Event run and number
  unsigned int runId;
  unsigned long int eventId;

  //Generator event weight
  float mcWeight;

  ///Weight used to modify the pt shape.
  float ptReWeight, ptReWeightSUSY;

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

  //MET filter decision
  unsigned int metFilterDecision;

};

class HTTParticle{

  public:

  HTTParticle(){ clear();}

  ~HTTParticle(){}

  void clear();

  ///Data member setters.
  void setP4(const TLorentzVector &aP4) { p4 = aP4;}

  void setChargedP4(const TLorentzVector &aP4) { chargedP4 = aP4;}

  void setNeutralP4(const TLorentzVector &aP4) { neutralP4 = aP4;}

  void setPCA(const TVector3 &aV3) {pca = aV3;}

  void setPCARefitPV(const TVector3 &aV3) {pcaRefitPV = aV3;}

  void setPCAGenPV(const TVector3 &aV3) {pcaGenPV = aV3;}

  void setSV(const TVector3 &aSV) {sv = aSV;}

  void setProperties(const std::vector<Double_t> & aProperties) { properties = aProperties;}

  ///Data member getters.
  const TLorentzVector & getP4(HTTAnalysis::sysEffects type=HTTAnalysis::NOMINAL) const {return getSystScaleP4(type);}

  const TLorentzVector & getChargedP4() const {return chargedP4;}

  const TLorentzVector getNeutralP4() const {return neutralP4;}

  const TVector3 & getPCA() const {return pca;}

  const TVector3 & getPCARefitPV() const {return pcaRefitPV;}

  const TVector3 & getPCAGenPV() const {return pcaGenPV;}

  const TVector3 & getSV() const {return sv;}

  int getPDGid() const {return getProperty(PropertyEnum::PDGId);}

  int getCharge() const {return getProperty(PropertyEnum::charge);}

  Double_t getProperty(PropertyEnum index) const {return (unsigned int)index<properties.size()?  properties[(unsigned int)index]: -999;}

  bool hasTriggerMatch(TriggerEnum index) const {return (unsigned int)getProperty(PropertyEnum::isGoodTriggerType)& (1<<(unsigned int)index) &&
                                                        (unsigned int)getProperty(PropertyEnum::FilterFired)& (1<<(unsigned int)index);}
 private:

  ///Return four-momentum modified according DATA/MC energy scale factors.
  const TLorentzVector & getNominalShiftedP4() const;

  ///Return four-momentum modified according to given systematic effect.
  ///The method recognises particle type, e.g. muons are not affected by
  ///TES variations etc.
  const TLorentzVector & getSystScaleP4(HTTAnalysis::sysEffects type=HTTAnalysis::NOMINAL) const;

  ///Return four-momentum shifted with scale.
  ///Shift modifies three-momentum transverse part only, leaving mass constant.
  const TLorentzVector & getShiftedP4(float scale, bool preserveMass=true) const;

  ///Nominal (as recontructed) four-momentum
  TLorentzVector p4;

  ///Scaled four-momentum cache;
  mutable TLorentzVector p4Cache;
  mutable HTTAnalysis::sysEffects lastSystEffect;

  ///Charged and neutral components four-momentum
  TLorentzVector chargedP4, neutralP4;

  ///Vectors from primary vertex to point of closest approach (PCA)
  ///calculated with respect to AOD vertex, refitted and generated vertex.
  TVector3 pca, pcaRefitPV, pcaGenPV;

  ///Decay vertex
  TVector3 sv;

  ///Vector of vaious particle properties.
  ///Index generated automatically during conversion from
  ///LLR ntuple format
  std::vector<Double_t> properties;

  static constexpr float TES_1p=-0.005;
  static constexpr float TES_1ppi0=0.011;
  static constexpr float TES_3p=0.006;
  static constexpr float TES = 0.012;
  static constexpr float EES = 0.03;
  static constexpr float MES = 0.03;
  
};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTPair{

 public:

  HTTPair(){ clear();}

  ~HTTPair(){}

  void clear();

  ///Data member setters.
  void setP4(const TLorentzVector &aP4, HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) {p4Vector[(unsigned int)type] = aP4;}

  void setLeg1P4(const TLorentzVector &aP4, HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) {leg1p4Vector[(unsigned int)type] = aP4;}

  void setLeg2P4(const TLorentzVector &aP4, HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) {leg2p4Vector[(unsigned int)type] = aP4;}
  
  void setMET(const TVector2 &aVector) {met = aVector;}

  void setSVMET(const TVector2 &aVector, HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) {svMetVector[(unsigned int)type] = aVector;}

  void setMTLeg1(const float & aMT) {mtLeg1 = aMT;}

  void setMTLeg2(const float & aMT) {mtLeg2 = aMT;}

  void setLeg1(const HTTParticle &aParticle){leg1 = aParticle;}

  void setLeg2(const HTTParticle &aParticle){leg2 = aParticle;}

  void setMETMatrix(float m00, float m01, float m10, float m11) {metMatrix.push_back(m00); metMatrix.push_back(m01); metMatrix.push_back(m10); metMatrix.push_back(m11);}

  ///Data member getters.
  const TLorentzVector & getP4(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const;

  const TLorentzVector & getLeg1P4(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const;

  const TLorentzVector & getLeg2P4(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const;

  const TVector2 & getMET(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const {return getSystScaleMET(type);}

  const TVector2 & getSVMET(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const {return svMetVector[(unsigned int)type];}

  float getMTLeg1(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const {return getSystScaleMT(leg1, type);}

  float getMTLeg2(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const {return getSystScaleMT(leg2, type);}

  const HTTParticle & getLeg1() const {return leg1;}

  const HTTParticle & getLeg2() const {return leg2;}

  const HTTParticle & getMuon() const {return abs(leg1.getPDGid())==13 ? leg1 : leg2; }

  const HTTParticle & getTau() const {return abs(leg1.getPDGid())==15 ? leg1 : leg2; }

  float getMTMuon(HTTAnalysis::sysEffects type = HTTAnalysis::NOMINAL) const {return abs(leg1.getPDGid())==13 ? getMTLeg1(type) : getMTLeg2(type); }

  std::vector<float> getMETMatrix() const {return metMatrix;}

 private:

  ///Return MET modified according to given systematic effect.
  ///The MET is corrected for accorging leptons corrections.
  ///The recoil correctino is not updated.
  const TVector2 & getSystScaleMET(HTTAnalysis::sysEffects type=HTTAnalysis::NOMINAL) const;

  ///Return transverse mass caluculated according to the scale shifts.
  float getSystScaleMT(const HTTParticle &aPerticle,
		       HTTAnalysis::sysEffects type=HTTAnalysis::NOMINAL) const;

  ///Nominal met as calculated from PF.
  ///Includes recoil corrections.
  TVector2 met;

  ///Scaled four-momentum cache;
  mutable TVector2 metCache;
  mutable HTTAnalysis::sysEffects lastSystEffect;
  mutable float mtCache;

  ///Vectors holding p4 and MET for
  ///for various scale variances.
  std::vector<TLorentzVector> p4Vector;
  std::vector<TLorentzVector> leg1p4Vector;
  std::vector<TLorentzVector> leg2p4Vector;
  std::vector<TVector2> svMetVector;

  //MVAMET covariance matrix in order 00,01,10,11
  std::vector<float> metMatrix;

  ///MT calculated for (leg1,MET) and (leg2,MET)
  float mtLeg1, mtLeg2;

  ///Pair legs
  HTTParticle leg1, leg2;

};

#endif
