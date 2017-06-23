#ifndef RootAnalysis_HTTSynchNTuple_H
#define RootAnalysis_HTTSynchNTuple_H

#include <string>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TRandom3.h"

#include "Analyzer.h"

class EventProxyHTT;
class HTTHistograms;

//For btag calibration
class BTagCalibration;
class BTagCalibrationReader;

class HTTSynchNTuple: public Analyzer{

 public:
  
  HTTSynchNTuple(const std::string & aName, const std::string & aDecayMode = "None");

  virtual ~HTTSynchNTuple();
  
  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *);
  
  void clearTTreeVariables();

  virtual void fillEventID(const HTTEvent &event);
  virtual void fillLegs(const HTTParticle &leg1, const HTTParticle &leg2);
  virtual void fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2);
  virtual void fillPair(const HTTEvent &event, HTTPair &pair);
  virtual void fillJets(const std::vector<HTTParticle> &jets);
  virtual void fillVetoes(const HTTEvent &event);

  virtual bool selectEvent(const HTTEvent &event, HTTPair &pair);

  virtual void initializeCorrections();
  virtual float getPUWeight(float nPU);

  virtual void initializeBTagCorrections();
  virtual bool promoteBJet(const HTTParticle &jet);
  virtual Analyzer* clone() const;

  bool filter() const { return filterEvent_;};

  const std::string & decayMode() const { return decayMode_; };

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  HTTHistograms *myHistos_;
  
  //should this HTTSynchNTuple be able to filter events
  bool filterEvent_;

  std::string decayMode_;

  //event counter
  //Int_t i_ = 0;

  ///Reconstructed objects selected for given event.
  //HTTEvent aEvent;
  //HTTPair aPair;

  HTTParticle leg1, leg2;

  ///Histograms with lepton corrections
  TH2F *h2DMuonIdCorrections;
  TH3F *h3DMuonIsoCorrections, *h3DMuonTrgCorrections, *h3DMuonXTrgCorrections;
  TH1F *h1DMuonTrkCorrections;
  TH3F *h3DTauCorrections;
  TH2F *h2DTauTrgGenuineCorrections, *h2DTauTrgFakeCorrections;
  TH2F *h2DTauXTrgGenuineCorrections, *h2DTauXTrgFakeCorrections;

  //For PU
  TFile *puDataFile_, *puMCFile_;
  std::vector<TH1F*> hPUVec_;

  //For btag calibration
  BTagCalibration *calib;
  BTagCalibrationReader *reader;
  TFile *btagEffFile_;
  TH2F *btag_eff_b_, *btag_eff_c_, *btag_eff_oth_;
  TRandom3 *rand_;

  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  //variables stored in TTree
  //Float_t muonPt;
  //event ID variables
  ULong64_t run;
  ULong64_t lumi;
  ULong64_t evt;
  //Pielup
  ULong64_t npv;
  Float_t npu;
  Float_t rho;
  Float_t puweight;
  //Weights
  Float_t trackingweight_1;
  Float_t trackingweight_2;
  Float_t effweight;
  Float_t weight;
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  Float_t pt_1;
  Float_t phi_1;
  Float_t eta_1;
  Float_t m_1;
  Float_t q_1;
  Float_t d0_1;
  Float_t dZ_1;
  Float_t mt_1;
  Float_t pfmt_1;
  Float_t puppimt_1;
  Float_t iso_1;
  Float_t id_e_mva_nt_loose_1;
  Int_t gen_match_1;
  Float_t againstElectronLooseMVA6_1;
  Float_t againstElectronMediumMVA6_1;
  Float_t againstElectronTightMVA6_1;
  Float_t againstElectronVLooseMVA6_1;
  Float_t againstElectronVTightMVA6_1;
  Float_t againstMuonLoose3_1;
  Float_t againstMuonTight3_1;
  Float_t byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t byIsolationMVA3newDMwoLTraw_1;
  Float_t byIsolationMVA3oldDMwoLTraw_1;
  Float_t byIsolationMVA3newDMwLTraw_1;
  Float_t byIsolationMVA3oldDMwLTraw_1;
  Float_t chargedIsoPtSum_1;
  Bool_t decayModeFindingOldDMs_1;
  Float_t neutralIsoPtSum_1;
  Float_t puCorrPtSum_1;
  Float_t trigweight_1;
  Float_t idisoweight_1;
  Int_t tau_decay_mode_1;
  Bool_t mva_olddm_medium_1;
  Bool_t mva_olddm_tight_1;
  Bool_t mva_olddm_vtight_1;
  //Leg 2 (trailing tau for tt, tau for et,mt, muon for em)
  Float_t pt_2;
  Float_t phi_2;
  Float_t eta_2;
  Float_t m_2;
  Float_t q_2;
  Float_t d0_2;
  Float_t dZ_2;
  Float_t mt_2;
  Float_t pfmt_2;
  Float_t puppimt_2;
  Float_t iso_2;
  Float_t id_e_mva_nt_loose_2;
  Int_t gen_match_2;
  Float_t againstElectronLooseMVA6_2;
  Float_t againstElectronMediumMVA6_2;
  Float_t againstElectronTightMVA6_2;
  Float_t againstElectronVLooseMVA6_2;
  Float_t againstElectronVTightMVA6_2;
  Float_t againstMuonLoose3_2;
  Float_t againstMuonTight3_2;
  Float_t byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t byIsolationMVA3newDMwoLTraw_2;
  Float_t byIsolationMVA3oldDMwoLTraw_2;
  Float_t byIsolationMVA3newDMwLTraw_2;
  Float_t byIsolationMVA3oldDMwLTraw_2;
  Float_t chargedIsoPtSum_2;
  Bool_t decayModeFindingOldDMs_2;
  Float_t neutralIsoPtSum_2;
  Float_t puCorrPtSum_2;
  Float_t trigweight_2;
  Float_t idisoweight_2;
  Int_t tau_decay_mode_2;
  Bool_t mva_olddm_medium_2;
  Bool_t mva_olddm_tight_2;
  Bool_t mva_olddm_vtight_2;
  //di-tau system
  Float_t pt_tt;
  Float_t mt_tot;
  Float_t m_vis;
  Float_t m_sv;
  Float_t mt_sv;
  Float_t pt_sv;
  Float_t eta_sv;
  Float_t phi_sv;
  Bool_t os;
  //MET
  Float_t met;
  Float_t metphi;
  Float_t pfmet;
  Float_t pfmetphi;
  Float_t puppimet;
  Float_t puppimetphi;
  Float_t mvamet;
  Float_t mvametphi;
  Float_t pzetavis;
  Float_t pzetamiss;
  Float_t pfpzetamiss;
  Float_t puppipzetamiss;
  Float_t mvacov00;
  Float_t mvacov01;
  Float_t mvacov10;
  Float_t mvacov11;
  Float_t metcov00;
  Float_t metcov01;
  Float_t metcov10;
  Float_t metcov11;
  //VBF system
  Float_t mjj;
  Float_t jdeta;
  Float_t njetingap;
  Float_t njetingap20;
  Float_t jdphi;
  //additional jets
  Float_t nbtag;
  Float_t njets;
  Float_t njetspt20;
  //leading jet sorted by pt
  Float_t jpt_1;
  Float_t jeta_1;
  Float_t jphi_1;
  Float_t jrawf_1;
  Float_t jmva_1;
  //trailing jet sorted by pt
  Float_t jpt_2;
  Float_t jeta_2;
  Float_t jphi_2;
  Float_t jrawf_2;
  Float_t jmva_2;
  //leading b-jet sorted by pt
  Float_t bpt_1;
  Float_t beta_1;
  Float_t bphi_1;
  Float_t brawf_1;
  Float_t bmva_1;
  Float_t bcsv_1;
  //trailing b-jet sorted by pt
  Float_t bpt_2;
  Float_t beta_2;
  Float_t bphi_2;
  Float_t brawf_2;
  Float_t bmva_2;
  Float_t bcsv_2;
  //Extra lepton vetos
  Bool_t dilepton_veto;
  Bool_t extraelec_veto;
  Bool_t extramuon_veto;
  //Trigger
  Bool_t trg_singlemuon;
  Bool_t trg_singleelectron;
  Bool_t trg_singletau_1;
  Bool_t trg_singletau_2;
  Bool_t trg_doubletau;
  Bool_t trg_muonelectron;
  
 private:

  std::map<std::string,int> *idMasks_;
	
};

#endif
