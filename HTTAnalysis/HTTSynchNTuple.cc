#include <sstream>
#include <algorithm>
#include <cstdlib>

#include "HTTSynchNTuple.h"
#include "HTTHistograms.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooBinning.h"

#include "BTagCalibrationStandalone.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::HTTSynchNTuple(const std::string & aName, const std::string & aDecayMode):Analyzer(aName){
  decayMode_ = aDecayMode;
  std::string cacheDirPath = "/tmp/"+std::string(std::getenv("USER"))+"/";
  TFile::SetCacheFileDir(cacheDirPath.c_str());
  //PU
  std::string dataPUFileName = "http://akalinow.web.cern.ch/akalinow/Data_Pileup_2016_271036-284044_13TeVMoriond17_23Sep2016ReReco_69p2mbMinBiasXS.root";
  puDataFile_ = TFile::Open(dataPUFileName.c_str(),"CACHEREAD");
  std::string mcPUFileName = "http://akalinow.web.cern.ch/akalinow/MC_Moriond17_PU25ns_V1.root";
  puMCFile_ = TFile::Open(mcPUFileName.c_str(),"CACHEREAD");
  //Corrections
  //initializeCorrections();//MB needed??
  h2DMuonIdCorrections = 0;
  h3DMuonIsoCorrections = 0;
  h3DMuonTrgCorrections = 0;
  h3DMuonXTrgCorrections = 0;
  h1DMuonTrkCorrections = 0;
  h2DTauTrgGenuineCorrections = 0;
  h2DTauTrgFakeCorrections = 0;
  h2DTauXTrgGenuineCorrections = 0;
  h2DTauXTrgFakeCorrections = 0;
  h3DTauCorrections = 0;
  initializeBTagCorrections();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::~HTTSynchNTuple(){
  if(myHistos_) delete myHistos_;
  if(idMasks_) delete idMasks_;
  //PU
  if(puDataFile_) delete puDataFile_;
  if(puMCFile_) delete puMCFile_;
  //Corrections
  if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
  if(h3DMuonIsoCorrections) delete h3DMuonIsoCorrections;
  if(h3DMuonTrgCorrections) delete h3DMuonTrgCorrections;
  if(h3DMuonXTrgCorrections) delete h3DMuonXTrgCorrections;
  if(h1DMuonTrkCorrections) delete h1DMuonTrkCorrections;
  if(h2DTauTrgGenuineCorrections) delete h2DTauTrgGenuineCorrections;
  if(h2DTauTrgFakeCorrections) delete h2DTauTrgFakeCorrections;
  if(h2DTauXTrgGenuineCorrections) delete h2DTauXTrgGenuineCorrections;
  if(h2DTauXTrgFakeCorrections) delete h2DTauXTrgFakeCorrections;
  if(h3DTauCorrections) delete h3DTauCorrections;
  //btagging
  if(calib) delete calib;
  if(reader) delete reader;
  if(btagEffFile_) delete btagEffFile_;
  if(rand_) delete rand_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTuple::clone() const{

  HTTSynchNTuple* clone = new HTTSynchNTuple(name(),decayMode());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::initialize(TDirectory* aDir,
			      pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "HTTSynchNTuple"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(aDir, selectionFlavours_, decayMode_);

  idMasks_ = new std::map<std::string,int>;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::finalize(){

  myHistos_->finalizeHistograms(std::vector<const HTTAnalysis::eventCategory*>());

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::clearTTreeVariables(){

  //event ID variables
  run = -999;
  lumi = -999;
  evt = -999;
  //weights
  effweight = -999;
  weight = -999;
  //Pielup
  npv = -999;
  npu = -999;
  rho = -999;
  puweight = -999;
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  pt_1 = -999;
  phi_1 = -999;
  eta_1 = -999;
  m_1 = -999;
  q_1 = -999;
  d0_1 = -999;
  dZ_1 = -999;
  mt_1 = -999;
  pfmt_1 = -999;
  puppimt_1 = -999;
  iso_1 = -999;
  id_e_mva_nt_loose_1 = -999;
  gen_match_1 = -999;
  againstElectronLooseMVA6_1 = -999;
  againstElectronMediumMVA6_1 = -999;
  againstElectronTightMVA6_1 = -999;
  againstElectronVLooseMVA6_1 = -999;
  againstElectronVTightMVA6_1 = -999;
  againstMuonLoose3_1 = -999;
  againstMuonTight3_1 = -999;
  byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -999;
  byIsolationMVA3newDMwoLTraw_1 = -999;
  byIsolationMVA3oldDMwoLTraw_1 = -999;
  byIsolationMVA3newDMwLTraw_1 = -999;
  byIsolationMVA3oldDMwLTraw_1 = -999;
  chargedIsoPtSum_1 = -999;
  decayModeFindingOldDMs_1 = 0;
  neutralIsoPtSum_1 = -999;
  puCorrPtSum_1 = -999;
  trigweight_1 = -999;
  idisoweight_1 = -999;
  tau_decay_mode_1 = -999 ;
  mva_olddm_medium_1 = 0;
  mva_olddm_tight_1 = 0;
  mva_olddm_vtight_1 = 0;
  //Leg 2 (trailing tau for tt; tau for et,mt; muon for em)
  pt_2 = -999;
  phi_2 = -999;
  eta_2 = -999;
  m_2 = -999;
  q_2 = -999;
  d0_2 = -999;
  dZ_2 = -999;
  mt_2 = -999;
  pfmt_2 = -999;
  puppimt_2 = -999;
  iso_2 = -999;
  id_e_mva_nt_loose_2 = -999;
  gen_match_2 = -999;
  againstElectronLooseMVA6_2 = -999;
  againstElectronMediumMVA6_2 = -999;
  againstElectronTightMVA6_2 = -999;
  againstElectronVLooseMVA6_2 = -999;
  againstElectronVTightMVA6_2 = -999;
  againstMuonLoose3_2 = -999;
  againstMuonTight3_2 = -999;
  byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = -999;
  byIsolationMVA3newDMwoLTraw_2 = -999;
  byIsolationMVA3oldDMwoLTraw_2 = -999;
  byIsolationMVA3newDMwLTraw_2 = -999;
  byIsolationMVA3oldDMwLTraw_2 = -999;
  chargedIsoPtSum_2 = -999;
  decayModeFindingOldDMs_2 = 0;
  neutralIsoPtSum_2 = -999;
  puCorrPtSum_2 = -999;
  trigweight_2 = -999;
  idisoweight_2 = -999;
  tau_decay_mode_2 = -999 ;
  mva_olddm_medium_2 = 0;
  mva_olddm_tight_2 = 0;
  mva_olddm_vtight_2 = 0;
  //di-tau system
  pt_tt = -999;
  mt_tot = -999;
  m_vis = -999;
  m_sv = -999;
  mt_sv = -999;
  mt_sv = -999;
  pt_sv = -999;
  eta_sv = -999;
  phi_sv = -999;
  os = 0;
  //MET
  met = -999;
  metphi = -999;
  pfmet = -999;
  pfmetphi = -999;
  puppimet = -999;
  puppimetphi = -999;
  mvamet = -999;
  mvametphi = -999;
  pzetavis = -999;
  pzetamiss = -999;
  pfpzetamiss = -999;
  puppipzetamiss = -999;
  mvacov00 = -999;
  mvacov01 = -999;
  mvacov10 = -999;
  mvacov11 = -999;
  metcov00 = -999;
  metcov01 = -999;
  metcov10 = -999;
  metcov11 = -999;
  //VBF system
  mjj = -999;
  jdeta = -999;
  njetingap = -999;
  njetingap20 = -999;
  jdphi = -999;
  //additional jets
  nbtag = -999;
  njets = -999;
  njetspt20 = -999;
  //leading jet sorted by pt
  jpt_1 = -999;
  jeta_1 = -999;
  jphi_1 = -999;
  jrawf_1 = -999;
  jmva_1 = -999;
  //trailing jet sorted by pt
  jpt_2 = -999;
  jeta_2 = -999;
  jphi_2 = -999;
  jrawf_2 = -999;
  jmva_2 = -999;
  //leading b-jet sorted by pt
  bpt_1 = -999;
  beta_1 = -999;
  bphi_1 = -999;
  brawf_1 = -999;
  bmva_1 = -999;
  bcsv_1 = -999;
  //trailing b-jet sorted by pt
  bpt_2 = -999;
  beta_2 = -999;
  bphi_2 = -999;
  brawf_2 = -999;
  bmva_2 = -999;
  bcsv_2 = -999;
  //Extra lepton vetos
  dilepton_veto = 0;
  extraelec_veto = 0;
  extramuon_veto = 0;
  //Trigger
  trg_singlemuon = 0;
  trg_singleelectron = 0;
  trg_singletau_1 = 0;
  trg_singletau_2 = 0;
  trg_doubletau = 0;
  trg_muonelectron = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::addBranch(TTree *tree){

  //event ID variables
  tree->Branch("run",&run);
  tree->Branch("lumi",&lumi);
  tree->Branch("evt",&evt);
  //weights
  tree->Branch("effweight",&effweight);
  tree->Branch("weight",&weight);
  //Pielup
  tree->Branch("npv",&npv);
  tree->Branch("npu",&npu);
  tree->Branch("rho",&rho);
  tree->Branch("puweight",&puweight);
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  tree->Branch("pt_1",&pt_1);
  tree->Branch("phi_1",&phi_1);
  tree->Branch("eta_1",&eta_1);
  tree->Branch("m_1",&m_1);
  tree->Branch("q_1",&q_1);
  tree->Branch("d0_1",&d0_1);
  tree->Branch("dZ_1",&dZ_1);
  tree->Branch("mt_1",&mt_1);
  tree->Branch("pfmt_1",&pfmt_1);
  tree->Branch("puppimt_1",&puppimt_1);
  tree->Branch("iso_1",&iso_1);
  tree->Branch("id_e_mva_nt_loose_1",&id_e_mva_nt_loose_1);
  tree->Branch("gen_match_1",&gen_match_1);
  tree->Branch("againstElectronLooseMVA6_1",&againstElectronLooseMVA6_1);
  tree->Branch("againstElectronMediumMVA6_1",&againstElectronMediumMVA6_1);
  tree->Branch("againstElectronTightMVA6_1",&againstElectronTightMVA6_1);
  tree->Branch("againstElectronVLooseMVA6_1",&againstElectronVLooseMVA6_1);
  tree->Branch("againstElectronVTightMVA6_1",&againstElectronVTightMVA6_1);
  tree->Branch("againstMuonLoose3_1",&againstMuonLoose3_1);
  tree->Branch("againstMuonTight3_1",&againstMuonTight3_1);
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
  tree->Branch("byIsolationMVA3newDMwoLTraw_1",&byIsolationMVA3newDMwoLTraw_1);
  tree->Branch("byIsolationMVA3oldDMwoLTraw_1",&byIsolationMVA3oldDMwoLTraw_1);
  tree->Branch("byIsolationMVA3newDMwLTraw_1",&byIsolationMVA3newDMwLTraw_1);
  tree->Branch("byIsolationMVA3oldDMwLTraw_1",&byIsolationMVA3oldDMwLTraw_1);
  tree->Branch("chargedIsoPtSum_1",&chargedIsoPtSum_1);
  tree->Branch("decayModeFindingOldDMs_1",&decayModeFindingOldDMs_1);
  tree->Branch("neutralIsoPtSum_1",&neutralIsoPtSum_1);
  tree->Branch("puCorrPtSum_1",&puCorrPtSum_1);
  tree->Branch("trigweight_1",&trigweight_1);
  tree->Branch("idisoweight_1",&idisoweight_1);
  tree->Branch("trackingweight_1",&trackingweight_1);
  tree->Branch("tau_decay_mode_1",&tau_decay_mode_1);
  tree->Branch("mva_olddm_medium_1",&mva_olddm_medium_1);
  tree->Branch("mva_olddm_tight_1",&mva_olddm_tight_1);
  tree->Branch("mva_olddm_vtight_1",&mva_olddm_vtight_1);
  //Leg 2 (trailing tau for tt, tau for et,mt, muon for em)
  tree->Branch("pt_2",&pt_2);
  tree->Branch("phi_2",&phi_2);
  tree->Branch("eta_2",&eta_2);
  tree->Branch("m_2",&m_2);
  tree->Branch("q_2",&q_2);
  tree->Branch("d0_2",&d0_2);
  tree->Branch("dZ_2",&dZ_2);
  tree->Branch("mt_2",&mt_2);
  tree->Branch("pfmt_2",&pfmt_2);
  tree->Branch("puppimt_2",&puppimt_2);
  tree->Branch("iso_2",&iso_2);
  tree->Branch("id_e_mva_nt_loose_2",&id_e_mva_nt_loose_2);
  tree->Branch("gen_match_2",&gen_match_2);
  tree->Branch("againstElectronLooseMVA6_2",&againstElectronLooseMVA6_2);
  tree->Branch("againstElectronMediumMVA6_2",&againstElectronMediumMVA6_2);
  tree->Branch("againstElectronTightMVA6_2",&againstElectronTightMVA6_2);
  tree->Branch("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2);
  tree->Branch("againstElectronVTightMVA6_2",&againstElectronVTightMVA6_2);
  tree->Branch("againstMuonLoose3_2",&againstMuonLoose3_2);
  tree->Branch("againstMuonTight3_2",&againstMuonTight3_2);
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
  tree->Branch("byIsolationMVA3newDMwoLTraw_2",&byIsolationMVA3newDMwoLTraw_2);
  tree->Branch("byIsolationMVA3oldDMwoLTraw_2",&byIsolationMVA3oldDMwoLTraw_2);
  tree->Branch("byIsolationMVA3newDMwLTraw_2",&byIsolationMVA3newDMwLTraw_2);
  tree->Branch("byIsolationMVA3oldDMwLTraw_2",&byIsolationMVA3oldDMwLTraw_2);
  tree->Branch("chargedIsoPtSum_2",&chargedIsoPtSum_2);
  tree->Branch("decayModeFindingOldDMs_2",&decayModeFindingOldDMs_2);
  tree->Branch("neutralIsoPtSum_2",&neutralIsoPtSum_2);
  tree->Branch("puCorrPtSum_2",&puCorrPtSum_2);
  tree->Branch("trigweight_2",&trigweight_2);
  tree->Branch("idisoweight_2",&idisoweight_2);
  tree->Branch("trackingweight_2",&trackingweight_2);
  tree->Branch("tau_decay_mode_2",&tau_decay_mode_2);
  tree->Branch("mva_olddm_medium_2",&mva_olddm_medium_2);
  tree->Branch("mva_olddm_tight_2",&mva_olddm_tight_2);
  tree->Branch("mva_olddm_vtight_2",&mva_olddm_vtight_2);
  //di-tau system
  tree->Branch("pt_tt",&pt_tt);
  tree->Branch("mt_tot",&mt_tot);
  tree->Branch("m_vis",&m_vis);
  tree->Branch("m_sv",&m_sv);
  tree->Branch("mt_sv",&mt_sv);
  tree->Branch("pt_sv",&pt_sv);
  tree->Branch("eta_sv",&eta_sv);
  tree->Branch("phi_sv",&phi_sv);
  tree->Branch("os",&os);
  //MET
  tree->Branch("met",&met);
  tree->Branch("metphi",&metphi);
  tree->Branch("pfmet",&pfmet);
  tree->Branch("pfmetphi",&pfmetphi);
  tree->Branch("puppimet",&puppimet);
  tree->Branch("puppimetphi",&puppimetphi);
  tree->Branch("mvamet",&mvamet);
  tree->Branch("mvametphi",&mvametphi);
  tree->Branch("pzetavis",&pzetavis);
  tree->Branch("pzetamiss",&pzetamiss);
  tree->Branch("pfpzetamiss",&pfpzetamiss);
  tree->Branch("puppipzetamiss",&puppipzetamiss);
  tree->Branch("mvacov00",&mvacov00);
  tree->Branch("mvacov01",&mvacov01);
  tree->Branch("mvacov10",&mvacov10);
  tree->Branch("mvacov11",&mvacov11);
  tree->Branch("metcov00",&metcov00);
  tree->Branch("metcov01",&metcov01);
  tree->Branch("metcov10",&metcov10);
  tree->Branch("metcov11",&metcov11);
  //VBF system
  tree->Branch("mjj",&mjj);
  tree->Branch("jdeta",&jdeta);
  tree->Branch("njetingap",&njetingap);
  tree->Branch("njetingap20",&njetingap20);
  tree->Branch("jdphi",&jdphi);
  //additional jets
  tree->Branch("nbtag",&nbtag);
  tree->Branch("njets",&njets);
  tree->Branch("njetspt20",&njetspt20);
  //leading jet sorted by pt
  tree->Branch("jpt_1",&jpt_1);
  tree->Branch("jeta_1",&jeta_1);
  tree->Branch("jphi_1",&jphi_1);
  tree->Branch("jrawf_1",&jrawf_1);
  tree->Branch("jmva_1",&jmva_1);
  //trailing jet sorted by pt
  tree->Branch("jpt_2",&jpt_2);
  tree->Branch("jeta_2",&jeta_2);
  tree->Branch("jphi_2",&jphi_2);
  tree->Branch("jrawf_2",&jrawf_2);
  tree->Branch("jmva_2",&jmva_2);
  //leading b-jet sorted by pt
  tree->Branch("bpt_1",&bpt_1);
  tree->Branch("beta_1",&beta_1);
  tree->Branch("bphi_1",&bphi_1);
  tree->Branch("brawf_1",&brawf_1);
  tree->Branch("bmva_1",&bmva_1);
  tree->Branch("bcsv_1",&bcsv_1);
  //trailing b-jet sorted by pt
  tree->Branch("bpt_2",&bpt_2);
  tree->Branch("beta_2",&beta_2);
  tree->Branch("bphi_2",&bphi_2);
  tree->Branch("brawf_2",&brawf_2);
  tree->Branch("bmva_2",&bmva_2);
  tree->Branch("bcsv_2",&bcsv_2);
  //Extra lepton vetos
  tree->Branch("dilepton_veto",&dilepton_veto);
  tree->Branch("extraelec_veto",&extraelec_veto);
  tree->Branch("extramuon_veto",&extramuon_veto);
  //Trigger
  tree->Branch("trg_singlemuon",&trg_singlemuon);
  tree->Branch("trg_singleelectron",&trg_singleelectron);
  tree->Branch("trg_singletau_1",&trg_singletau_1);
  tree->Branch("trg_singletau_2",&trg_singletau_2);
  tree->Branch("trg_doubletau",&trg_doubletau);
  tree->Branch("trg_muonelectron",&trg_muonelectron);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTSynchNTuple::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  HTTEvent &aEvent = *myEventProxy.event;
  HTTPair &aPair = (*myEventProxy.pairs)[0];
  std::vector<HTTParticle> &aJets = *myEventProxy.jets;

  //i_++;
  //if(i_%1000==0){
    //std::cout<<i_<<std::endl;
  //}

  //Set map of tau-Id masks (only ones)
  //for decay channels with taus
  if(idMasks_->size()==0 && decayMode_.find("Tau")!=std::string::npos){
    idMasks_->insert(std::pair<std::string,int>("againstMuonLoose3",0) );
    idMasks_->insert(std::pair<std::string,int>("againstMuonTight3",0) );
    idMasks_->insert(std::pair<std::string,int>("againstElectronVLooseMVA6",0) );
    idMasks_->insert(std::pair<std::string,int>("againstElectronLooseMVA6",0) );
    idMasks_->insert(std::pair<std::string,int>("againstElectronMediumMVA6",0) );
    idMasks_->insert(std::pair<std::string,int>("againstElectronTightMVA6",0) );
    idMasks_->insert(std::pair<std::string,int>("againstElectronVTightMVA6",0) );
    idMasks_->insert(std::pair<std::string,int>("byMediumIsolationMVArun2v1DBoldDMwLT",0) );
    idMasks_->insert(std::pair<std::string,int>("byTightIsolationMVArun2v1DBoldDMwLT",0) );
    idMasks_->insert(std::pair<std::string,int>("byVTightIsolationMVArun2v1DBoldDMwLT",0) );
    for(unsigned int iBit=0; iBit<aEvent.ntauIds; iBit++) {
      for(std::map<std::string,int>::iterator it=idMasks_->begin(); it!=idMasks_->end(); ++it){
	if(aEvent.tauIDStrings[iBit]==it->first) {
	  it->second |= (1<<iBit);
	  //debug
	  ///std::cout<<"\t"<<it->first<<" : "<<it->second<<" : "<<iBit<<std::endl;
	  break;
	}
      }
    }
  }
  //Apply an additional selection
  if( !selectEvent(aEvent,aPair) )
      return false;

  //Set default values to tree variables
  clearTTreeVariables();

  ///Filling TTree
  // Event ID variables
  fillEventID(aEvent);
  // Fill legs
  fillLegs(aPair.getLeg1(), aPair.getLeg2());
  // Fill di-tau system, but also MET which can be pairwise
  fillPair(aEvent, aPair);
  // Fill jet variables includig VBF system
  fillJets(aJets);
  // Fill extra lepton vetoes info
  fillVetoes(aEvent);
  // Weights
  Float_t w_leg1 = trackingweight_1*trigweight_1*idisoweight_1;
  Float_t w_leg2 = trackingweight_2*trigweight_2*idisoweight_2;
  effweight = w_leg1*w_leg2;
  weight = puweight*effweight;
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillEventID(const HTTEvent &event){

  run = event.getRunId();
  lumi = 0; //FIXME, not in ntuples
  evt = event.getEventId();
  //Pielup
  npv = event.getNPV();
  npu = event.getNPU();
  rho = 0; //FIXME, not in ntuples
  puweight = getPUWeight(npu);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillLegs(const HTTParticle &leg1, const HTTParticle &leg2){

  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  pt_1 = leg1.getP4().Pt();
  phi_1 = leg1.getP4().Phi();
  eta_1 = leg1.getP4().Eta();
  m_1 = leg1.getP4().M();
  if(std::abs(leg1.getPDGid())==13)
    m_1 = 0.1057; //muon mass
  else if(std::abs(leg1.getPDGid())==15 && leg1.getProperty(PropertyEnum::decayMode)==0)
    m_1 = 0.13957; //pi+/- mass
  if(m_1 < 0) m_1 = 0; //protect against numerical instability
  q_1 = leg1.getCharge();
  d0_1 = leg1.getProperty(PropertyEnum::dxy);
  dZ_1 = leg1.getProperty(PropertyEnum::dz);
  gen_match_1 = leg1.getProperty(PropertyEnum::mc_match);

  //Leg 2 (trailing tau for tt, electon for et,em muon for mt)
  pt_2 = leg2.getP4().Pt();
  phi_2 = leg2.getP4().Phi();
  eta_2 = leg2.getP4().Eta();
  m_2 = leg2.getP4().M();
  if(std::abs(leg2.getPDGid())==13)
    m_2 = 0.1057; //muon mass
  else if(std::abs(leg2.getPDGid())==15 && leg2.getProperty(PropertyEnum::decayMode)==0)
    m_2 = 0.13957; //pi+/- mass
  if(m_2 < 0) m_2 = 0; //protect against numerical instability
  q_2 = leg2.getCharge();
  d0_2 = leg2.getProperty(PropertyEnum::dxy);
  dZ_2 = leg2.getProperty(PropertyEnum::dz);
  gen_match_2 = leg2.getProperty(PropertyEnum::mc_match); //according to: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching

  os = q_1*q_2<0;

  //Decay channel specific
  fillLegsSpecific(leg1,leg2);
 }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2){

  if(!h2DMuonIdCorrections) initializeCorrections();//MB why it is needed??
  if(decayMode_=="MuTau"){
    //Specific implementation for the mu+tau decay channel
    //Leg1: muon
    iso_1 =  leg1.getProperty(PropertyEnum::combreliso);
    //weights
    int iBin;
    iBin = h3DMuonTrgCorrections->FindBin(std::min(pt_1,(Float_t)999.9), eta_1, std::min(iso_1,(Float_t)0.499));
    trigweight_1 = h3DMuonTrgCorrections->GetBinContent(iBin);
    if(pt_1<23){//xtrigger mu-tau
      iBin = h3DMuonXTrgCorrections->FindBin(std::min(pt_1,(Float_t)999.9), eta_1, std::min(iso_1,(Float_t)0.299));
      trigweight_1 = h3DMuonXTrgCorrections->GetBinContent(iBin);      
    }
    iBin = h2DMuonIdCorrections->FindBin(std::min(pt_1,(Float_t)999.9), eta_1);
    idisoweight_1 = h2DMuonIdCorrections->GetBinContent(iBin);
    iBin = h3DMuonIsoCorrections->FindBin(std::min(pt_1,(Float_t)999.9), eta_1, std::min(iso_1,(Float_t)0.499));
    idisoweight_1 *= h3DMuonIsoCorrections->GetBinContent(iBin);
    iBin = h1DMuonTrkCorrections->FindBin(eta_1);
    trackingweight_1 = h1DMuonTrkCorrections->GetBinContent(iBin);
    //trigger
    trg_singlemuon = ( leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
		       leg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
		       leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
		       leg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1) );

    //Leg2: tau
    iso_2 = leg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
    byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
    int TauID = (int)leg2.getProperty(PropertyEnum::tauID);
    againstElectronLooseMVA6_2 = idMasks_->find("againstElectronLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronLooseMVA6")->second) == idMasks_->find("againstElectronLooseMVA6")->second :
      0;
    againstElectronMediumMVA6_2 = idMasks_->find("againstElectronMediumMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronMediumMVA6")->second) == idMasks_->find("againstElectronMediumMVA6")->second :
      0;
    againstElectronTightMVA6_2 = idMasks_->find("againstElectronTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronTightMVA6")->second) == idMasks_->find("againstElectronTightMVA6")->second :
      0;
    againstElectronVLooseMVA6_2 = idMasks_->find("againstElectronVLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVLooseMVA6")->second) == idMasks_->find("againstElectronVLooseMVA6")->second :
      0;
    againstElectronVTightMVA6_2 = idMasks_->find("againstElectronVTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVTightMVA6")->second) == idMasks_->find("againstElectronVTightMVA6")->second :
      0;
    againstMuonLoose3_2 = idMasks_->find("againstMuonLoose3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonLoose3")->second) == idMasks_->find("againstMuonLoose3")->second :
      0;
    againstMuonTight3_2 = idMasks_->find("againstMuonTight3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonTight3")->second) == idMasks_->find("againstMuonTight3")->second :
      0;
    mva_olddm_medium_2 = idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_tight_2 = idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_vtight_2 = idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    /*
      chargedIsoPtSum_2;
      neutralIsoPtSum_2;
      puCorrPtSum_2;
    */
    tau_decay_mode_2 = leg2.getProperty(PropertyEnum::decayMode);
    decayModeFindingOldDMs_2 = (tau_decay_mode_2==0 || tau_decay_mode_2==1 || tau_decay_mode_2==2 || tau_decay_mode_2==10); //FIXME: is it possible to take ID directly?
    trigweight_2 = 1; //1, single mu
    if(pt_1<23){//xtrigger mu-tau
      if(gen_match_2==5){ //genuine tau
	iBin = h2DTauXTrgGenuineCorrections->FindBin(std::min(pt_2,(Float_t)999.9),eta_2);
	trigweight_2 = h2DTauXTrgGenuineCorrections->GetBinContent(iBin);
      }
      else{ //fake tau
	iBin = h2DTauXTrgFakeCorrections->FindBin(std::min(pt_2,(Float_t)999.9), eta_2);
	trigweight_2 = h2DTauXTrgFakeCorrections->GetBinContent(iBin);
      }  
    }
    idisoweight_2 = (gen_match_2==5 ? 0.95 : 1); //0.95 for genuine tau, otherwise 1
    trackingweight_2 = 1; //1 for tau
    trg_singletau_1 = leg2.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau120_Trk50_eta2p1);
    trg_singletau_2 = leg2.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau140_Trk50_eta2p1);

    return;
  }
  else if(decayMode_=="TauTau"){
    //Specific implementation for the tau+tau decay channel
    //Leg1: leading tau
    iso_1 = leg1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
    byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = leg1.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
    int TauID = (int)leg1.getProperty(PropertyEnum::tauID);
    againstElectronLooseMVA6_1 = idMasks_->find("againstElectronLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronLooseMVA6")->second) == idMasks_->find("againstElectronLooseMVA6")->second :
      0;
    againstElectronMediumMVA6_1 = idMasks_->find("againstElectronMediumMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronMediumMVA6")->second) == idMasks_->find("againstElectronMediumMVA6")->second :
      0;
    againstElectronTightMVA6_1 = idMasks_->find("againstElectronTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronTightMVA6")->second) == idMasks_->find("againstElectronTightMVA6")->second :
      0;
    againstElectronVLooseMVA6_1 = idMasks_->find("againstElectronVLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVLooseMVA6")->second) == idMasks_->find("againstElectronVLooseMVA6")->second :
      0;
    againstElectronVTightMVA6_1 = idMasks_->find("againstElectronVTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVTightMVA6")->second) == idMasks_->find("againstElectronVTightMVA6")->second :
      0;
    againstMuonLoose3_1 = idMasks_->find("againstMuonLoose3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonLoose3")->second) == idMasks_->find("againstMuonLoose3")->second :
      0;
    againstMuonTight3_1 = idMasks_->find("againstMuonTight3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonTight3")->second) == idMasks_->find("againstMuonTight3")->second :
      0;
    mva_olddm_medium_1 = idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_tight_1 = idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_vtight_1 = idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    /*
      chargedIsoPtSum_1;
      neutralIsoPtSum_1;
      puCorrPtSum_1;
    */
    tau_decay_mode_1 = leg1.getProperty(PropertyEnum::decayMode);
    decayModeFindingOldDMs_1 = (tau_decay_mode_1==0 || tau_decay_mode_1==1 || tau_decay_mode_1==2 || tau_decay_mode_1==10); //FIXME: is it possible to take ID directly?
    //Weights
    int iBin;
    int tau_decay_mode_tmp;
    tau_decay_mode_tmp  = tau_decay_mode_1==2 ? 1 : tau_decay_mode_1; //to handle very rare dm==2
    if(gen_match_1==5){ //genuine tau
      iBin = h2DTauTrgGenuineCorrections->FindBin(std::min(pt_1,(Float_t)999.9), tau_decay_mode_tmp);
      trigweight_1 = h2DTauTrgGenuineCorrections->GetBinContent(iBin);
      idisoweight_1 = 0.95;
    }
    else{ //fake tau
      iBin = h2DTauTrgFakeCorrections->FindBin(std::min(pt_1,(Float_t)999.9), tau_decay_mode_tmp);
      trigweight_1 = h2DTauTrgFakeCorrections->GetBinContent(iBin);
      idisoweight_1 = 1;
    }
    trackingweight_1 = 1;//1 for tau

    //Leg2: trailing tau
    iso_2 = leg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw);
    byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);
    TauID = (int)leg2.getProperty(PropertyEnum::tauID);
    againstElectronLooseMVA6_2 = idMasks_->find("againstElectronLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronLooseMVA6")->second) == idMasks_->find("againstElectronLooseMVA6")->second :
      0;
    againstElectronMediumMVA6_2 = idMasks_->find("againstElectronMediumMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronMediumMVA6")->second) == idMasks_->find("againstElectronMediumMVA6")->second :
      0;
    againstElectronTightMVA6_2 = idMasks_->find("againstElectronTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronTightMVA6")->second) == idMasks_->find("againstElectronTightMVA6")->second :
      0;
    againstElectronVLooseMVA6_2 = idMasks_->find("againstElectronVLooseMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVLooseMVA6")->second) == idMasks_->find("againstElectronVLooseMVA6")->second :
      0;
    againstElectronVTightMVA6_2 = idMasks_->find("againstElectronVTightMVA6")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstElectronVTightMVA6")->second) == idMasks_->find("againstElectronVTightMVA6")->second :
      0;
    againstMuonLoose3_2 = idMasks_->find("againstMuonLoose3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonLoose3")->second) == idMasks_->find("againstMuonLoose3")->second :
      0;
    againstMuonTight3_2 = idMasks_->find("againstMuonTight3")!=idMasks_->end() ?
      (TauID & idMasks_->find("againstMuonTight3")->second) == idMasks_->find("againstMuonTight3")->second :
      0;
    mva_olddm_medium_2 = idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byMediumIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_tight_2 = idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    mva_olddm_vtight_2 = idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")!=idMasks_->end() ?
      (TauID & idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second) == idMasks_->find("byVTightIsolationMVArun2v1DBoldDMwLT")->second :
      0;
    /*
      chargedIsoPtSum_2;
      neutralIsoPtSum_2;
      puCorrPtSum_2;
    */
    tau_decay_mode_2 = leg2.getProperty(PropertyEnum::decayMode);
    decayModeFindingOldDMs_2 = (tau_decay_mode_2==0 || tau_decay_mode_2==1 || tau_decay_mode_2==2 || tau_decay_mode_2==10); //FIXME: is it possible to take ID directly?
    //Weights
    tau_decay_mode_tmp  = tau_decay_mode_2==2 ? 1 : tau_decay_mode_2; //to handle very rare dm==2
    if(gen_match_2==5){ //genuine tau
      iBin = h2DTauTrgGenuineCorrections->FindBin(std::min(pt_2,(Float_t)999.9), tau_decay_mode_tmp);
      trigweight_2 = h2DTauTrgGenuineCorrections->GetBinContent(iBin);
      idisoweight_2 = 0.95;
    }
    else{ //fake tau
      iBin = h2DTauTrgFakeCorrections->FindBin(std::min(pt_2,(Float_t)999.9), tau_decay_mode_tmp);
      trigweight_2 = h2DTauTrgFakeCorrections->GetBinContent(iBin);
      idisoweight_2 = 1;
    }
    trackingweight_2 = 1;//1 for tau
    //triggers
    trg_singletau_1 = (leg1.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau120_Trk50_eta2p1) ||
		       leg2.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau120_Trk50_eta2p1) );
    trg_singletau_2 = (leg1.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau140_Trk50_eta2p1) ||
		       leg2.hasTriggerMatch(TriggerEnum::HLT_VLooseIsoPFTau140_Trk50_eta2p1) );
    trg_doubletau = ( (leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) &&
		       leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) ) ||
		      (leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg) &&
		       leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg) ) );


    return;
  }
  else if(decayMode_=="MuMu"){
    //Specific implementation for the mu+mu decay channel
    //Leg1: leading muon
    iso_1 =  leg1.getProperty(PropertyEnum::combreliso);
    trigweight_1 = 1; //FIXME, to be computed
    idisoweight_1 = 1; //FIXME
    trackingweight_1 = 1; //FIXME
    //Leg2: trailing muon
    iso_2 =  leg2.getProperty(PropertyEnum::combreliso);
    trigweight_2 = 1; //FIXME, to be computed
    idisoweight_2 = 1; //FIXME
    trackingweight_2 = 1; //FIXME

    return;
  }
  //Do nothing for undefined decayMode
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillPair(const HTTEvent &event, HTTPair &pair){

  //Legs
  TLorentzVector leg1P4 = pair.getLeg1().getP4();
  TLorentzVector leg2P4 = pair.getLeg2().getP4();

  //MET
  TLorentzVector pfmetP4(event.getMET().X(), event.getMET().Y(), 0, event.getMET().Mod());
  pfmet = event.getMET().Mod();
  pfmetphi = event.getMET().Phi_mpi_pi(event.getMET().Phi());
  pfmt_1 = TMath::Sqrt(2.*leg1P4.Pt()*pfmetP4.Pt()*(1.-TMath::Cos(leg1P4.Phi()-pfmetP4.Phi())));
  pfmt_2 = TMath::Sqrt(2.*leg2P4.Pt()*pfmetP4.Pt()*(1.-TMath::Cos(leg2P4.Phi()-pfmetP4.Phi())));
  //Pair met: by default called MVAMET which can be wrong...
  TLorentzVector mvametP4(pair.getMET().X(), pair.getMET().Y(), 0, pair.getMET().Mod());
  mvamet = pair.getMET().Mod();
  mvametphi = pair.getMET().Phi_mpi_pi(pair.getMET().Phi());
  met = mvamet;
  metphi = mvametphi;
  mt_1 = TMath::Sqrt(2.*leg1P4.Pt()*mvametP4.Pt()*(1.-TMath::Cos(leg1P4.Phi()-mvametP4.Phi())));
  mt_2 = TMath::Sqrt(2.*leg2P4.Pt()*mvametP4.Pt()*(1.-TMath::Cos(leg2P4.Phi()-mvametP4.Phi())));
  /*
    puppimet;
    puppimetphi;
    puppimt_1;
    puppimt_2;
  */

  //di-tau system
  pt_tt = (mvametP4 + leg1P4 + leg2P4).Pt();
  mt_tot = TMath::Sqrt( 2. * leg1P4.Pt() * mvametP4.Pt() * (1. - TMath::Cos( leg1P4.Phi() - mvametP4.Phi() ) ) + 2. * leg2P4.Pt() * mvametP4.Pt() * (1. - TMath::Cos( leg2P4.Phi() - mvametP4.Phi() ) ) + 2. * leg1P4.Pt() * leg2P4.Pt() * (1. - TMath::Cos(leg1P4.Phi() - leg2P4.Phi() ) ) );
  m_vis = (leg1P4 + leg2P4).M();
  m_sv = pair.getP4().M();

  float leg1Px = leg1P4.Px(), leg1Py = leg1P4.Py(), leg1Phi = leg1P4.Phi();
  float leg2Px = leg2P4.Px(), leg2Py = leg2P4.Py(), leg2Phi = leg2P4.Phi();
  pzetavis = ( (leg1Px+leg2Px)*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) + (leg1Py + leg2Py)*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  pzetamiss = ( pair.getMET().X()*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) +  pair.getMET().Y()*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  pfpzetamiss = ( event.getMET().X()*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) +  event.getMET().Y()*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  //puppipzetamiss;

  //The following does not like that pair is a const reference
  mvacov00 = pair.getMETMatrix().at(0);
  mvacov01 = pair.getMETMatrix().at(1);
  mvacov10 = pair.getMETMatrix().at(2);
  mvacov11 = pair.getMETMatrix().at(3);
  metcov00 = mvacov00;
  metcov01 = mvacov01;
  metcov10 = mvacov10;
  metcov11 = mvacov11;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillJets(const std::vector<HTTParticle> &jets){

  njetspt20 = jets.size();
  njets = 0;
  nbtag = 0;
  std::vector<HTTParticle> bjets;
  for(unsigned int iJet=0; iJet<jets.size(); ++iJet){
    if(jets.at(iJet).getP4().Pt()>30)
      njets++;
    if(std::abs(jets.at(iJet).getP4().Eta())<2.4 &&
       jets.at(iJet).getP4().Pt()>20 &&
       jets.at(iJet).getProperty(PropertyEnum::bCSVscore)>0.8484 && //FIXME 0.8484 Correct??
       promoteBJet(jets.at(iJet)) &&
       true){
      nbtag++;
      bjets.push_back(jets.at(iJet));
    }
  }
  if(jets.size()==0) return;

  HTTParticle aLeadingJet = jets.at(0);
  //leading jet sorted by pt
  jpt_1 = aLeadingJet.getP4().Pt();
  jeta_1 = aLeadingJet.getP4().Eta();
  jphi_1 = aLeadingJet.getP4().Phi();
  jrawf_1 =  aLeadingJet.getProperty(PropertyEnum::rawPt)/aLeadingJet.getP4().Pt();
  jmva_1 = aLeadingJet.getProperty(PropertyEnum::PUJetID);

  if(jets.size()>1){
    HTTParticle aTrailingJet = jets.at(1);

    //trailing jet sorted by pt
    jpt_2 = aTrailingJet.getP4().Pt();
    jeta_2 = aTrailingJet.getP4().Eta();
    jphi_2 = aTrailingJet.getP4().Phi();
    jrawf_2 = aTrailingJet.getProperty(PropertyEnum::rawPt)/aTrailingJet.getP4().Pt();
    jmva_2 = aTrailingJet.getProperty(PropertyEnum::PUJetID);

    //VBF system
    mjj = (aLeadingJet.getP4()+aTrailingJet.getP4()).M();
    jdeta = std::abs(aLeadingJet.getP4().Eta()-aTrailingJet.getP4().Eta());
    jdphi = aLeadingJet.getP4().Phi()-aTrailingJet.getP4().Phi();
    while(jdphi>TMath::Pi()) jdphi -= 2.*TMath::Pi();
    while(jdphi<=-TMath::Pi()) jdphi += 2.*TMath::Pi();
    jdphi = std::abs(jdphi);
    //jets in eta gap
    for(unsigned int iJet=2; iJet<jets.size(); ++iJet){
      if( (jets.at(iJet).getP4().Eta()>aLeadingJet.getP4().Eta()&&jets.at(iJet).getP4().Eta()<aTrailingJet.getP4().Eta()) ||
	  (jets.at(iJet).getP4().Eta()<aLeadingJet.getP4().Eta()&&jets.at(iJet).getP4().Eta()>aTrailingJet.getP4().Eta()) ){
	njetingap20++;
	if(jets.at(iJet).getP4().Pt()>30)
	  njetingap++;
      }
    }
  }
  //b-jets
  if(bjets.size()>0){
    //leading b-jet sorted by pt
    bpt_1 = bjets.at(0).getP4().Pt();
    beta_1 = bjets.at(0).getP4().Eta();
    bphi_1 = bjets.at(0).getP4().Phi();
    brawf_1 = bjets.at(0).getProperty(PropertyEnum::rawPt)/bjets.at(0).getP4().Pt();
    bmva_1 = bjets.at(0).getProperty(PropertyEnum::PUJetID);
    bcsv_1 = bjets.at(0).getProperty(PropertyEnum::bCSVscore);
    if(bjets.size()>1){
      //trailing b-jet sorted by pt
      bpt_2 = bjets.at(1).getP4().Pt();
      beta_2 = bjets.at(1).getP4().Eta();
      bphi_2 = bjets.at(1).getP4().Phi();
      brawf_2 = bjets.at(1).getProperty(PropertyEnum::rawPt)/bjets.at(1).getP4().Pt();
      bmva_2 = bjets.at(1).getProperty(PropertyEnum::PUJetID);
      bcsv_2 = bjets.at(1).getProperty(PropertyEnum::bCSVscore);
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::fillVetoes(const HTTEvent &event){

  dilepton_veto = event.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
  extraelec_veto = event.checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
  extramuon_veto = event.checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTSynchNTuple::selectEvent(const HTTEvent &event, HTTPair &pair){
  //MET filters
  int metFilters = (int)event.getMETFilterDecision();
  bool metFiltersPass = true;
  for(int i=0;i<6;++i){//0-5 classic filters
    metFiltersPass &= ( (metFilters & 1<<i) == 1<<i );
    //std::cout<<i<<". flt: "<<( (metFilters & 1<<i) == 1<<i )<<", all: "<<metFiltersPass<<std::endl;
  }
  for(int i=6;i<8;++i){//6-badChCand, 7-badPFMuon
    metFiltersPass &= ( (metFilters & 1<<i) == 1<<i );
    //std::cout<<i<<". flt: "<<!( (metFilters & 1<<i) == 1<<i )<<", all: "<<metFiltersPass<<std::endl;
  }
  metFiltersPass = true; //AP 08.06.2017 there are no MET filters in SM HTT
  if(!metFiltersPass) return false;

  //Add selection to tighten what is in ntuples
  if(decayMode_=="MuTau"){
    if( !(std::abs(pair.getLeg1().getP4().Eta())<2.4) )
      return false;
    if( !(pair.getLeg1().getP4().Pt()>23) )
      return false;
    if( !(std::abs(pair.getLeg2().getP4().Eta())<2.3) )
      return false;
    if( !(pair.getLeg2().getP4().Pt()>20) )
      return false;
  }
  else if(decayMode_=="TauTau"){
    if( !(std::abs(pair.getLeg1().getP4().Eta())<2.1 && std::abs(pair.getLeg2().getP4().Eta())<2.1) )
      return false;
    if( !(pair.getLeg1().getP4().Pt()>40 && pair.getLeg2().getP4().Pt()>40) )
      return false;
  }
  //
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::initializeCorrections(){
  //MB: basically a copy of ChannelSpecifics::initializeCorrections()

  std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v16_4.root";
  TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");

  RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");

  RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
  RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_binned_ratio");
  RooAbsReal *muon_trg_scalefactor = scaleWorkspace->function("m_trgOR4_binned_ratio");//MB 24->22
  RooAbsReal *muon_xtrg_scalefactor_iso  = scaleWorkspace->function("m_trgMu19leg_eta2p1_desy_ratio");//mu-iso<0.15?
  RooAbsReal *muon_xtrg_scalefactor_aiso = scaleWorkspace->function("m_trgMu19leg_eta2p1_aiso0p15to0p3_desy_ratio");//0.15<mu-iso<0.3?
  RooAbsReal *muon_trk_scalefactor = scaleWorkspace->function("m_trk_ratio");//MB not in HTTAnalysis
  RooAbsReal *tau_trg_genuine_efficiency = scaleWorkspace->function("t_genuine_TightIso_tt_ratio");//MB data->ratio
  RooAbsReal *tau_trg_fake_efficiency = scaleWorkspace->function("t_fake_TightIso_tt_ratio");//MB data->ratio

  h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
								     *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
								     RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
								     RooFit::Extended(kFALSE),
								     RooFit::Scaling(kFALSE));
  //
  RooArgSet dependentVarsForMu(*scaleWorkspace->var("m_pt"),*scaleWorkspace->var("m_eta"),*scaleWorkspace->var("m_iso"));
  RooArgSet projectedVarsForMu;
  const RooAbsReal *muon_iso_scalefactor_proj = muon_iso_scalefactor->createPlotProjection(dependentVarsForMu,projectedVarsForMu);
  h3DMuonIsoCorrections = (TH3F*)muon_iso_scalefactor_proj->createHistogram("h3DMuonIsoCorrections",
								       *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
								       RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),
									    RooFit::ZVar(*scaleWorkspace->var("m_iso"),RooFit::Binning(12,-0.05,0.55)),
											RooFit::Extended(kFALSE),
											RooFit::Scaling(kFALSE));
  const RooAbsReal *muon_trg_scalefactor_proj = muon_trg_scalefactor->createPlotProjection(dependentVarsForMu,projectedVarsForMu);
  h3DMuonTrgCorrections = (TH3F*)muon_trg_scalefactor_proj->createHistogram("h3DMuonTrgCorrections",
									    *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
									    RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),
									    RooFit::ZVar(*scaleWorkspace->var("m_iso"),RooFit::Binning(12,-0.05,0.55)),
									    RooFit::Extended(kFALSE),
									    RooFit::Scaling(kFALSE));
  //Uwaga: ponizsze korekcje sa 2d, do wlozenia w 3d
  TH2F *h2DMuonXTrgCorrections_iso = (TH2F*)muon_xtrg_scalefactor_iso->createHistogram("h2DMuonXTrgCorrections_iso",
										       *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
										       RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
										       RooFit::Extended(kFALSE),
										       RooFit::Scaling(kFALSE));
  TH2F *h2DMuonXTrgCorrections_aiso = (TH2F*)muon_xtrg_scalefactor_iso->createHistogram("h2DMuonXTrgCorrections_aiso",
										       *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
										       RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
										       RooFit::Extended(kFALSE),
										       RooFit::Scaling(kFALSE));
  h3DMuonXTrgCorrections = new TH3F("h3DMuonXTrgCorrections","",1980,10,1000,48,-2.4,2.4,8,-0.05,0.35);//less iso-bins than for id/single-trg due to different range
  for(int iX=0; iX<1980; ++iX){
    for(int iY=0; iY<48; ++iY){
      float binContent = h2DMuonXTrgCorrections_iso->GetBinContent(iX+1,iY+1);
      for(int iZ=0; iZ<4; ++iZ)//bins 1-4 correspond with iso<0.15
	h3DMuonXTrgCorrections->SetBinContent(iX+1,iY+1,iZ+1,binContent);
      binContent = h2DMuonXTrgCorrections_aiso->GetBinContent(iX+1,iY+1);
      for(int iZ=4; iZ<8; ++iZ)//bins 4-8 correspond with 0.15<iso<0.3
	h3DMuonXTrgCorrections->SetBinContent(iX+1,iY+1,iZ+1,binContent);
    }
  }
  h1DMuonTrkCorrections = (TH1F*)muon_trk_scalefactor->createHistogram("h1DMuonTrkCorrections",
								       *scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4),//MB m_abs_eta->m_eta
								       RooFit::Extended(kFALSE),
								       RooFit::Scaling(kFALSE));

  ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
  RooArgSet dependentVars(*scaleWorkspace->var("t_pt"),*scaleWorkspace->var("t_dm"));
  RooArgSet projectedVars;

  const RooAbsReal * tau_trg_genuine_efficiency_proj = tau_trg_genuine_efficiency->createPlotProjection(dependentVars,projectedVars);
  const RooAbsReal * tau_trg_fake_efficiency_proj = tau_trg_fake_efficiency->createPlotProjection(dependentVars,projectedVars);
  /*
  RooBinning binsForTauTrg(0,1000);
  binsForTauTrg.addUniform(1, 0, 30);
  binsForTauTrg.addUniform(2000, 30, 130);
  binsForTauTrg.addUniform(1000, 130, 330);
  binsForTauTrg.addUniform(1340, 330, 1000);
  */
  h2DTauTrgGenuineCorrections = (TH2F*)tau_trg_genuine_efficiency_proj->createHistogram("h2DTauTrgGenuineCorrections",
											*scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
											//*scaleWorkspace->var("t_pt"),RooFit::Binning(binsForTauTrg),
											RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(11,-0.5,10.5)),//MB proper binning (3->11) for DMs==0,1(2),10
											RooFit::Extended(kFALSE),
											RooFit::Scaling(kFALSE));

  h2DTauTrgFakeCorrections = (TH2F*)tau_trg_fake_efficiency_proj->createHistogram("h2DTauTrgFakeCorrections",
										  *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
										  //*scaleWorkspace->var("t_pt"),RooFit::Binning(binsForTauTrg),
										  RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(11,-0.5,10.5)),//MB proper binning (3->11) for DMs==0,1(2),10
										  RooFit::Extended(kFALSE),
										  RooFit::Scaling(kFALSE));

  delete aFile;
  if(h2DMuonXTrgCorrections_iso) delete h2DMuonXTrgCorrections_iso;
  if(h2DMuonXTrgCorrections_aiso) delete h2DMuonXTrgCorrections_aiso;

  ////
  //Open different RooWorkspace for missig SFs - ugly:(
				   
  std::string correctionFileName_b = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_sm_moriond_v2.root";
  TFile *bFile = TFile::Open(correctionFileName_b.c_str(),"CACHEREAD");
  
  RooWorkspace *scaleWorkspace_b = (RooWorkspace*)bFile->Get("w");

  RooAbsReal *tau_xtrg_genuine_efficiency = scaleWorkspace_b->function("t_genuine_TightIso_mt_ratio");//MB data->ratio
  RooAbsReal *tau_xtrg_fake_efficiency = scaleWorkspace_b->function("t_fake_TightIso_mt_ratio");//MB data->ratio

  ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
  RooArgSet dependentVarsForMT(*scaleWorkspace_b->var("t_pt"),*scaleWorkspace_b->var("t_eta"));
  RooArgSet projectedVarsForMT;

  const RooAbsReal * tau_xtrg_genuine_efficiency_proj = tau_xtrg_genuine_efficiency->createPlotProjection(dependentVarsForMT,projectedVarsForMT);
  const RooAbsReal * tau_xtrg_fake_efficiency_proj = tau_xtrg_fake_efficiency->createPlotProjection(dependentVarsForMT,projectedVarsForMT);

  h2DTauXTrgGenuineCorrections = (TH2F*)tau_xtrg_genuine_efficiency_proj->createHistogram("h2DTauXTrgGenuineCorrections",
											  *scaleWorkspace_b->var("t_pt"),RooFit::Binning(5000,0,1000),
											  RooFit::YVar(*scaleWorkspace_b->var("t_eta"),RooFit::Binning(4,-2.5,2.5)),//MB binning in eta: barrel/endcaps, i.e. <1.5
											  RooFit::Extended(kFALSE),
											  RooFit::Scaling(kFALSE));

  h2DTauXTrgFakeCorrections = (TH2F*)tau_xtrg_fake_efficiency_proj->createHistogram("h2DTauXTrgFakeCorrections",
										    *scaleWorkspace_b->var("t_pt"),RooFit::Binning(5000,0,1000),
										    RooFit::YVar(*scaleWorkspace_b->var("t_eta"),RooFit::Binning(4,-2.5,2.5)),//MB binning in eta: barrel/endcaps, i.e. <1.5
										    RooFit::Extended(kFALSE),
										    RooFit::Scaling(kFALSE));

  delete bFile;

  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTSynchNTuple::getPUWeight(float nPU){
  //MB: basically a copy of HTTAnalyzer::getPUWeight

  if(nPU<0) return 1.0;//not defined, e.g. for data

  if(!puDataFile_ || !puMCFile_ ||
     puDataFile_->IsZombie() ||
     puMCFile_->IsZombie()) { return 1.0; }

  if(!hPUVec_.size()) hPUVec_.resize(1);

  if(!hPUVec_[0]) {
    std::string hName = "pileup";
    TH1F *hPUData = (TH1F*)puDataFile_->Get(hName.c_str());
    TH1F *hPUSample = (TH1F*)puMCFile_->Get(hName.c_str());
    ///Normalise both histograms.
    hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
    ///
    hPUData->SetDirectory(0);
    hPUSample->SetDirectory(0);
    hPUData->SetName("h1DPUWeightSynch");
    hPUVec_[0] =  hPUData;
  }

  int iBinPU = hPUVec_[0]->FindBin(nPU);
  return hPUVec_[0]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::initializeBTagCorrections(){

  std::string csvFileName =  "CSVv2_Moriond17_B_H.csv";
  std::string weightsFileName = "http://akalinow.web.cern.ch/akalinow/"+csvFileName;
  TFile::Open(weightsFileName.c_str(),"CACHEREAD");
  std::string userName(std::getenv("USER"));
  std::string correctionFileName = "/tmp/"+userName+"/akalinow/"+csvFileName;

  calib = new BTagCalibration("CSVv2", correctionFileName);
  reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
				     "central",             // central sys type
				     {"up", "down"});       // other sys types

  reader->load(*calib,               // calibration instance
	       BTagEntry::FLAV_B,    // btag flavour
	       "comb");              // measurement type
  reader->load(*calib,               // calibration instance
	       BTagEntry::FLAV_C,    // btag flavour
	       "comb");              // measurement type
  reader->load(*calib,               // calibration instance
	       BTagEntry::FLAV_UDSG, // btag flavour
	       "incl");              // measurement type

  std::string efficiencyFileName = "http://akalinow.web.cern.ch/akalinow/tagging_efficiencies_Moriond2017.root";
  btagEffFile_ = TFile::Open(efficiencyFileName.c_str(),"CACHEREAD");

  btag_eff_b_ = (TH2F*)btagEffFile_->Get("btag_eff_b")->Clone("btag_eff_b");
  btag_eff_c_ = (TH2F*)btagEffFile_->Get("btag_eff_c")->Clone("btag_eff_c");
  btag_eff_oth_ = (TH2F*)btagEffFile_->Get("btag_eff_oth")->Clone("btag_eff_oth");

  rand_ = new TRandom3();

  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTSynchNTuple::promoteBJet(const HTTParticle &jet){
  //MB: https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Standalone
  bool decision = false;

  if(!reader) initializeBTagCorrections();

  BTagEntry::JetFlavor jetFlavour;
  if(std::abs(jet.getProperty(PropertyEnum::Flavour))==5)//b-quark
    jetFlavour = BTagEntry::FLAV_B;
  else if(std::abs(jet.getProperty(PropertyEnum::Flavour))==4)//c-quark
    jetFlavour = BTagEntry::FLAV_C;
  else //light quark, gluon or undefined
    jetFlavour = BTagEntry::FLAV_UDSG;
  // Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
  double btag_SF = reader->eval_auto_bounds("central",
					    jetFlavour,
					    jet.getP4().Eta(),
					    jet.getP4().Pt()
					    //,jet.getProperty(PropertyEnum::bCSVscore) //MB: it is not needed when WP is definied
					    );
  rand_->SetSeed((int)((jet.getP4().Eta()+5)*100000));
  double rand_num = rand_->Rndm();
  //debug
  //std::cout<<"\tbtag_SF(flav,CSVv2): "<<btag_SF
  //	   <<"("<<jetFlavour<<","
  //	   <<jet.getProperty(PropertyEnum::bCSVscore)<<")"<<std::endl;
  //std::cout<<"\tbtag_rand_num: "<<rand_num<<std::endl;
  if(btag_SF>1){
    double tagging_efficiency = 1;
    TH2F *histo_eff = btag_eff_oth_;
    if(jetFlavour == BTagEntry::FLAV_B)
      histo_eff = btag_eff_b_;
    else if(jetFlavour == BTagEntry::FLAV_C)
      histo_eff = btag_eff_c_;
    if( jet.getP4().Pt() > histo_eff->GetXaxis()->GetBinLowEdge(histo_eff->GetNbinsX()+1) ){
      tagging_efficiency = histo_eff->GetBinContent( histo_eff->GetNbinsX(),histo_eff->GetYaxis()->FindBin(std::abs(jet.getP4().Eta())) );
    }
    else{
      tagging_efficiency = histo_eff->GetBinContent( histo_eff->GetXaxis()->FindBin(jet.getP4().Pt()),histo_eff->GetYaxis()->FindBin(std::abs(jet.getP4().Eta())) );
    }
    //debug
    //std::cout<<"\tbtag_eff: "<<tagging_efficiency<<std::endl;
    if(tagging_efficiency < 1e-9)//protection
      decision = false;
    else if(tagging_efficiency > 1.-1e-9)//protection
      decision = true;
    else
      decision = (rand_num < (1. - btag_SF)/(1. - 1./tagging_efficiency) );
  }
  else{
    decision = (rand_num < 1. - btag_SF);
  }
  //debug
  //std::cout<<"\tbtag_decision: "<<decision<<std::endl;
  return !decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
