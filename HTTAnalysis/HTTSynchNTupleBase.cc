#include <sstream>

#include "HTTSynchNTupleBase.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTupleBase::HTTSynchNTupleBase(const std::string & aName):Analyzer(aName){ tmpName = "h1DXSignal";}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTupleBase::~HTTSynchNTupleBase(){ if(myHistos_) delete myHistos_; }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTupleBase::clone() const{

  HTTSynchNTupleBase* clone = new HTTSynchNTupleBase(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::initialize(TDirectory* aDir,
			      pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "HTTSynchNTupleBase"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::clearTTreeVariables(){ 
  
  //event ID variables
  run = -999;
  lumi = -999;
  evt = -999;
  //Pielup
  npv = -999;
  npu = -999;
  rho = -999;
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
  decayModeFindingOldDMs_1 = -999;
  neutralIsoPtSum_1 = -999;
  puCorrPtSum_1 = -999;
  trigweight_1 = -999;
  idisoweight_1 = -999;
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
  decayModeFindingOldDMs_2 = -999;
  neutralIsoPtSum_2 = -999;
  puCorrPtSum_2 = -999;
  trigweight_2 = -999;
  idisoweight_2 = -999;
  //di-tau system
  pt_tt = -999;
  mt_tot = -999;
  m_vis = -999;
  m_sv = -999;
  mt_sv = -999;
  //MET
  met = -999;
  metphi = -999;
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
  dilepton_veto = -999;
  extraelec_veto = -999;
  extramuon_veto = -999;
  puweight = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::addBranch(TTree *tree){

  //event ID variables
  tree->Branch("run",&run);
  tree->Branch("lumi",&lumi);
  tree->Branch("evt",&evt);
  //Pielup
  tree->Branch("npv",&npv);
  tree->Branch("npu",&npu);
  tree->Branch("rho",&rho);
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
  //di-tau system
  tree->Branch("pt_tt",&pt_tt);
  tree->Branch("mt_tot",&mt_tot);
  tree->Branch("m_vis",&m_vis);
  tree->Branch("m_sv",&m_sv);
  tree->Branch("mt_sv",&mt_sv);
  //MET
  tree->Branch("met",&met);
  tree->Branch("metphi",&metphi);
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
  tree->Branch("puweight",&puweight);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTSynchNTupleBase::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  HTTEvent &aEvent = *myEventProxy.event;  
  HTTPair &aPair = (*myEventProxy.pairs)[0];
  std::vector<HTTParticle> &aJets = *myEventProxy.jets;
  
  clearTTreeVariables();

  i_++;
  if(i_%1000==0){
    //std::cout<<i_<<std::endl;
  }
	
  ///Filling TTree
  //event ID variables
  fillEventID(aEvent);
  // Fill legs
  fillLegs(aPair.getLeg1(), aPair.getLeg2());
  // Fill di-tau system, but also MET which can be pairwise
  fillPair(aEvent, aPair);
  // Fill jet variables includig VBF system
  fillJets(aJets);
  //Fill extra lepton vetoes info
  fillVetoes(aEvent);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::fillEventID(const HTTEvent &event){

  run = event.getRunId();
  lumi = 0; //NEED TO FIX, not in ntuples
  evt = event.getEventId();
  //Pielup
  npv = event.getNPV();
  npu = event.getNPU();
  rho = 0; //NEED TO FIX, not in ntuples
  //puweight;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::fillLegs(const HTTParticle &leg1, const HTTParticle &leg2){

  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  pt_1 = leg1.getP4().Pt();
  phi_1 = leg1.getP4().Phi();
  eta_1 = leg1.getP4().Eta();
  m_1 = leg1.getP4().M();
  q_1 = leg1.getCharge();
  d0_1 = leg1.getProperty(PropertyEnum::dxy);
  dZ_1 = leg1.getProperty(PropertyEnum::dz);
  gen_match_1 = leg1.getProperty(PropertyEnum::mc_match); 

  //Leg 2 (trailing tau for tt, electon for et,em muon for mt)
  pt_2 = leg2.getP4().Pt();
  phi_2 = leg2.getP4().Phi();
  eta_2 = leg2.getP4().Eta();
  m_2 = leg2.getP4().M();
  q_2 = leg2.getCharge();
  d0_2 = leg2.getProperty(PropertyEnum::dxy);
  dZ_2 = leg2.getProperty(PropertyEnum::dz);
  gen_match_2 = leg2.getProperty(PropertyEnum::mc_match); //according to: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching

  //Decay channel specific
  fillLegsSpecific(leg1,leg2);
 }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2){
  //Do nothing, to be implemented for derivied class for specific decay channels
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::fillPair(const HTTEvent &event, HTTPair &pair){

  //Legs
  TLorentzVector leg1P4 = pair.getLeg1().getP4();
  TLorentzVector leg2P4 = pair.getLeg2().getP4();
	
  //MET
  TLorentzVector metP4(event.getMET().X(), event.getMET().Y(), 0, event.getMET().Mod());
  met = event.getMET().Mod();
  metphi = event.getMET().Phi_mpi_pi(event.getMET().Phi());
  pfmt_1 = TMath::Sqrt(2.*leg1P4.Pt()*metP4.Pt()*(1.-TMath::Cos(leg1P4.Phi()-metP4.Phi()))); 
  pfmt_2 = TMath::Sqrt(2.*leg2P4.Pt()*metP4.Pt()*(1.-TMath::Cos(leg2P4.Phi()-metP4.Phi()))); 
  TLorentzVector mvametP4(pair.getMET().X(), pair.getMET().Y(), 0, pair.getMET().Mod());
  mvamet = pair.getMET().Mod();
  mvametphi = pair.getMET().Phi_mpi_pi(pair.getMET().Phi());
  mt_1 = pair.getMTLeg1();
  mt_2 = pair.getMTLeg2();
  /*
    puppimet;
    puppimetphi;
    puppimt_1;
    puppimt_2;
  */

  //di-tau system
  pt_tt = (mvametP4 + leg1P4 + leg2P4).Pt();
  mt_tot = TMath::Sqrt( 2. * leg1P4.Pt() * mvametP4.Pt() * ( 1. - TMath::Cos( leg1P4.Phi() - mvametP4.Phi() ) ) + 2. * leg2P4.Pt() * mvametP4.Pt() * (1. - TMath::Cos( leg2P4.Phi() - mvametP4.Phi() )) + 2. * leg1P4.Pt() * leg2P4.Pt() * TMath::Cos(leg1P4.Phi() - leg2P4.Phi()) );
  m_vis = (leg1P4 + leg2P4).M();
  m_sv = pair.getP4SVFit().M();

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
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTupleBase::fillJets(const std::vector<HTTParticle> &jets){
  
  njetspt20 = jets.size();
  njets = 0;
  std::vector<HTTParticle> bjets;
  for(unsigned int iJet=0; iJet<jets.size(); ++iJet){
    if(jets.at(iJet).getP4().Pt()>30)
      njets++;
    if(std::abs(jets.at(iJet).getP4().Eta())>2.4 && 
       jets.at(iJet).getProperty(PropertyEnum::bCSVscore)>0.8){//FIXME Correct??
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
void HTTSynchNTupleBase::fillVetoes(const HTTEvent &event){

  dilepton_veto = event.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
  extraelec_veto = event.checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
  extramuon_veto = event.checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
