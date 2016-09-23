#include <sstream>

#include "HTTSynchNTuple.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::HTTSynchNTuple(const std::string & aName):Analyzer(aName){ tmpName = "h1DXSignal";}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::~HTTSynchNTuple(){ if(myHistos_) delete myHistos_; }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTuple::clone() const{

  HTTSynchNTuple* clone = new HTTSynchNTuple(name());
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
  myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::clearTTreeVariables(){ 

  
	//event ID variables
	run = -999;
	lumi = -999;
	evt = -999;
	//Pielup
	npv = -999;
	npu = -999;
	rho = -999;
	//Leg 1 (leading tau for tt, electon for et,em muon for mt) - here Leg 1 is always muon!
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
	//Leg 2 (trailing tau for tt, electon for et,em muon for mt)
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
	
	aLeadingJet.clear();
	aTrailingJet.clear();
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::addBranch(TTree *tree){

	//event ID variables
	tree->Branch("run",&run);
	tree->Branch("lumi",&lumi);
	tree->Branch("evt",&evt);
	//Pielup
	tree->Branch("npv",&npv);
	tree->Branch("npu",&npu);
	tree->Branch("rho",&rho);
	//Leg 1 (leading tau for tt, electon for et,em muon for mt)
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
	//Leg 2 (trailing tau for tt, electon for et,em muon for mt)
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
bool HTTSynchNTuple::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  aEvent = *myEventProxy.event;  
  aPair = (*myEventProxy.pairs)[0];
  aTau = aPair.getTau();
  aMuon = aPair.getMuon();
  aJets = myEventProxy.jets;
  
  clearTTreeVariables();
  if(aJets->size()!=0) {aLeadingJet = aJets->at(0);} 
  if(aJets->size()>1) {aTrailingJet = aJets->at(1);}

	i_++;
	if(i_%1000==0){
		//std::cout<<i_<<std::endl;
		}
	
	//Filling TTree
	//event ID variables
	run = aEvent.getRunId();
	lumi = 0;						//NEED TO FIX, not in ntuples
	evt = aEvent.getEventId();
	//Pielup
	npv = aEvent.getNPV();
	npu = aEvent.getNPU();
	rho = 0;						//NEED TO FIX, not in ntuples
	//Leg 1 (leading tau for tt, electon for et,em muon for mt) - here Leg 1 is always muon!
	leg1 = aMuon;
	pt_1 = leg1.getP4().Perp();
	phi_1 = leg1.getP4().Phi();
	eta_1 = leg1.getP4().Eta();
	m_1 = leg1.getP4().M();
	q_1 = leg1.getCharge();
	d0_1 = leg1.getProperty(PropertyEnum::dxy);
	dZ_1 = leg1.getProperty(PropertyEnum::dz);
	mt_1 = aPair.getMTMuon();/*
	pfmt_1;
	puppimt_1;*/
	iso_1 = leg1.getProperty(PropertyEnum::combreliso);/*
	id_e_mva_nt_loose_1;*/
	gen_match_1 = leg2.getProperty(PropertyEnum::PDGId);	/*	
	againstElectronLooseMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronLooseMVA6);				//FIX: following 7 properites need to be added
	againstElectronMediumMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronMediumMVA6);
	againstElectronTightMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronTightMVA6);
	againstElectronVLooseMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronVLooseMVA6);
	againstElectronVTightMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronVTightMVA6);
	againstMuonLoose3_1 = leg1.getProperty(PropertyEnum::againstMuonLoose3);
	againstMuonTight3_1 = leg1.getProperty(PropertyEnum::againstMuonLoose3);*/
	byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = leg1.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);/*
	byIsolationMVA3newDMwoLTraw_1;
	byIsolationMVA3oldDMwoLTraw_1;
	byIsolationMVA3newDMwLTraw_1;
	byIsolationMVA3oldDMwLTraw_1;
	chargedIsoPtSum_1;
	decayModeFindingOldDMs_1;
	neutralIsoPtSum_1;
	puCorrPtSum_1;
	trigweight_1;
	idisoweight_1;*/
	//Leg 2 (trailing tau for tt, electon for et,em muon for mt)
	leg2 = aTau;
	pt_2 = leg2.getP4().Perp();
	phi_2 = leg2.getP4().Phi();
	eta_2 = leg2.getP4().Eta();
	m_2 = leg2.getP4().M();
	q_2 = leg2.getCharge();
	d0_2 = leg2.getProperty(PropertyEnum::dxy);
	dZ_2 = leg2.getProperty(PropertyEnum::dz);
	mt_2 = abs(leg2.getPDGid())==15 ? aPair.getMTLeg2() : aPair.getMTLeg1();/*
	pfmt_2;
	puppimt_2;*/
	iso_2 = leg2.getProperty(PropertyEnum::combreliso);/*
	id_e_mva_nt_loose_2;*/
	gen_match_2 = leg2.getProperty(PropertyEnum::PDGId);		/*//according to: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching
	againstElectronLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronLooseMVA6);				//FIX: following 7 properites need to be added
	againstElectronMediumMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronMediumMVA6);
	againstElectronTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronTightMVA6);
	againstElectronVLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVLooseMVA6);
	againstElectronVTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVTightMVA6);
	againstMuonLoose3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);
	againstMuonTight3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);*/
	byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);/*
	byIsolationMVA3newDMwoLTraw_2;
	byIsolationMVA3oldDMwoLTraw_2;
	byIsolationMVA3newDMwLTraw_2;
	byIsolationMVA3oldDMwLTraw_2;
	chargedIsoPtSum_2;
	decayModeFindingOldDMs_2;
	neutralIsoPtSum_2;
	puCorrPtSum_2;
	trigweight_2;
	idisoweight_2;*/
	
	//di-tau system
	TLorentzVector metP4(aEvent.getMET().X(), aEvent.getMET().Y(), 0, aEvent.getMET().Mod());
	TLorentzVector mvametP4(aPair.getMET().X(), aPair.getMET().Y(), 0, aPair.getMET().Mod());
	pt_tt = (metP4 + aTau.getP4() + aMuon.getP4()).Perp();
	mt_tot = TMath::Sqrt( 2 * pt_1 * mvametP4.Perp() * ( 1 - TMath::Cos( phi_1 - mvametP4.Phi() ) ) + 2 * pt_2 * mvametP4.Perp() * (1 - TMath::Cos( phi_2 - mvametP4.Phi() )) + 2 * pt_1 * pt_2 * TMath::Cos(phi_1 - phi_2) );
	m_vis = (aMuon.getP4() + aTau.getP4()).M();
	m_sv = aPair.getP4SVFit().M();
	//mt_sv;
	
	//MET
	met = aEvent.getMET().Mod();
	metphi = aEvent.getMET().Phi_mpi_pi(aEvent.getMET().Phi());
/*
	puppimet;
	puppimetphi;
*/
	mvamet = aPair.getMET().Mod();
	mvametphi = aPair.getMET().Phi_mpi_pi(aPair.getMET().Phi());
	float muPx = aMuon.getP4().X(), muPy = aMuon.getP4().Y(), tauPx = aTau.getP4().X(), tauPy = aTau.getP4().Y();
	pzetavis = ( (muPx+tauPx)*(TMath::Cos(phi_1) + TMath::Cos(phi_2)) + (muPy + tauPy)*(TMath::Sin(phi_1) + TMath::Sin(phi_2)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(phi_1) + TMath::Cos(phi_2),2) + TMath::Power(TMath::Sin(phi_1) + TMath::Sin(phi_2),2) ) );
	pzetamiss = ( aPair.getMET().X()*(TMath::Cos(phi_1) + TMath::Cos(phi_2)) +  aPair.getMET().Y()*(TMath::Sin(phi_1) + TMath::Sin(phi_2)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(phi_1) + TMath::Cos(phi_2),2) + TMath::Power(TMath::Sin(phi_1) + TMath::Sin(phi_2),2) ) );
	pfpzetamiss = ( aEvent.getMET().X()*(TMath::Cos(phi_1) + TMath::Cos(phi_2)) +  aEvent.getMET().Y()*(TMath::Sin(phi_1) + TMath::Sin(phi_2)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(phi_1) + TMath::Cos(phi_2),2) + TMath::Power(TMath::Sin(phi_1) + TMath::Sin(phi_2),2) ) );
	//puppipzetamiss;
	mvacov00 = aPair.getMETMatrix().at(0);
	mvacov01 = aPair.getMETMatrix().at(1);
	mvacov10 = aPair.getMETMatrix().at(2);
	mvacov11 = aPair.getMETMatrix().at(3);/*
	metcov00;
	metcov01;
	metcov10;
	metcov11;*/
	
	//VBF system
	if(aJets->size()>=2){
		mjj = (aLeadingJet.getP4()+aTrailingJet.getP4()).M();
		jdeta = std::abs(aLeadingJet.getP4().Eta()-aTrailingJet.getP4().Eta());/*
		njetingap;
		njetingap20;*/
		jdphi = aLeadingJet.getP4().Phi()-aTrailingJet.getP4().Phi();
		while (jdphi > TMath:Pi() ) jdphi -= 2.*TMath:Pi();
	        while (jdphi <= -TMath:Pi() ) jdphi += 2.*TMath:Pi();
		}/*
		
	//additional jets
	nbtag;
	njets;
	njetspt20;*/
	
	//leading jet sorted by pt
	if(aJets->size()>0){
		jpt_1 = aLeadingJet.getP4().Perp();
		//std::cout<<aLeadingJet.getP4().Perp()<<"\n";
		jeta_1 = aLeadingJet.getP4().Eta();
		jphi_1 = aLeadingJet.getP4().Phi();
//		jrawf_1;
		jmva_1 = aLeadingJet.getProperty(PropertyEnum::PUJetID);
		}
		
	//trailing jet sorted by pt
	if(aJets->size()>1){
		jpt_2 = aTrailingJet.getP4().Perp();
		jeta_2 = aTrailingJet.getP4().Eta();
		jphi_2 = aTrailingJet.getP4().Phi();
		//jrawf_2;
		jmva_2 = aTrailingJet.getProperty(PropertyEnum::PUJetID);
		}/*
	//leading b-jet sorted by pt
	bpt_1;
	beta_1;
	bphi_1;
	brawf_1;
	bmva_1;
	bcsv_1;
	//trailing b-jet sorted by pt
	bpt_2;
	beta_2;
	bphi_2;
	brawf_2;
	bmva_2;
	bcsv_2;
	//Extra lepton vetos
	dilepton_veto = aEvent.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
/*	extraelec_veto;
	extramuon_veto;
	puweight;*/

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

