#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName),nPCAMin_(0.003){

  //pileupCalc.py -i lumiSummary_Run2016BCDE_PromptReco_v12.json
  //--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
  //--calcMode true --minBiasXsec 69200 --maxPileupBin 60 --numPileupBins 600 Data_Pileup_Cert_271036-277148.root

  ///Load ROOT file with PU histograms.
  std::string filePath = "Data_Pileup_2016_BCDEFG_v26.root";
  //filePath = "Data_Pileup_2016_July22.root";
  puDataFile_ = new TFile(filePath.c_str());

  filePath = "MC_Spring16_PU25ns_V1.root";
  puMCFile_ = new TFile(filePath.c_str());

  ntupleFile_ = 0;
  hStatsFromFile = 0;

  h2DMuonIdCorrections = 0;
  h2DMuonIsoCorrections = 0;
  h2DMuonTrgCorrections = 0;
  h3DTauCorrections = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
  if(puDataFile_) delete puDataFile_;
  if(puMCFile_) delete puMCFile_;

if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
if(h2DMuonIsoCorrections) delete h2DMuonIsoCorrections;
if(h2DMuonTrgCorrections) delete h2DMuonTrgCorrections;
if(h3DTauCorrections) delete h3DTauCorrections;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTAnalyzer::clone() const{

  HTTAnalyzer* clone = new HTTAnalyzer(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initialize(TDirectory* aDir,
			     pat::strbitset *aSelections){

  mySelections_ = aSelections;

  myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initializeCorrections(){

#pragma omp critical
{
  std::string filePath = "htt_scalefactors_v5.root";
  TFile aFile(filePath.c_str());
  RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile.Get("w");

  RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
  RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio");
  RooAbsReal *muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data");//OR of the HLT_IsoMu22 and HLT_IsoTkMu22
  //RooAbsReal *tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_m_pt30_sf");

  h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
  *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
  RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
  RooFit::Scaling(kFALSE));

  h2DMuonIsoCorrections = (TH2F*)muon_iso_scalefactor->createHistogram("h2DMuonIsoCorrections",
  *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
  RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
  RooFit::Scaling(kFALSE));

  h2DMuonTrgCorrections = (TH2F*)muon_trg_efficiency->createHistogram("h2DMuonTrgCorrections",
  *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
  RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
  RooFit::Scaling(kFALSE));
}
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){ myHistos_->finalizeHistograms(0,1.0); }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> HTTAnalyzer::getSeparatedJets(const EventProxyHTT & myEventProxy,
						float deltaR){

  std::vector<HTTParticle> separatedJets;

  for(auto aJet: *myEventProxy.jets){
    float dRTau = aJet.getP4().DeltaR(aTau.getP4());
    float dRMu = aJet.getP4().DeltaR(aMuon.getP4());
    bool loosePFJetID = aJet.getProperty(PropertyEnum::PFjetID)>=1;
    if(dRTau>deltaR && dRMu>deltaR && loosePFJetID) separatedJets.push_back(aJet);
  }

  return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  aEvent = *myEventProxy.event;
  aPair = (*myEventProxy.pairs)[0];
  aTau = aPair.getTau();
  aMuon = aPair.getMuon();

  TLorentzVector met4v(aPair.getMET().X(),
		       aPair.getMET().Y(),
		       0,
		       aPair.getMET().Mod());

  aMET = HTTParticle();
  aMET.setP4(met4v);

  aGenMuonTau = HTTParticle();
  aGenHadTau = HTTParticle();

  if(myEventProxy.genLeptons && myEventProxy.genLeptons->size()){
    HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
    if(aGenTau.getProperty(PropertyEnum::decayMode)==tauDecayMuon) aGenMuonTau = aGenTau;
    else if(aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecaysElectron) aGenHadTau = aGenTau;
    if(myEventProxy.genLeptons->size()>1){
      aGenTau =  myEventProxy.genLeptons->at(1);
      if(aGenTau.getProperty(PropertyEnum::decayMode)==tauDecayMuon) aGenMuonTau = aGenTau;
      else if(aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecaysElectron) aGenHadTau = aGenTau;
    }
  }

  aSeparatedJets = getSeparatedJets(myEventProxy, 0.5);
  aJet1 = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
  aJet2 = aSeparatedJets.size()>1 ? aSeparatedJets[1] : HTTParticle();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<bool, bool> HTTAnalyzer::checkTauDecayMode(const EventProxyHTT & myEventProxy){

  bool goodGenDecayMode = false;
  bool goodRecoDecayMode = false;

  std::vector<std::string> decayNamesGen = getTauDecayName(aGenHadTau.getProperty(PropertyEnum::decayMode),
							   aGenMuonTau.getProperty(PropertyEnum::decayMode));
  std::vector<std::string> decayNamesReco = getTauDecayName(aTau.getProperty(PropertyEnum::decayMode),HTTAnalyzer::tauDecayMuon);

  for(auto it: decayNamesGen) if(it.find("Lepton1Prong")!=std::string::npos) goodGenDecayMode = true;
  for(auto it: decayNamesReco) if(it.find("Lepton1Prong")!=std::string::npos) goodRecoDecayMode = true;

  return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::addBranch(TTree *tree){/*tree->Branch("muonPt",&muonPt);*/}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillControlHistos(const std::string & hNameSuffix, float eventWeight,
				    const sysEffects::sysEffectsEnum & aSysEffect){

  ///Histograms filled for each systematic effect
  ///Fill SVfit and visible masses
  TLorentzVector aVisSum = aMuon.getP4(aSysEffect) + aTau.getP4(aSysEffect);
  float visMass = aVisSum.M();
  float higgsPt =  (aVisSum + aMET.getP4(aSysEffect)).Pt();
  float jetsMass = 0;
  if(nJets30>1) jetsMass = (aJet1.getP4(aSysEffect)+aJet2.getP4(aSysEffect)).M();

  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(aSysEffect),eventWeight);
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4(aSysEffect).M(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix, visMass, eventWeight);

  ///Unrolled distributions for 2D fit
  myHistos_->fill2DUnrolledHistogram("h1DUnRollTauPtMassVis"+hNameSuffix, visMass, aTau.getP4(aSysEffect).Pt(),eventWeight);
  myHistos_->fill2DUnrolledHistogram("h1DUnRollHiggsPtMassSV"+hNameSuffix, aPair.getP4(aSysEffect).M(), higgsPt, eventWeight);
  myHistos_->fill2DUnrolledHistogram("h1DUnRollMjjMassSV"+hNameSuffix, aPair.getP4(aSysEffect).M(), jetsMass, eventWeight);

  fillDecayPlaneAngle(hNameSuffix, eventWeight);
  if(aSysEffect!=sysEffects::NOMINAL_SVFIT) return;

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.getNPV(),eventWeight);

  ///Fill muon
  myHistos_->fill1DHistogram("h1DPtMuon"+hNameSuffix,aMuon.getP4(aSysEffect).Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+hNameSuffix,aMuon.getP4(aSysEffect).Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMuon"+hNameSuffix,aMuon.getP4(aSysEffect).Phi(),eventWeight);
  myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
  myHistos_->fill1DHistogram("h1DnPCAMuon"+hNameSuffix,aMuon.getPCARefitPV().Mag(),eventWeight);

  ///Fill tau
  myHistos_->fill1DHistogram("h1DPtTau"+hNameSuffix,aTau.getP4(aSysEffect).Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+hNameSuffix,aTau.getP4(aSysEffect).Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+hNameSuffix,aTau.getP4(aSysEffect).Phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DIDTau"+hNameSuffix,aTau.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw) ,eventWeight);
  myHistos_->fill1DHistogram("h1DStatsDecayMode"+hNameSuffix, aTau.getProperty(PropertyEnum::decayMode), eventWeight);
  myHistos_->fill1DHistogram("h1DnPCATau"+hNameSuffix,aTau.getPCARefitPV().Mag(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtTauLeadingTk"+hNameSuffix,aTau.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);
  myHistos_->fill1DHistogram("h1DPtMuTauMET"+hNameSuffix,higgsPt,eventWeight);

  ///Fill jets info
  myHistos_->fill1DHistogram("h1DStatsNJ30"+hNameSuffix,nJets30,eventWeight);
  if(nJets30>0){
    myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet1.getP4(aSysEffect).Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet1.getP4(aSysEffect).Eta(),eventWeight);
    myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet1.getProperty(PropertyEnum::bCSVscore),eventWeight);
  }
  if(nJets30>1) myHistos_->fill1DHistogram("h1DWideMass2J"+hNameSuffix,jetsMass,eventWeight);
  myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aMET.getP4(aSysEffect).Pt(),eventWeight);

  if(aJet1.getProperty(PropertyEnum::bCSVscore)>0.8){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet1.getP4(aSysEffect).Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet1.getP4(aSysEffect).Eta(),eventWeight);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::fillVertices(const std::string & sysType){

  TVector3 aVertexGen = aEvent.getGenPV();
  TVector3 aVertex = aEvent.getRefittedPV();
  if(sysType.find("AODPV")!=std::string::npos) aVertex = aEvent.getAODPV();
  if(sysType.find("RefitPV")!=std::string::npos) aVertex = aEvent.getRefittedPV();

  float pullX = aVertexGen.X() - aVertex.X();
  float pullY = aVertexGen.Y() - aVertex.Y();
  float pullZ = aVertexGen.Z() - aVertex.Z();

  myHistos_->fill1DHistogram("h1DVxPullX"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullY"+sysType,pullY);
  myHistos_->fill1DHistogram("h1DVxPullZ"+sysType,pullZ);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::passCategory(const HTTAnalyzer::muTauCategory & aCategory,
			       const sysEffects::sysEffectsEnum & aSysEffect){

  if(categoryDecisions.size()) categoryDecisions.clear();
  else categoryDecisions = std::vector<bool>((int)HTTAnalyzer::DUMMY);

  nJets30 = 0;
  for(auto itJet:aSeparatedJets){
    if(itJet.getP4(aSysEffect).Pt()>30) ++nJets30;
  }

  nJetsInGap30 = 0;
  if(nJets30>=2){
    for(unsigned int iJet=2; iJet<aSeparatedJets.size(); ++iJet){
      if( (aSeparatedJets.at(iJet).getP4().Eta()>aJet1.getP4().Eta()&&aSeparatedJets.at(iJet).getP4().Eta()<aJet2.getP4().Eta()) ||
          (aSeparatedJets.at(iJet).getP4().Eta()<aJet1.getP4().Eta()&&aSeparatedJets.at(iJet).getP4().Eta()>aJet2.getP4().Eta()) ){
        if(aSeparatedJets.at(iJet).getP4(aSysEffect).Pt()>30) nJetsInGap30++;
      }
    }
  }

  float jetsMass = (aJet1.getP4(aSysEffect)+aJet2.getP4(aSysEffect)).M();
  float higgsPt =  (aTau.getP4(aSysEffect) + aTau.getP4(aSysEffect) + aMET.getP4(aSysEffect)).Pt();
  bool mtSelection = aPair.getMTMuon(aSysEffect)<50;

  bool jet0_low =  aTau.getP4(aSysEffect).Pt()>20  && aTau.getP4(aSysEffect).Pt()<50 && nJets30==0;
  bool jet0_high = aTau.getP4(aSysEffect).Pt()>50 && nJets30==0;

  bool jet1_low = (nJets30==1 || (nJets30==2 && jetsMass<500)) &&
    (aTau.getP4(aSysEffect).Pt()>30 && aTau.getP4(aSysEffect).Pt()<40 ||
     aTau.getP4(aSysEffect).Pt()>40 && higgsPt<140);

  bool jet1_high = (nJets30==1 || (nJets30==2 && jetsMass<500)) &&
    (aTau.getP4(aSysEffect).Pt()>40 && higgsPt>140);

  bool vbf_low = aTau.getP4(aSysEffect).Pt()>20 &&
    nJets30==2 && jetsMass>500 &&
                 (jetsMass<800 || higgsPt<100);

  bool vbf_high = aTau.getP4(aSysEffect).Pt()>20 &&
    nJets30==2 && jetsMass>800 && higgsPt>100;

  bool cpMuonSelection = aMuon.getPCARefitPV().Perp()>nPCAMin_;
  bool cpTauSelection =  aTau.getPCARefitPV().Mag()>nPCAMin_;
  bool cpPi = cpMuonSelection && cpTauSelection && aTau.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero;
  bool cpRho = cpMuonSelection && aTau.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero &&
    isOneProng(aTau.getProperty(PropertyEnum::decayMode));

  //2D categories
  bool jet0 = aTau.getP4(aSysEffect).Perp()>30 && nJets30 == 0;
  bool boosted = aTau.getP4(aSysEffect).Perp()>30 && (nJets30==1 || (nJets30==2 && jetsMass < 300) || nJets30 > 2);
  bool vbf = aTau.getP4(aSysEffect).Perp()>30 && nJets30==2 && jetsMass>300;

  bool wSelection = aPair.getMTMuon(aSysEffect)>80 && aMuon.getProperty(PropertyEnum::combreliso)<0.15;
  bool ttSelection =  aPair.getMTMuon(aSysEffect)>150;
  bool antiiso = aMuon.getProperty(PropertyEnum::combreliso)>0.15 && aMuon.getProperty(PropertyEnum::combreliso)<0.30;

  categoryDecisions[(int)HTTAnalyzer::jet0_low] = mtSelection && jet0_low;
  categoryDecisions[(int)HTTAnalyzer::jet0_high] = mtSelection && jet0_high;

  categoryDecisions[(int)HTTAnalyzer::jet1_low] = mtSelection && jet1_low;
  categoryDecisions[(int)HTTAnalyzer::jet1_high] = mtSelection && jet1_high;

  categoryDecisions[(int)HTTAnalyzer::vbf_low] = mtSelection && vbf_low;
  categoryDecisions[(int)HTTAnalyzer::vbf_high] = mtSelection && vbf_high;

  categoryDecisions[(int)HTTAnalyzer::jet0] = mtSelection && jet0;
  categoryDecisions[(int)HTTAnalyzer::CP_Pi] = mtSelection && cpPi;
  categoryDecisions[(int)HTTAnalyzer::CP_Rho] = mtSelection && cpRho;
  categoryDecisions[(int)HTTAnalyzer::boosted] = mtSelection && boosted;
  categoryDecisions[(int)HTTAnalyzer::vbf] = mtSelection && vbf;

  categoryDecisions[(int)HTTAnalyzer::wjets_jet0] = wSelection && jet0;
  categoryDecisions[(int)HTTAnalyzer::wjets_boosted] = wSelection && boosted;
  categoryDecisions[(int)HTTAnalyzer::wjets_vbf] = wSelection && vbf;

  categoryDecisions[(int)HTTAnalyzer::antiiso_jet0] = antiiso && mtSelection && jet0;
  categoryDecisions[(int)HTTAnalyzer::antiiso_boosted] = antiiso && mtSelection && boosted;
  categoryDecisions[(int)HTTAnalyzer::antiiso_vbf] = antiiso && mtSelection && vbf;

  categoryDecisions[(int)HTTAnalyzer::W] = wSelection;
  categoryDecisions[(int)HTTAnalyzer::TT] = ttSelection;

  return categoryDecisions[(int)aCategory];
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  sampleName = getSampleName(myEventProxy);

  std::string hNameSuffix = sampleName;
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);
  float ptReweight = 1.0;
  if(sampleName.find("DYJets")!=std::string::npos ||
     sampleName.find("TTbar")!=std::string::npos)
    ptReweight = myEventProxy.event->getPtReWeight();
  float eventWeight = puWeight*genWeight*ptReweight;

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,0,eventWeight);
  myHistos_->fill1DHistogram("h1DNPartons"+hNameSuffix,myEventProxy.event->getLHEnOutPartons(),eventWeight);
  getPreselectionEff(myEventProxy);

  bool postSynchTau = aEvent.checkSelectionBit(SelectionBitsEnum::postSynchTau);
  bool postSynchMuon = aEvent.checkSelectionBit(SelectionBitsEnum::postSynchMuon);
  bool diMuonVeto = aEvent.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
  bool extraMuonVeto = aEvent.checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
  bool extraElectronVeto = aEvent.checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
  bool postSynch = postSynchTau && postSynchMuon && !diMuonVeto && !extraMuonVeto && !extraElectronVeto;
  if(!myEventProxy.pairs->size()) return true;

  setAnalysisObjects(myEventProxy);

  std::pair<bool, bool> goodDecayModes = checkTauDecayMode(myEventProxy);
  bool goodGenDecayMode = goodDecayModes.first;
  bool goodRecoDecayMode = goodDecayModes.second;

  if(goodGenDecayMode) fillGenDecayPlaneAngle(sampleName+"GenNoOfflineSel", eventWeight);

  int tauIDmask = 0;
  for(unsigned int iBit=0;iBit<aEvent.ntauIds;iBit++){
    if(aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
  }
  bool tauID = ( (int)aTau.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
  bool muonKinematics = aMuon.getP4().Pt()>24 && fabs(aMuon.getP4().Eta())<2.1;
  bool trigger = aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
		 aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
		 aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
		 aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1);

  if(sampleName!="Data") trigger = true; //MC trigger included in muon SF
  if(!muonKinematics || !tauID || !trigger) return true;

  ///Note: parts of the signal/control region selection are applied in the following code.
  bool SS = aTau.getCharge()*aMuon.getCharge() == 1;
  bool OS = aTau.getCharge()*aMuon.getCharge() == -1;
  bool muonIso = aMuon.getProperty(PropertyEnum::combreliso)<0.15;

  std::string categorySuffix = "";
  std::string systEffectName = "";
  for(unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT;
    iSystEffect<(unsigned int)sysEffects::DUMMY;++iSystEffect){

    sysEffects::sysEffectsEnum aSystEffect = static_cast<sysEffects::sysEffectsEnum>(iSystEffect);

    float muonScaleFactor = getLeptonCorrection(aMuon.getP4(aSystEffect).Eta(), aMuon.getP4(aSystEffect).Pt(), hadronicTauDecayModes::tauDecayMuon);
    float tauScaleFactor = getLeptonCorrection(aTau.getP4(aSystEffect).Eta(), aTau.getP4(aSystEffect).Pt(),
    					       static_cast<hadronicTauDecayModes>(aTau.getProperty(PropertyEnum::decayMode)));
    float weightSyst = getSystWeight(aSystEffect);
    float eventWeightWithSyst=eventWeight*weightSyst*muonScaleFactor*tauScaleFactor;

      TLorentzVector met4v(aPair.getMET(aSystEffect).X(),
			   aPair.getMET(aSystEffect).Y(), 0 ,
			   aPair.getMET(aSystEffect).Mod());
    aMET.setP4(met4v);

    for(unsigned int iCategory = HTTAnalyzer::jet0;iCategory<HTTAnalyzer::CP_Pi;++iCategory){
      HTTAnalyzer::muTauCategory categoryType = static_cast<HTTAnalyzer::muTauCategory>(iCategory);

      if(!passCategory(categoryType, aSystEffect)) continue;

      categorySuffix = std::to_string(iCategory);
      systEffectName = HTTAnalyzer::systEffectName(iSystEffect);
      if(systEffectName.find("CAT")!=std::string::npos){
        std::string categoryName = HTTAnalyzer::categoryName(iCategory);
        systEffectName.replace(systEffectName.find("CAT"),3,categoryName);
      }

      if(OS && muonIso){
	hNameSuffix = sampleName+"_OS_"+categorySuffix+systEffectName;

	fillControlHistos(hNameSuffix, eventWeightWithSyst, aSystEffect);
	if(aSystEffect==sysEffects::NOMINAL_SVFIT){
	  fillGenDecayPlaneAngle(hNameSuffix+"_Gen", eventWeightWithSyst);
	  fillDecayPlaneAngle(hNameSuffix+"_RefitPV", eventWeightWithSyst);
	  fillDecayPlaneAngle(hNameSuffix+"_AODPV", eventWeightWithSyst);
	  fillDecayPlaneAngle(hNameSuffix+"_GenPV", eventWeightWithSyst);
	  fillVertices(hNameSuffix+"_RefitPV");
	  fillVertices(hNameSuffix+"_AODPV");
	}
      }
      if(SS && muonIso){
	hNameSuffix = sampleName+"_SS_"+categorySuffix+systEffectName;
	fillControlHistos(hNameSuffix, eventWeightWithSyst);
      }
      if(OS){
	hNameSuffix = sampleName+"_OSnoMuIso_"+categorySuffix+systEffectName;
	myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeightWithSyst);
	myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(aSystEffect),eventWeightWithSyst);
  //needed for antiiso control region//
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4(sysEffects::NOMINAL_SVFIT).M(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.getP4(sysEffects::NOMINAL_SVFIT).M(),eventWeight);
      }
      if(SS){
	hNameSuffix = sampleName+"_SSnoMuIso_"+categorySuffix+systEffectName;
	myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeightWithSyst);
	myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(aSystEffect),eventWeightWithSyst);
  //needed for antiiso control region//
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4(sysEffects::NOMINAL_SVFIT).M(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.getP4(sysEffects::NOMINAL_SVFIT).M(),eventWeight);
      }
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
