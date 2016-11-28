#include <sstream>
#include <bitset>

#include "HTTAnalyzerTT.h"
#include "HTTHistogramsTT.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzerTT::HTTAnalyzerTT(const std::string & aName):Analyzer(aName){

  //pileupCalc.py -i lumiSummary_Run2016BCDE_PromptReco_v12.json
  //--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
  //--calcMode true --minBiasXsec 69200 --maxPileupBin 60 --numPileupBins 600 Data_Pileup_Cert_271036-277148.root
  
  ///Load ROOT file with PU histograms.
  std::string filePath = "Data_Pileup_2016_BCDEFG_v26.root";
  filePath = "Data_Pileup_2016_July22.root";
  puDataFile_ = new TFile(filePath.c_str());

  filePath = "MC_Spring16_PU25ns_V1.root";
  puMCFile_ = new TFile(filePath.c_str());

#pragma omp critical
  {
    filePath = "htt_scalefactors_v5.root";
    TFile aFile(filePath.c_str());
    scaleWorkspace = (RooWorkspace*)aFile.Get("w")->Clone("w");
    aFile.Close();
  }
  
  ntupleFile_ = 0;
  hStatsFromFile = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzerTT::~HTTAnalyzerTT(){

  if(myHistos_) delete myHistos_;
  if(puDataFile_) delete puDataFile_;
  if(puMCFile_) delete puMCFile_;
  if(scaleWorkspace) delete scaleWorkspace;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTAnalyzerTT::clone() const{

  HTTAnalyzerTT* clone = new HTTAnalyzerTT(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::initialize(TDirectory* aDir,
			     pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  myHistos_ = new HTTHistogramsTT(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::finalize(){ myHistos_->finalizeHistograms(0,1.0); }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> HTTAnalyzerTT::getSeparatedJets(const EventProxyHTT & myEventProxy,
						float deltaR){

  std::vector<HTTParticle> separatedJets;
  
  for(auto aJet: *myEventProxy.jets){
    float dRTau1 = aJet.getP4().DeltaR(aTau1.getP4());
    float dRTau2 = aJet.getP4().DeltaR(aTau2.getP4());
    bool loosePFJetID = aJet.getProperty(PropertyEnum::PFjetID)>=1;
    if(dRTau1>deltaR && dRTau2>deltaR && loosePFJetID) separatedJets.push_back(aJet);
  }

  return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  aEvent = *myEventProxy.event;  
  aPair = (*myEventProxy.pairs)[0];
  aTau1 = aPair.getLeg1();
  aTau2 = aPair.getLeg2();

  TLorentzVector met4v(aPair.getMET().X(),
		       aPair.getMET().Y(),
		       0,
		       aPair.getMET().Mod());
  
  aMET = HTTParticle();
  aMET.setP4(met4v);
  
  aGenTau1 = HTTParticle();
  aGenTau2 = HTTParticle();

  if(myEventProxy.genLeptons && myEventProxy.genLeptons->size()){
    HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);    
    if(aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecaysElectron) {
      if(aGenTau.getP4().DeltaR(aTau1.getP4())<0.5)
	aGenTau1 = aGenTau;
      else if(aGenTau.getP4().DeltaR(aTau2.getP4())<0.5)
	aGenTau2 = aGenTau;
    }
    if(myEventProxy.genLeptons->size()>1){
      aGenTau =  myEventProxy.genLeptons->at(1);
      if(aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=tauDecaysElectron) {
	if(aGenTau.getP4().DeltaR(aTau1.getP4())<0.5)
	  aGenTau1 = aGenTau;
	else if(aGenTau.getP4().DeltaR(aTau2.getP4())<0.5)
	  aGenTau2 = aGenTau;
    }
    }
  }
  
  aSeparatedJets = getSeparatedJets(myEventProxy, 0.5);
  aJet1 = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
  aJet2 = aSeparatedJets.size()>1 ? aSeparatedJets[1] : HTTParticle();
  nJets30 = count_if(aSeparatedJets.begin(), aSeparatedJets.end(),[](const HTTParticle & aJet){return aJet.getP4().Pt()>30;});
  nJetsInGap30 = 0;
  if(nJets30>=2){
    for(unsigned int iJet=2; iJet<aSeparatedJets.size(); ++iJet){
      if( (aSeparatedJets.at(iJet).getP4().Eta()>aJet1.getP4().Eta()&&aSeparatedJets.at(iJet).getP4().Eta()<aJet2.getP4().Eta()) ||
          (aSeparatedJets.at(iJet).getP4().Eta()<aJet1.getP4().Eta()&&aSeparatedJets.at(iJet).getP4().Eta()>aJet2.getP4().Eta()) ){
        if(aSeparatedJets.at(iJet).getP4().Pt()>30) nJetsInGap30++;
      }                           
    }
  }

  categoryDecisions.clear();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<bool, bool> HTTAnalyzerTT::checkTauDecayMode(const EventProxyHTT & myEventProxy){

  bool goodGenDecayMode = false;
  bool goodRecoDecayMode = false;

  std::vector<std::string> decayNamesGen = getTauDecayName(aGenTau1.getProperty(PropertyEnum::decayMode),
							   aGenTau2.getProperty(PropertyEnum::decayMode));
  std::vector<std::string> decayNamesReco = getTauDecayName(aTau1.getProperty(PropertyEnum::decayMode),aTau2.getProperty(PropertyEnum::decayMode));

  for(auto it: decayNamesGen) if(it.find("1Prong1Prong")!=std::string::npos) goodGenDecayMode = true;
  for(auto it: decayNamesReco) if(it.find("1Prong1Prong")!=std::string::npos) goodRecoDecayMode = true;

  return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::addBranch(TTree *tree){/*tree->Branch("muonPt",&muonPt);*/}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::fillControlHistos(const std::string & hNameSuffix, float eventWeight){

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.getNPV(),eventWeight);

  ///Fill SVfit and visible masses
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4SVFit().M(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.getP4().M(),eventWeight);
  
  ///Fill leading tau
  myHistos_->fill1DHistogram("h1DPtLeadingTau"+hNameSuffix,aTau1.getP4().Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaLeadingTau"+hNameSuffix,aTau1.getP4().Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiLeadingTau"+hNameSuffix,aTau1.getP4().Phi(),eventWeight);
  myHistos_->fill1DHistogram("h1DIDLeadingTau"+hNameSuffix,aTau1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw) ,eventWeight);  
  myHistos_->fill1DHistogram("h1DStatsDecayModeLeadingTau"+hNameSuffix,aTau1.getProperty(PropertyEnum::decayMode), eventWeight);
  myHistos_->fill1DHistogram("h1DnPCALeadingTau"+hNameSuffix,aTau1.getPCARefitPV().Mag(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtLeadingTauLeadingTk"+hNameSuffix,aTau1.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);

  ///Fill trailing tau
  myHistos_->fill1DHistogram("h1DPtTrailingTau"+hNameSuffix,aTau2.getP4().Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTrailingTau"+hNameSuffix,aTau2.getP4().Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTrailingTau"+hNameSuffix,aTau2.getP4().Phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DIDTrailingTau"+hNameSuffix,aTau2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw) ,eventWeight);  
  myHistos_->fill1DHistogram("h1DStatsDecayMode"+hNameSuffix,aTau2.getProperty(PropertyEnum::decayMode), eventWeight);
  myHistos_->fill1DHistogram("h1DnPCATrailingTau"+hNameSuffix,aTau2.getPCARefitPV().Mag(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtTrailingTauLeadingTk"+hNameSuffix,aTau2.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);
  float higgsPt =  (aTau1.getP4() + aTau2.getP4() + aMET.getP4()).Pt();
  myHistos_->fill1DHistogram("h1DPtTauTauMET"+hNameSuffix,higgsPt,eventWeight);
  
  ///Fill jets info           
  myHistos_->fill1DHistogram("h1DStatsNJ30"+hNameSuffix,nJets30,eventWeight);
  if(nJets30>0){
    myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet1.getP4().Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet1.getP4().Eta(),eventWeight);
    myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet1.getProperty(PropertyEnum::bCSVscore),eventWeight);
  }
  if(nJets30>1){  
    myHistos_->fill1DHistogram("h1DStatsNJGap30"+hNameSuffix,nJetsInGap30,eventWeight);
    float jetsMass = (aJet1.getP4() + aJet2.getP4()).M();
    float jetsEta = std::abs(aJet1.getP4().Eta() - aJet2.getP4().Eta());
    myHistos_->fill1DHistogram("h1DWideMass2J"+hNameSuffix,jetsMass,eventWeight);
    myHistos_->fill1DHistogram("h1DEta2J"+hNameSuffix,jetsMass,eventWeight);
    myHistos_->fill1DHistogram("h1DPtTrailingJet"+hNameSuffix,aJet2.getP4().Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaTrailingJet"+hNameSuffix,aJet2.getP4().Eta(),eventWeight);
    myHistos_->fill1DHistogram("h1DCSVBtagTrailingJet"+hNameSuffix,aJet2.getProperty(PropertyEnum::bCSVscore),eventWeight);
  }
  
  myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aMET.getP4().Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMET"+hNameSuffix,aMET.getP4().Phi(),eventWeight);

  fillDecayPlaneAngle(hNameSuffix, eventWeight);

  if(aJet1.getProperty(PropertyEnum::bCSVscore)>0.8){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet1.getP4().Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet1.getP4().Eta(),eventWeight);
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzerTT::fillVertices(const std::string & sysType){
  
  TVector3 aVertexGen = aEvent.getGenPV();
  TVector3 aVertex = aEvent.getRefittedPV();
  if(sysType.find("AODPV")!=std::string::npos) aVertex = aEvent.getAODPV();
  if(sysType.find("RefitPV")!=std::string::npos) aVertex = aEvent.getRefittedPV();
  
  float pullX = aVertexGen.X() - aVertex.X();
  float pullY = aVertexGen.Y() - aVertex.Y();
  float pullZ = aVertexGen.Z() - aVertex.Z();

  myHistos_->fill1DHistogram("h1DVxPullX_"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullY_"+sysType,pullY);
  myHistos_->fill1DHistogram("h1DVxPullZ_"+sysType,pullZ);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzerTT::passCategory(const HTTAnalyzerTT::tauTauCategory & aCategory){

  if(categoryDecisions.size()) return  categoryDecisions[(int)aCategory];
  else categoryDecisions = std::vector<bool>((int)HTTAnalyzerTT::DUMMY);
  
  float jetsMass = (aJet1.getP4() + aJet2.getP4()).M();
  float jetsEta = std::abs(aJet1.getP4().Eta() - aJet2.getP4().Eta());
  float higgsPt = (aTau1.getP4() + aTau2.getP4() + aMET.getP4()).Pt();
  
  bool jet0 = nJets30==0;

  bool jet1 = (nJets30==1 || (nJets30>=2 && !(jetsMass>300 && jetsEta>2.5 && nJetsInGap30<1)));
  bool jet1_low = jet1 && (higgsPt>100 && higgsPt<170);
  
  bool jet1_high = jet1 && higgsPt>170;

  bool vbf_1d = nJets30>=2 && jetsEta>2.5 && nJetsInGap30<1;
  bool vbf_low =  vbf_1d &&
    ((higgsPt<100 && jetsMass>300) || (higgsPt>100 && jetsMass>300 && jetsMass<500));

  bool vbf_high = vbf_1d &&
    (higgsPt>100 && jetsMass>500);

  bool booted = nJets30==1 || ( nJets30>=2 && !(jetsEta>2.5 && nJetsInGap30<1 && higgsPt>100) );
  
  bool vbf_2d = nJets30>=2 && (jetsEta>2.5 && nJetsInGap30<1 && higgsPt>100);

  //////////
  // categories by tau decay modes for CP
  bool isPi1 = (aTau1.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero && aTau1.getPCARefitPV().Mag()>0.003);
  bool isPi2 = (aTau2.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero && aTau2.getPCARefitPV().Mag()>0.003);
  bool isRho1 = (aTau1.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero && isOneProng(aTau1.getProperty(PropertyEnum::decayMode)) );
  bool isRho2 = (aTau2.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero && isOneProng(aTau2.getProperty(PropertyEnum::decayMode)) );

  bool piPi = isPi1 && isPi2;
  bool piRho = (isPi1 && isRho2) || (isPi2 && isRho1);
  bool rhoRho = isRho1 && isRho2;
  
  //////////

  categoryDecisions[(int)HTTAnalyzerTT::jet0] = jet0;
  
  categoryDecisions[(int)HTTAnalyzerTT::jet1_low] = jet1_low;
  categoryDecisions[(int)HTTAnalyzerTT::jet1_high] = jet1_high;

  categoryDecisions[(int)HTTAnalyzerTT::vbf_low] = vbf_low;
  categoryDecisions[(int)HTTAnalyzerTT::vbf_high] = vbf_high;

  categoryDecisions[(int)HTTAnalyzerTT::boosted] = boosted;
  categoryDecisions[(int)HTTAnalyzerTT::vbf] = vbf_2d;

  categoryDecisions[(int)HTTAnalyzerTT::pipi] = piPi;
  categoryDecisions[(int)HTTAnalyzerTT::pirho] = piRho;
  categoryDecisions[(int)HTTAnalyzerTT::rhorho] = rhoRho;

  return categoryDecisions[(int)aCategory];
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzerTT::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  sampleName = getSampleName(myEventProxy);

  std::string hNameSuffix = sampleName;
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);
  float eventWeight = puWeight*genWeight;

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,0,eventWeight);
  myHistos_->fill1DHistogram("h1DNPartons"+hNameSuffix,myEventProxy.event->getLHEnOutPartons(),eventWeight);
  getPreselectionEff(myEventProxy);

  bool postSynchTau1 = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::postSynchMuon);
  bool postSynchTau2 = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::postSynchTau);
  bool thirdLeptonVeto = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::thirdLeptonVeto);
  if(!myEventProxy.pairs->size()) return true;

  setAnalysisObjects(myEventProxy);
  float tau1ScaleFactor = getLeptonCorrection(aTau1.getP4().Eta(), aTau1.getP4().Pt(),
					      static_cast<hadronicTauDecayModes>(aTau1.getProperty(PropertyEnum::decayMode)));
  float tau2ScaleFactor = getLeptonCorrection(aTau2.getP4().Eta(), aTau2.getP4().Pt(),
					      static_cast<hadronicTauDecayModes>(aTau2.getProperty(PropertyEnum::decayMode)));
  eventWeight*=tau1ScaleFactor*tau2ScaleFactor;


  std::pair<bool, bool> goodDecayModes = checkTauDecayMode(myEventProxy);
  bool goodGenDecayMode = goodDecayModes.first;
  bool goodRecoDecayMode = goodDecayModes.second;

  if(goodGenDecayMode) fillGenDecayPlaneAngle(sampleName+"GenNoOfflineSel", eventWeight);
  
  ///This stands for core selection, that is common to all regions.
  bool tau1Kinematics = aTau1.getP4().Pt()>50 && std::abs(aTau1.getP4().Eta())<2.1;
  bool tau2Kinematics = aTau2.getP4().Pt()>40 && std::abs(aTau2.getP4().Eta())<2.1;

  int tauIDmask=0, tauIsoTmask=0, tauIsoMmask=0, tauIsoLmask=0;
  
  for(unsigned int iBit=0;iBit<aEvent.ntauIds;iBit++){
    if(aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIsoTmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="byMediumIsolationMVArun2v1DBoldDMwLT") tauIsoMmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="byLooseIsolationMVArun2v1DBoldDMwLT") tauIsoLmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstMuonLoose3") tauIDmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
  }

  bool tau1ID = ( (int)aTau1.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
  bool tau1IsoT = ( (int)aTau1.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
  bool tau1IsoM = ( (int)aTau1.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
  bool tau1IsoL = ( (int)aTau1.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;
  bool tau2ID = ( (int)aTau2.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
  bool tau2IsoT = ( (int)aTau2.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
  bool tau2IsoM = ( (int)aTau2.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
  bool tau2IsoL = ( (int)aTau2.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;

  bool fullIso = tau1IsoT && tau2IsoT;
  bool relaxedIso = (tau1IsoM && tau2IsoL) || (tau2IsoM && tau1IsoL);
  bool antiIso = (tau1IsoM && tau2IsoL && !tau2IsoT) || (tau2IsoM && tau1IsoL && !tau1IsoT);
  

  bool trigger = aTau1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) &&
    aTau2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg); //FIXME: does not work, wrong numbering??

  //if(sampleName!="Data") trigger = true; //MC trigger included in muon SF
  if(true) trigger = true; //FIXME: problem with trigger matching??
									   
  bool cpTau1Selection = (aTau1.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero && aTau1.getPCARefitPV().Mag()>0.003) ||
                        (aTau1.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero &&
			 isOneProng(aTau1.getProperty(PropertyEnum::decayMode))); 
  bool cpTau2Selection = (aTau2.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero && aTau2.getPCARefitPV().Mag()>0.003) ||
                        (aTau2.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero &&
			 isOneProng(aTau2.getProperty(PropertyEnum::decayMode))); 
  bool cpSelection = cpTau1Selection && cpTau2Selection;
  
  /*MB TEST
  std::cout<<"Selection "
    //<<" ,"<<tau1Kinematics
    //<<" ,"<<tau1ID 
    //<<" ,"<<tau2Kinematics
    //<<" ,"<<tau2ID
    //<<" ,"<<relaxedIso
    	   <<" trigger: "<<trigger
	   <<" loose: "<< (tau1Kinematics&&tau1ID&&tau2Kinematics&&tau2ID&&relaxedIso&&trigger)
	   <<" full: "<< (tau1Kinematics&&tau1ID&&tau2Kinematics&&tau2ID&&fullIso&&trigger)
	   <<" anti: "<< (tau1Kinematics&&tau1ID&&tau2Kinematics&&tau2ID&&antiIso&&trigger)
	   <<std::endl;
  */
  if(!tau1Kinematics || !tau1ID || !tau2Kinematics || !tau2ID || !relaxedIso || !trigger) return true;
  //if(!cpSelection) return true;

 
  ///Note: parts of the signal/control region selection are applied in the following code.
  bool SS = aTau2.getCharge()*aTau1.getCharge() == 1;
  bool OS = aTau2.getCharge()*aTau1.getCharge() == -1;
  
  std::string categorySuffix = "";
  for(unsigned int iCategory = HTTAnalyzerTT::jet0;
      iCategory<HTTAnalyzerTT::DUMMY;++iCategory){

    HTTAnalyzerTT::tauTauCategory categoryType = static_cast<HTTAnalyzerTT::tauTauCategory>(iCategory);
    
    if(!passCategory(categoryType)) continue;    
    categorySuffix = std::to_string(iCategory);

    if(OS && fullIso){      
      hNameSuffix = sampleName+"_OS_"+categorySuffix;
      fillControlHistos(hNameSuffix, eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"RefitPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"AODPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"GenPV", eventWeight);
      fillVertices(hNameSuffix+"RefitPV");
      fillVertices(hNameSuffix+"AODPV");
    }
    if(OS && antiIso){
      hNameSuffix = sampleName+"_OSantiIso_"+categorySuffix;
      fillControlHistos(hNameSuffix, eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"RefitPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"AODPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"GenPV", eventWeight);
      fillVertices(hNameSuffix+"RefitPV");
      fillVertices(hNameSuffix+"AODPV");
    }
    if(SS && fullIso){
      hNameSuffix = sampleName+"_SS_"+categorySuffix;
      fillControlHistos(hNameSuffix, eventWeight);
    }      
    if(SS && antiIso){
      hNameSuffix = sampleName+"_SSantiIso_"+categorySuffix;
      fillControlHistos(hNameSuffix, eventWeight);
    }      
    if(OS && relaxedIso){
      hNameSuffix = sampleName+"_OSrelaxedIso_"+categorySuffix;
      if(!tau1IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
      if(!tau2IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
      if(tau1IsoM&&tau2IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
    }
    if(SS && relaxedIso){
      hNameSuffix = sampleName+"_SSrelaxedIso_"+categorySuffix;
      if(!tau1IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
      if(!tau2IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
      if(tau1IsoM&&tau2IsoM){
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau1.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
	myHistos_->fill1DHistogram("h1DID"+hNameSuffix,aTau2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),eventWeight);
      }
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
