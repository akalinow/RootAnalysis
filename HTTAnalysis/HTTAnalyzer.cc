#include <sstream>
#include <bitset>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName){

  //pileupCalc.py -i lumiSummary_Run2016BCDE_PromptReco_v12.json
  //--inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
  //--calcMode true --minBiasXsec 69200 --maxPileupBin 60 --numPileupBins 600 Data_Pileup_Cert_271036-277148.root
  
  ///Load ROOT file with PU histograms.
  std::string filePath = "Data_Pileup_2016_BCDEFG_v26.root";
  filePath = "Data_Pileup_2016_July21.root";
  puDataFile_ = new TFile(filePath.c_str());

  filePath = "MC_Spring16_PU25ns_V1.root";
  puMCFile_ = new TFile(filePath.c_str());

#pragma omp critical
  {
    filePath = "htt_scalefactors_v4.root";
    TFile aFile(filePath.c_str());
    scaleWorkspace = (RooWorkspace*)aFile.Get("w")->Clone("w");
    aFile.Close();
  }
  
  ntupleFile_ = 0;
  hStatsFromFile = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
  if(puDataFile_) delete puDataFile_;
  if(puMCFile_) delete puMCFile_;
  if(scaleWorkspace) delete scaleWorkspace;

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
void HTTAnalyzer::finalize(){ myHistos_->finalizeHistograms(0,1.0); }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> HTTAnalyzer::getSeparatedJets(const EventProxyHTT & myEventProxy,
						float deltaR){

  std::vector<HTTParticle> separatedJets;
  
  for(auto aJet: *myEventProxy.jets){
    float dRTau = aJet.getP4().DeltaR(aTau.getP4());
    float dRMu = aMuon.getP4().DeltaR(aTau.getP4());
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
  aJet = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
  nJets30 = count_if(aSeparatedJets.begin(), aSeparatedJets.end(),[](const HTTParticle & aJet){return aJet.getP4().Pt()>30;});
  
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
void HTTAnalyzer::fillControlHistos(const std::string & hNameSuffix, float eventWeight){

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.getNPV(),eventWeight);

  ///Fill SVfit and visible masses
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4SVFit().M(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.getP4().M(),eventWeight);
  
  ///Fill muon
  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtMuon"+hNameSuffix,aMuon.getP4().Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+hNameSuffix,aMuon.getP4().Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMuon"+hNameSuffix,aMuon.getP4().Phi(),eventWeight);
  myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
  myHistos_->fill1DHistogram("h1DnPCAMuon"+hNameSuffix,aMuon.getPCARefitPV().Mag(),eventWeight);

  ///Fill tau
  myHistos_->fill1DHistogram("h1DPtTau"+hNameSuffix,aTau.getP4().Pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+hNameSuffix,aTau.getP4().Eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+hNameSuffix,aTau.getP4().Phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DIDTau"+hNameSuffix,aTau.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw) ,eventWeight);  
  myHistos_->fill1DHistogram("h1DStatsDecayMode"+hNameSuffix, aTau.getProperty(PropertyEnum::decayMode), eventWeight);
  myHistos_->fill1DHistogram("h1DnPCATau"+hNameSuffix,aTau.getPCARefitPV().Mag(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtTauLeadingTk"+hNameSuffix,aTau.getProperty(PropertyEnum::leadChargedParticlePt),eventWeight);
  
  ///Fill jets info           
  myHistos_->fill1DHistogram("h1DStatsNJ30"+hNameSuffix,nJets30,eventWeight);
  if(nJets30>0){
    myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet.getP4().Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet.getP4().Eta(),eventWeight);
    myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet.getProperty(PropertyEnum::bCSVscore),eventWeight);
  }
  
  myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aPair.getMET().Mod(),eventWeight);

  fillDecayPlaneAngle(hNameSuffix, eventWeight);

  if(aJet.getProperty(PropertyEnum::bCSVscore)>0.8){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet.getP4().Pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet.getP4().Eta(),eventWeight);
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

  myHistos_->fill1DHistogram("h1DVxPullX_"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullY_"+sysType,pullY);
  myHistos_->fill1DHistogram("h1DVxPullZ_"+sysType,pullZ);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

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

  bool postSynchTau = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::postSynchTau);
  bool postSynchMuon = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::postSynchMuon);
  bool diMuonVeto = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::diMuonVeto);
  bool thirdLeptonVeto = myEventProxy.event->checkSelectionBit(SelectionBitsEnum::thirdLeptonVeto);
  //if(diMuonVeto || thirdLeptonVeto) return true;  
  //if(!myEventProxy.pairs->size()) return true;

  setAnalysisObjects(myEventProxy);
  float muonScaleFactor = getLeptonCorrection(aMuon.getP4().Eta(), aMuon.getP4().Pt(), hadronicTauDecayModes::tauDecayMuon);
  float tauScaleFactor = getLeptonCorrection(aTau.getP4().Eta(), aTau.getP4().Pt(),
					     static_cast<hadronicTauDecayModes>(aTau.getProperty(PropertyEnum::decayMode)));
  eventWeight*=muonScaleFactor*tauScaleFactor;


  /////TEST
  /*
    if(aTau.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero){
    std::cout<<"Tau pt: "<<aTau.getP4().Perp()<<" mass: "<<aTau.getP4().M()
    <<" leading tk: "<<aTau.getProperty(PropertyEnum::leadChargedParticlePt)
    <<std::endl;
    }
  */
  
  std::pair<bool, bool> goodDecayModes = checkTauDecayMode(myEventProxy);
  bool goodGenDecayMode = goodDecayModes.first;
  bool goodRecoDecayMode = goodDecayModes.second;

  if(goodGenDecayMode) fillGenDecayPlaneAngle(sampleName+"GenNoOfflineSel", eventWeight);
  
  ///This stands for core selection, that is common to all regions.
  bool tauKinematics = aTau.getP4().Pt()>30 && fabs(aTau.getP4().Eta())<2.3;
  int tauIDmask = 0;
  
  for(unsigned int iBit=0;iBit<aEvent.ntauIds;iBit++){
    if(aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
    if(aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
  }

  bool tauID = ( (int)aTau.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
  bool muonKinematics = aMuon.getP4().Pt()>24 && fabs(aMuon.getP4().Eta())<2.1;
  bool trigger = aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) || aMuon.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22);
  if(sampleName!="Data") trigger = true; //MC trigger included in muon SF
  bool zeroJets = (nJets30==0);
									   
  bool cpMuonSelection = aMuon.getPCARefitPV().Perp()>0.003;    
  bool cpTauSelection = (aTau.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero && aTau.getPCARefitPV().Mag()>0.003) ||
                        (aTau.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero &&
			 isOneProng(aTau.getProperty(PropertyEnum::decayMode))); 
  bool cpSelection = cpMuonSelection && cpTauSelection;
  /*
  std::cout<<" tauKinematics: "<<tauKinematics
           <<" tauID: "<<tauID
           <<" muonKinematics: "<<muonKinematics
           <<" trigger: "<<trigger
           <<" cpSelection: "<<cpSelection
           <<std::endl;
  */
  if(!tauKinematics || !tauID || !muonKinematics || !trigger) return true;
  //if(!cpSelection) return true;
  
  ///Note: parts of the signal/control region selection are applied in the following code.
  bool SS = aTau.getCharge()*aMuon.getCharge() == 1;
  bool OS = aTau.getCharge()*aMuon.getCharge() == -1;
  bool baselineSelection = OS && aPair.getMTMuon()<50 && aMuon.getProperty(PropertyEnum::combreliso)<0.15;
  bool baselineSelectionNoMT = OS && aMuon.getProperty(PropertyEnum::combreliso)<0.15;
  bool wSelection = aPair.getMTMuon()>80 && aMuon.getProperty(PropertyEnum::combreliso)<0.15;
  bool qcdSelectionSS = SS;
  bool qcdSelectionOS = OS;
  //bool ttSelection = aJet.getProperty(PropertyEnum::bDiscriminator)>0.5 && nJets30>1;
  bool ttSelection = nJets30>1;
  bool mumuSelection =  aPair.getMTMuon()<40 &&  aMuon.getProperty(PropertyEnum::combreliso)<0.15 && aPair.getP4().M()>85 && aPair.getP4().M()<95;
  /*
  bool postSynchMuon = aEvent.checkSelectionBit(SelectionBitsEnum::postSynchMuon);
  if(postSynchMuon!=(aMuon.getProperty(PropertyEnum::combreliso)<0.15))
    std::cout<<"postSynchMuon: "<<postSynchMuon
	     <<" aMuon.getProperty(PropertyEnum::combreliso)<0.15: "<<(aMuon.getProperty(PropertyEnum::combreliso)<0.15)<<std::endl;
  */
  ///Histograms for the baseline selection  
  if(baselineSelection){
    fillControlHistos(hNameSuffix, eventWeight);
    fillDecayPlaneAngle(hNameSuffix+"RefitPV", eventWeight);
    fillDecayPlaneAngle(hNameSuffix+"AODPV", eventWeight);
    fillDecayPlaneAngle(hNameSuffix+"GenPV", eventWeight);

    fillVertices(hNameSuffix+"RefitPV");
    fillVertices(hNameSuffix+"AODPV");
  }
  if(baselineSelectionNoMT){
    hNameSuffix = sampleName+"fullMt";
    myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon(),eventWeight); 
  }

  ///Histograms for the QCD control region
  if(qcdSelectionSS){
    hNameSuffix = sampleName+"qcdselSS";
    ///SS ans OS isolation histograms are filled only for mT<40 to remove possible contamination
    //from TT in high mT region.
    if(aPair.getMTMuon()<50) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
    ///Fill SS histos in signal mu isolation region. Those histograms
    ///provide shapes for QCD estimate in signal region and in various control regions.
    ///If control region has OS we still use SS QCD estimate.
    if(aMuon.getProperty(PropertyEnum::combreliso)<0.15) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"fullMt",aPair.getMTMuon(),eventWeight); 
    if(aPair.getMTMuon()<50 && aMuon.getProperty(PropertyEnum::combreliso)<0.15) fillControlHistos(hNameSuffix, eventWeight);
    if(aPair.getMTMuon()>80 && aMuon.getProperty(PropertyEnum::combreliso)<0.15){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS",aPair.getMTMuon(),eventWeight);    
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS",aPair.getMTMuon(),eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"wselOS",eventWeight);
    }
    if(ttSelection){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"ttselOS",aPair.getMTMuon(),eventWeight);
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"ttselSS",aPair.getMTMuon(),eventWeight);
      myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix+"ttselOS",aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
    }
    if(mumuSelection){
      myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"mumuselSS",aPair.getMET().Mod(),eventWeight);
      myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"mumuselOS",aPair.getMET().Mod(),eventWeight);
    }
  }
  ///Make QCD shape histograms for specific selection.
  ///Using the same SS/OS scaling factor for now.    
  if(qcdSelectionOS){
    hNameSuffix = sampleName+"qcdselOS";
    if(aPair.getMTMuon()<50) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
  }

  ///Histograms for the WJet control region. 
  if(wSelection){
    hNameSuffix = sampleName+"wsel";
    if(OS){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aPair.getMTMuon(),eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"OS",eventWeight);
    }
    if(SS) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aPair.getMTMuon(),eventWeight);
  }

  ///Histograms for the tt control region
  if(ttSelection){
    hNameSuffix = sampleName+"ttsel";
    if(OS){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aPair.getMTMuon(),eventWeight);
      myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix+"OS",aMuon.getProperty(PropertyEnum::combreliso),eventWeight);
    }
    if(SS) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aPair.getMTMuon(),eventWeight);    
  }

  ///Histograms for the DY->mu mu
  if(mumuSelection){
    hNameSuffix = sampleName+"mumusel";
    if(OS) myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"OS",aPair.getMET().Mod(),eventWeight);
    if(SS) myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"SS",aPair.getMET().Mod(),eventWeight);
  }
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
