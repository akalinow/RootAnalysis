#include <sstream>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName){

  ///Load ROOT file with PU histograms.
  std::string filePath = "Data_Pileup_2015D_Feb02.root";
  puDataFile_ = new TFile(filePath.c_str());

  filePath = "MC_Spring15_PU25_Startup.root";
  puMCFile_ = new TFile(filePath.c_str());

  ntupleFile_ = 0;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
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
void HTTAnalyzer::initialize(TFileDirectory& aDir,
			     pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  myHistos_ = new HTTHistograms(&aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::getPreselectionEff(const EventProxyHTT & myEventProxy){

    TH1F *hStatsFromFile = (TH1F*)myEventProxy.getTTree()->GetCurrentFile()->Get("m2n/hStats");

    std::string hName = "h1DStats"+getSampleName(myEventProxy);
    TH1F *hStats = myHistos_->get1DHistogram(hName.c_str(),true);

    float genWeight = getGenWeight(myEventProxy);
    
    hStats->SetBinContent(2,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1))*genWeight);   
    hStats->SetBinContent(3,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3))*genWeight);
    delete hStatsFromFile;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return "Data";
  if(myEventProxy.wevent->sample()==1){
    int decayModeBoson = myEventProxy.wevent->decayModeBoson();
    if(decayModeBoson==7) return "DYJetsMuMu";
    else if(decayModeBoson==6) return "DYJetsEE";
    else if(decayModeBoson==0) return "DYJetsMuTau";
    else return "DYJetsOther";    
  }
  if(myEventProxy.wevent->sample()==11) return "DYJetsLowM";
  if(myEventProxy.wevent->sample()==2){
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsToLNu_HT-100To200")!=std::string::npos) return "WJetsHT100to200";
    if(fileName.find("WJetsToLNu_HT-200To400")!=std::string::npos) return "WJetsHT200to400";
    if(fileName.find("WJetsToLNu_HT-400To600")!=std::string::npos) return "WJetsHT400to600";
    if(fileName.find("WJetsToLNu_HT-600ToInf")!=std::string::npos) return "WJetsHT600toInf";
    return "WJetsHT0";
  }
  if(myEventProxy.wevent->sample()==3) return "TTbar";
  if(myEventProxy.wevent->sample()==5) return "H";
  if(myEventProxy.wevent->sample()==6) return "A";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getPUWeight(const EventProxyHTT & myEventProxy){

  ///Load histogram only once,later fetch it from vector<TH1F*>
  ///At the same time divide the histogram to get the weight.
  ///First load Data PU
  if(!hPUVec_.size())  hPUVec_.resize(64);

  if(!hPUVec_[myEventProxy.wevent->sample()]){
    std::string hName = "pileup";
    TH1F *hPUData = (TH1F*)puDataFile_->Get(hName.c_str());
    TH1F *hPUSample = (TH1F*)puMCFile_->Get(hName.c_str());
    ///Normalise both histograms.
    hPUData->Scale(1.0/hPUData->Integral(0,hPUData->GetNbinsX()+1));
    hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
    ///
    hPUData->SetDirectory(0);
    hPUSample->SetDirectory(0);
    hPUData->Divide(hPUSample);
    hPUData->SetName(("h1DPUWeight"+getSampleName(myEventProxy)).c_str());
    ///To get uniform treatment put weight=1.0 for under/overlow bins of
    ///data PU, as nPU for data has a dummy value.
    if(getSampleName(myEventProxy)=="Data"){
      hPUData->SetBinContent(0,1.0);
      hPUData->SetBinContent(hPUData->GetNbinsX()+1,1.0);
    }
    hPUVec_[myEventProxy.wevent->sample()] =  hPUData;
  }

  int iBinPU = hPUVec_[myEventProxy.wevent->sample()]->FindBin(myEventProxy.wevent->npu());
  return  hPUVec_[myEventProxy.wevent->sample()]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getGenWeight(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return 1.0;
  ///generator weight broken in miniAODv2
  /*
  if(myEventProxy.wevent->sample()==1) return myEventProxy.wevent->genevtweight()/23443.423;  
  if(myEventProxy.wevent->sample()==2) return myEventProxy.wevent->genevtweight()/225892.45;  
  if(myEventProxy.wevent->sample()==3) return myEventProxy.wevent->genevtweight()/6383;
  */
  
  if(myEventProxy.wevent->sample()==2){
    //https://twiki.cern.ch/twiki/pub/CMS/HiggsToTauTauWorking2015/WplusHtWeights.xls
    //NLO to LO scaling removed, as NLO cross section used in normalisation.
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsToLNu_HT-100To200")!=std::string::npos) return 0.1352710705/1.2137837838;
    if(fileName.find("WJetsToLNu_HT-200To400")!=std::string::npos) return 0.076142149/1.2137837838;
    if(fileName.find("WJetsToLNu_HT-400To600")!=std::string::npos) return 0.0326980819/1.2137837838;
    if(fileName.find("WJetsToLNu_HT-600ToInf")!=std::string::npos) return 0.0213743732/1.2137837838;
    return 0.8520862372/1.2137837838;
  }
  
  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillControlHistos(float eventWeight, std::string & hNameSuffix){

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.npv(),eventWeight);

  ///Fill SVfit and visible masses
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.svfit(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.m_vis(),eventWeight);
  
  ///Fill muon
  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aMuon.mt(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtMuon"+hNameSuffix,aMuon.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+hNameSuffix,aMuon.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix,aMuon.iso(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMuon"+hNameSuffix,aMuon.phi(),eventWeight);

  ///Fill tau
  myHistos_->fill1DHistogram("h1DPtTau"+hNameSuffix,aTau.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+hNameSuffix,aTau.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+hNameSuffix,aTau.phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DIDTau"+hNameSuffix,aTau.tauID(byCombinedIsolationDeltaBetaCorrRaw3Hits) ,eventWeight);  
  myHistos_->fill1DHistogram("h1DStatsDecayMode"+hNameSuffix, aTau.decayMode(), eventWeight);

  ///Fill leading tau track pt
  myHistos_->fill1DHistogram("h1DPtTauLeadingTk"+hNameSuffix,aTau.leadingTk().Pt(),eventWeight);
  ///Fill jets info           
  myHistos_->fill1DHistogram("h1DStatsNJets30"+hNameSuffix,nJets30,eventWeight);
  myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet.csvtag(),eventWeight);

  myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aMET.metpt(),eventWeight);
  
  if(aJet.bjet()){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet.pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet.eta(),eventWeight);
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////				      
void HTTAnalyzer::fillDecayPlaneAngle(float eventWeight, std::string & hNameSuffix){

  ///Method from http://arxiv.org/abs/1108.0670 (Berger)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles;

  TLorentzVector positiveLeadingTk, negativeLeadingTk;
  TLorentzVector positive_nPCA, negative_nPCA;

  if(aMuon.charge()>0){
    positiveLeadingTk = aMuon.leadingTk();
    if(hNameSuffix.find("AODPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCAGenvx(),0);

    negativeLeadingTk = aTau.leadingTk();
    if(hNameSuffix.find("AODPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCAGenvx(),0);
  }
  else{    
    positiveLeadingTk = aTau.leadingTk();
    if(hNameSuffix.find("AODPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCAGenvx(),0);
    
    negativeLeadingTk = aMuon.leadingTk();
    if(hNameSuffix.find("AODPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCAGenvx(),0);
  }

  angles = angleBetweenPlanes(negativeLeadingTk,negative_nPCA,
			      positiveLeadingTk,positive_nPCA);


  if(aEvent.nTracksInRefit()>3 && positive_nPCA.Vect().Mag()>0.01 && negative_nPCA.Vect().Mag()>0.01)
    myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);

  if(aEvent.nTracksInRefit()>3 && positive_nPCA.Vect().Mag()>0.002){
  float cosPositive =  positive_nPCA.Vect().Unit()*aGenPositiveTau.nPCA().Unit();
  myHistos_->fill1DHistogram("h1DCosPhi_CosPositive"+hNameSuffix,cosPositive,eventWeight);
  
  float cosNegative = negative_nPCA.Vect().Unit()*aGenNegativeTau.nPCA().Unit();
  myHistos_->fill1DHistogram("h1DCosPhi_CosNegative"+hNameSuffix,cosNegative,eventWeight);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillGenDecayPlaneAngle(float eventWeight, std::string & hNameSuffix){

  ///Method from http://arxiv.org/abs/1108.0670 (Berger)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame  
  std::pair<float,float>  angles;

  TLorentzVector positiveLeadingTk, negativeLeadingTk;
  TLorentzVector positive_nPCA, negative_nPCA;

  positiveLeadingTk = aGenPositiveTau.leadingTk();
  positive_nPCA = TLorentzVector(aGenPositiveTau.nPCA(),0);

  negativeLeadingTk = aGenNegativeTau.leadingTk();
  negative_nPCA = TLorentzVector(aGenNegativeTau.nPCA(),0);

  angles = angleBetweenPlanes(negativeLeadingTk,negative_nPCA,
			      positiveLeadingTk,positive_nPCA);
  
  myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<float,float> HTTAnalyzer::angleBetweenPlanes(const TLorentzVector &tau1, 
						       const TLorentzVector &tau1Daughter,
						       const TLorentzVector &tau2, 
						       const TLorentzVector &tau2Daughter){
  //Boost all 4v to (tau1+tau2) rest frame
  TVector3 boost = (tau1+tau2).BoostVector();

  TLorentzVector tau1Star = tau1;
  TLorentzVector tau2Star = tau2;
  tau1Star.Boost(-boost);
  tau2Star.Boost(-boost);
  
  TLorentzVector tau1DaughterStar = tau1Daughter;
  tau1DaughterStar.Boost(-boost);
  
  TLorentzVector tau2DaughterStar = tau2Daughter;
  tau2DaughterStar.Boost(-boost);

  //define common direction and normal vectors to decay planes
  TVector3 direction = tau1Star.Vect().Unit();
  TVector3 n1 = ( direction.Cross( tau1DaughterStar.Vect() ) ).Unit(); 
  TVector3 n2 = ( direction.Cross( tau2DaughterStar.Vect() ) ).Unit(); 

  ///angle between decay planes
  float phi=TMath::ACos(n1*n2);

  ///angle between tau1 and tau2 daughter momentums
  float rho=TMath::ACos( (tau1DaughterStar.Vect().Unit() )*(tau2DaughterStar.Vect().Unit() ) );

  return std::make_pair(phi,rho);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<std::string> HTTAnalyzer::getTauDecayName(int decModeMinus, int decModePlus){

  std::vector<std::string> types;

  if(decModeMinus==tauDecay1ChargedPion0PiZero && decModePlus==tauDecay1ChargedPion0PiZero) types.push_back("PiPi0Pi0");

  if(isOneProng(decModeMinus) && isOneProng(decModePlus) ) types.push_back("1Prong1Prong");

  if( (decModeMinus==tauDecay1ChargedPion0PiZero && isLepton(decModePlus) ) ||
      (isLepton(decModeMinus) && decModePlus==tauDecay1ChargedPion0PiZero)) types.push_back("Lepton1Prong0Pi0");
    
  if( (isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
      ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ) types.push_back("Lepton1Prong");

  if(decModeMinus==tauDecay1ChargedPion1PiZero && decModePlus==tauDecay1ChargedPion1PiZero ) types.push_back("PiPlusPiMinus2Pi0");


  if( isOneProng(decModeMinus) && decModeMinus!=tauDecay1ChargedPion0PiZero && 
      isOneProng(decModePlus) && decModePlus!=tauDecay1ChargedPion0PiZero )   types.push_back("1Prong1ProngXPi0");

  if(isLepton(decModePlus) && isLepton(decModeMinus)) types.push_back("LeptonLepton");

  return types;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::isOneProng(int decMode){
  if(decMode==tauDecay1ChargedPion0PiZero ||
     decMode==tauDecay1ChargedPion1PiZero ||
     decMode==tauDecay1ChargedPion2PiZero ||
     decMode==tauDecay1ChargedPion3PiZero ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::isLepton(int decMode){
  if(decMode==tauDecaysElectron || decMode==tauDecayMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<Wjet> HTTAnalyzer::getSeparatedJets(const EventProxyHTT & myEventProxy, const Wtau & aTau,
						const Wmu & aMuon, float deltaR){

  std::vector<Wjet> separatedJets;
  
  for(auto aJet: *myEventProxy.wjet){
    float dRTau = sqrt(pow(aJet.eta() - aTau.eta(),2) + pow(aJet.phi() - aTau.phi(),2));
    float dRMu = sqrt(pow(aJet.eta() - aMuon.eta(),2) + pow(aMuon.phi() - aTau.phi(),2));
    if(dRTau>deltaR && dRMu>deltaR) separatedJets.push_back(aJet);
  }

  return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  aEvent = *myEventProxy.wevent;  
  aPair = (*myEventProxy.wpair)[0];
  aTau = (*myEventProxy.wtau)[0];
  if(myEventProxy.wtauGen && myEventProxy.wtauGen->size()){
    aGenNegativeTau = (*myEventProxy.wtauGen)[0];
    aGenPositiveTau = (*myEventProxy.wtauGen)[1];
  }
  aMuon = (*myEventProxy.wmu)[0];
  aMET = (*myEventProxy.wmet)[0];
  aSeparatedJets = getSeparatedJets(myEventProxy, aTau, aMuon, 0.5);
  aJet = aSeparatedJets.size() ? aSeparatedJets[0] : Wjet();
  nJets30 = count_if(aSeparatedJets.begin(), aSeparatedJets.end(),[](const Wjet & aJet){return aJet.pt()>30;});

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<bool, bool> HTTAnalyzer::checkTauDecayMode(const EventProxyHTT & myEventProxy){

  bool goodGenDecayMode = false;
  bool goodRecoDecayMode = false;
  std::vector<std::string> decayNamesGen = getTauDecayName(myEventProxy.wevent->decModeMinus(), myEventProxy.wevent->decModePlus());
  std::vector<std::string> decayNamesReco = getTauDecayName(aTau.decayMode(), tauDecayMuon);
  for(auto it: decayNamesGen) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodGenDecayMode = true;
  for(auto it: decayNamesReco) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodRecoDecayMode = true;

  return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::addBranch(TTree *tree){

  //tree->Branch("muonPt",&muonPt);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  std::string sampleName = getSampleName(myEventProxy);
  std::string hNameSuffix = sampleName;
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);
  float eventWeight = puWeight*genWeight;

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,0,eventWeight);
  getPreselectionEff(myEventProxy);
  /////////////////////////////////////////////////////////////////////////////
  if(!myEventProxy.wpair->size() || !myEventProxy.wtau->size() || !myEventProxy.wmu->size()) return true;

  setAnalysisObjects(myEventProxy);

  std::pair<bool, bool> goodDecayModes = checkTauDecayMode(myEventProxy);
  bool goodGenDecayMode = goodDecayModes.first;
  bool goodRecoDecayMode = goodDecayModes.second;

  if(goodGenDecayMode && goodRecoDecayMode){
    std::string hNameSuffixCP1 = sampleName+"Gen";
    hNameSuffixCP1 = sampleName+"GenNoOfflineSel";
    fillGenDecayPlaneAngle(eventWeight, hNameSuffixCP1);
  }

  ///This stands for core selection, that is common to all regions.
  bool tauKinematics = aTau.pt()>30 && fabs(aTau.eta())<2.3;
  bool tauID = aTau.tauID(byMediumCombinedIsolationDeltaBetaCorr3Hits);
  bool muonKinematics = aMuon.pt()>19 && fabs(aMuon.eta())<2.1;
  bool trigger = aPair.trigger(HLT_IsoMu17_eta2p1);
  if(sampleName=="Data") trigger = aPair.trigger(HLT_IsoMu18);
  bool extraRequirements = aTau.decayMode()!=5 && aTau.decayMode()!=6 && nJets30==0;

  if(!myEventProxy.wpair->size()) return true;
  if(!tauKinematics || !muonKinematics || !trigger) return true;
  //if(!tauKinematics || !tauID || !muonKinematics || !trigger) return true;
  //if(!extraRequirements) return true;

  ///Note: parts of the signal/control region selection are applied in the following code.
  ///FIXME AK: this should be made in a more clear way.
  bool baselineSelection = aPair.diq()==-1 && aMuon.mt()<40 && aMuon.iso()<0.1;
  bool wSelection = aMuon.mt()>60 && aMuon.iso()<0.1;
  bool qcdSelectionSS = aPair.diq()==1;
  bool qcdSelectionOS = aPair.diq()==-1;
  bool ttSelection = aJet.csvtag()>0.9 && nJets30>1;
  bool mumuSelection =  aMuon.mt()<40 && aMuon.iso()<0.1 && aPair.m_vis()>85 && aPair.m_vis()<95;

  ///Fill variables stored in TTree
  muonPt = aMuon.pt();

  ///Histograms for the baseline selection  
  if(baselineSelection){
    fillControlHistos(eventWeight, hNameSuffix);
    if(goodGenDecayMode){
      std::string hNameSuffixCP = hNameSuffix+"RefitPV";    
      fillDecayPlaneAngle(eventWeight, hNameSuffixCP);
      hNameSuffixCP = hNameSuffix+"AODPV";
      fillDecayPlaneAngle(eventWeight, hNameSuffixCP);
      hNameSuffixCP = hNameSuffix+"GenPV";
      fillDecayPlaneAngle(eventWeight, hNameSuffixCP);
      hNameSuffixCP = hNameSuffix+"Gen";
      fillDecayPlaneAngle(eventWeight, hNameSuffixCP);
    }
  }

  ///Histograms for the QCD control region
  if(qcdSelectionSS){
    hNameSuffix = sampleName+"qcdselSS";
    ///SS ans OS isolation histograms are filled only for mT<40 to remove possible contamnation
    //from TT in high mT region.
    if(aMuon.mt()<40) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
    ///Fill SS histos in signal mu isolation region. Those histograms
    ///provide shapes for QCD estimate in signal region and in various control regions.
    ///If control region has OS we still use SS QCD estimate.
    if(aMuon.mt()<40 && aMuon.iso()<0.1) fillControlHistos(eventWeight, hNameSuffix);
    if(aMuon.mt()>60){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS",aMuon.mt(),eventWeight);    
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS",aMuon.mt(),eventWeight);    
    }
    if(ttSelection){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"ttselOS",aMuon.mt(),eventWeight);
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"ttselSS",aMuon.mt(),eventWeight);
      myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix+"ttselOS",aMuon.iso(),eventWeight);
    }
    if(mumuSelection){
      myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"mumuselSS",aMET.metpt(),eventWeight);
      myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"mumuselOS",aMET.metpt(),eventWeight);
    }
  }
  ///Make QCD shape histograms for specific selection.
  ///Using the same SS/OS scaling factor for now.    
  if(qcdSelectionOS){
    hNameSuffix = sampleName+"qcdselOS";
    if(aMuon.mt()<40) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
  }

  ///Histograms for the WJet control region. 
  if(wSelection){
    hNameSuffix = sampleName+"wsel";
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aMuon.mt(),eventWeight);
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aMuon.mt(),eventWeight);
  }

  ///Histograms for the tt control region
  if(ttSelection){
    hNameSuffix = sampleName+"ttsel";
    if(aPair.diq()==-1){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aMuon.mt(),eventWeight);
      myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix+"OS",aMuon.iso(),eventWeight);
    }
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aMuon.mt(),eventWeight);    
  }

  ///Histograms for the DY->mu mu
  if(mumuSelection){
    hNameSuffix = sampleName+"mumusel";
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"OS",aMET.metpt(),eventWeight);
    if(aPair.diq()==1) myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix+"SS",aMET.metpt(),eventWeight);
  }
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
