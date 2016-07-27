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
void HTTAnalyzer::initialize(TDirectory* aDir,
			     pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
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
    //TESK AK hStats->SetBinContent(3,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3))*genWeight);
    hStats->SetBinContent(3,hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(2))*genWeight);///buggy weights in DY in NTUPLES_20_06_2016
    delete hStatsFromFile;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

  ///Missing sample type assignement for VBF sample in v43 ntuples.
  if(myEventProxy.wevent->sample()==0){
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("VBFHToTauTau")!=std::string::npos) return "H";
  }
  
  if(myEventProxy.wevent->sample()==0) return "Data";
  if(myEventProxy.wevent->sample()==1){
    int decayModeBoson = myEventProxy.wevent->decayModeBoson();
    std::string decayName;
    if(decayModeBoson==7) decayName = "DYJetsMuMu";
    else if(decayModeBoson==6) decayName = "DYJetsEE";
    else if(decayModeBoson==0) decayName = "DYJetsMuTau";
    else decayName = "DYJetsOther";

    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    std::string HTName = "";
    /*
    if(fileName.find("DYJetsToLL_M-50_HT-100to200")!=std::string::npos) HTName ="HT100to200";
    else if(fileName.find("DYJetsToLL_M-50_HT-200to400")!=std::string::npos) HTName = "HT200to400";
    else if(fileName.find("DYJetsToLL_M-50_HT-400to600")!=std::string::npos) HTName = "HT400to600";
    else if(fileName.find("DYJetsToLL_M-50_HT-600toInf")!=std::string::npos) HTName = "HT600toInf";
    else HTName = "HT0";
    */
    return decayName+HTName;
    
  }
  if(myEventProxy.wevent->sample()==11) return "DYJetsLowM";
  if(myEventProxy.wevent->sample()==2){
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsToLNu_HT-100To200")!=std::string::npos) return "WJetsHT100to200";
    else if(fileName.find("WJetsToLNu_HT-200To400")!=std::string::npos) return "WJetsHT200to400";
    else if(fileName.find("WJetsToLNu_HT-400To600")!=std::string::npos) return "WJetsHT400to600";
    else if(fileName.find("WJetsToLNu_HT-600ToInf")!=std::string::npos) return "WJetsHT600toInf";
    else return "WJetsHT0";
  }
  if(myEventProxy.wevent->sample()==3) return "TTbar";
  if(myEventProxy.wevent->sample()==5) return "H";
  if(myEventProxy.wevent->sample()==6) return "A";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getPUWeight(const EventProxyHTT & myEventProxy){

  ///RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1 MINIAOD
  ///is produced with PU profile matching the Run2015 data. No nee to PU rescale.
  return 1.0;

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
  //if(myEventProxy.wevent->sample()==1) return myEventProxy.wevent->genevtweight()/23443.423;  
  //if(myEventProxy.wevent->sample()==2) return myEventProxy.wevent->genevtweight()/225892.45;  
  //if(myEventProxy.wevent->sample()==3) return myEventProxy.wevent->genevtweight()/6383;

  /* Not using the HT samples at the moment.
  if(myEventProxy.wevent->sample()==2){
    //https://twiki.cern.ch/twiki/pub/CMS/HiggsToTauTauWorking2015/WplusHtWeights.xls
    //NLO to LO scaling removed, as NLO cross section used in normalisation.
    std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
    if(fileName.find("WJetsToLNu_HT-100To200")!=std::string::npos) return 0.1352710705/1.2137837838;
    else if(fileName.find("WJetsToLNu_HT-200To400")!=std::string::npos) return 0.076142149/1.2137837838;
    else if(fileName.find("WJetsToLNu_HT-400To600")!=std::string::npos) return 0.0326980819/1.2137837838;
    else if(fileName.find("WJetsToLNu_HT-600ToInf")!=std::string::npos) return 0.0213743732/1.2137837838;
    else return 0.8520862372/1.2137837838;
  }
  */
  
  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillControlHistos(const std::string & hNameSuffix, float eventWeight){

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
  myHistos_->fill1DHistogram("h1DnPCAMuon"+hNameSuffix,aMuon.nPCARefitvx().Mag(),eventWeight);

  ///Fill tau
  myHistos_->fill1DHistogram("h1DPtTau"+hNameSuffix,aTau.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+hNameSuffix,aTau.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+hNameSuffix,aTau.phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DIDTau"+hNameSuffix,aTau.tauID(byCombinedIsolationDeltaBetaCorrRaw3Hits) ,eventWeight);  
  myHistos_->fill1DHistogram("h1DStatsDecayMode"+hNameSuffix, aTau.decayMode(), eventWeight);
  myHistos_->fill1DHistogram("h1DnPCATau"+hNameSuffix,aTau.nPCARefitvx().Mag(),eventWeight);

  ///Fill leading tau track pt
  myHistos_->fill1DHistogram("h1DPtTauLeadingTk"+hNameSuffix,aTau.leadingTk().Pt(),eventWeight);
  ///Fill jets info           
  myHistos_->fill1DHistogram("h1DStatsNJets30"+hNameSuffix,nJets30,eventWeight);
  myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet.eta(),eventWeight);
  myHistos_->fill1DHistogram("h1DCSVBtagLeadingJet"+hNameSuffix,aJet.csvtag(),eventWeight);

  myHistos_->fill1DHistogram("h1DPtMET"+hNameSuffix,aMET.metpt(),eventWeight);

  fillDecayPlaneAngle(hNameSuffix, eventWeight);
  
  if(aJet.bjet()){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet.pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet.eta(),eventWeight);
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::fillVertices(const std::string & sysType){

  TVector3 aVertexGen = aEvent.genPV();
  TVector3 aVertex;
  if(sysType.find("AODPV")!=std::string::npos) aVertex = aEvent.thePV();
  if(sysType.find("RefitPV")!=std::string::npos) aVertex = aEvent.refitPfPV();
  
  float pullX = aVertexGen.X() - aVertex.X();
  float pullY = aVertexGen.Y() - aVertex.Y();
  float pullZ = aVertexGen.Z() - aVertex.Z();

  myHistos_->fill1DHistogram("h1DVxPullX_"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullY_"+sysType,pullY);
  myHistos_->fill1DHistogram("h1DVxPullZ_"+sysType,pullZ);

  myHistos_->fill2DHistogram("h2DVxPullVsNTrackTrans_"+sysType, aEvent.nTracksInRefit(), sqrt(pullX*pullX + pullY*pullY));
  myHistos_->fill2DHistogram("h2DVxPullVsNTrackLong_"+sysType, aEvent.nTracksInRefit(), pullZ);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight){

  ///Method from http://arxiv.org/abs/1108.0670 (S. Berge)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles;

  ///Angles using rho decay of a hadronic leg
  ///
  std::pair<float,float>  anglesIPRho;

  TLorentzVector positiveLeadingTk, negativeLeadingTk;
  TLorentzVector positive_nPCA, negative_nPCA;
  float yTau;

  fillVertices(hNameSuffix);

  if(aMuon.charge()>0){
    positiveLeadingTk = aMuon.leadingTk();
    positive_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);

    if(hNameSuffix.find("AODPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) positive_nPCA = TLorentzVector(aMuon.nPCAGenvx(),0);

    negativeLeadingTk = aTau.leadingTk();
    negative_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);
    if(hNameSuffix.find("AODPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) negative_nPCA = TLorentzVector(aTau.nPCAGenvx(),0);

    anglesIPRho = angleBetweenPlanes(aTau.leadingTk(), aTau.p4()-aTau.leadingTk(),
				     aMuon.leadingTk(), positive_nPCA);

  }
  else{    
    positiveLeadingTk = aTau.leadingTk();
    positive_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);

    if(hNameSuffix.find("AODPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) positive_nPCA = TLorentzVector(aTau.nPCAGenvx(),0);
    
    negativeLeadingTk = aMuon.leadingTk();
    negative_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);
    if(hNameSuffix.find("AODPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCA(),0);
    if(hNameSuffix.find("RefitPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCARefitvx(),0);
    if(hNameSuffix.find("GenPV")!=std::string::npos) negative_nPCA = TLorentzVector(aMuon.nPCAGenvx(),0);

    anglesIPRho = angleBetweenPlanes(aMuon.leadingTk(), negative_nPCA,
				     aTau.leadingTk(), aTau.p4()-aTau.leadingTk() );

  }
    
  angles = angleBetweenPlanes(negativeLeadingTk,negative_nPCA,
			      positiveLeadingTk,positive_nPCA);

  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  yTau =  2.*aTau.leadingTk().Pt()/aTau.pt() - 1.;
  float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;
 
  if(aTau.decayMode()!=tauDecay1ChargedPion0PiZero && isOneProng( aTau.decayMode())){
     myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);
     myHistos_->fill1DHistogram("h1DPhi_nVecIP_"+hNameSuffix,shiftedIPrho,eventWeight);
     if(yTau>0){
       myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
     else{
       myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
  }

  if(aTau.decayMode()==tauDecay1ChargedPion0PiZero){

    float deltaPhi  = fabs(positive_nPCA.Vect().Phi() - aGenPositiveTau.nPCA().Phi());
    float magnitudeRatio = (positive_nPCA.Vect().Mag()/aGenPositiveTau.nPCA().Mag());

    myHistos_->fill1DHistogram("h1DPhi_PCA_deltaPhi_"+hNameSuffix, deltaPhi);
    myHistos_->fill1DHistogram("h1DPhi_PCA_magnitudeRatio_"+hNameSuffix, magnitudeRatio);

    
    myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);
    float cosPositive =  positive_nPCA.Vect().Unit()*aGenPositiveTau.nPCA().Unit();
    myHistos_->fill1DHistogram("h1DCosPhi_CosPositive"+hNameSuffix,cosPositive,eventWeight);
    myHistos_->fillProfile("hProfPhiVsMag_"+hNameSuffix,positive_nPCA.Vect().Perp(),cosPositive);
    
    float cosNegative = negative_nPCA.Vect().Unit()*aGenNegativeTau.nPCA().Unit();
    myHistos_->fill1DHistogram("h1DCosPhi_CosNegative"+hNameSuffix,cosNegative,eventWeight);
    myHistos_->fillProfile("hProfPhiVsMag_"+hNameSuffix,negative_nPCA.Vect().Perp(),cosNegative);
    
    myHistos_->fillProfile("hProfPtVsMag_"+hNameSuffix,aGenNegativeTau.nPCA().Mag(), negativeLeadingTk.Perp());
    myHistos_->fillProfile("hProfPtVsMag_"+hNameSuffix,aGenPositiveTau.nPCA().Mag(), positiveLeadingTk.Perp());
    
    myHistos_->fillProfile("hProfMagVsPt_"+hNameSuffix, negativeLeadingTk.Perp(), aGenNegativeTau.nPCA().Mag());
    myHistos_->fillProfile("hProfMagVsPt_"+hNameSuffix, positiveLeadingTk.Perp(), aGenPositiveTau.nPCA().Mag());

    float cosPhiNN =  negative_nPCA.Vect().Unit().Dot(positive_nPCA.Vect().Unit());
    myHistos_->fill1DHistogram("h1DCosPhiNN_"+hNameSuffix,cosPhiNN);
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight){

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

  if(aGenPositiveTau.decayMode()==tauDecay1ChargedPion0PiZero ||
     aGenNegativeTau.decayMode()==tauDecay1ChargedPion0PiZero){
  
    float cosPhiNN =  negative_nPCA.Vect().Unit().Dot(positive_nPCA.Vect().Unit());
    myHistos_->fill1DHistogram("h1DCosPhiNN_"+hNameSuffix,cosPhiNN);
    myHistos_->fillProfile("hProfPtVsMag_"+hNameSuffix,aGenNegativeTau.nPCA().Mag(),negativeLeadingTk.Perp());
    myHistos_->fillProfile("hProfPtVsMag_"+hNameSuffix,aGenPositiveTau.nPCA().Mag(),positiveLeadingTk.Perp());
    myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);
  }
  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  std::pair<float,float>  anglesIPRho(-99,-99);
  
  float yTau;
  if(aGenPositiveTau.decayMode()!=tauDecay1ChargedPion0PiZero && isOneProng( aGenPositiveTau.decayMode() ) ){
    anglesIPRho = angleBetweenPlanes(aGenNegativeTau.leadingTk(), negative_nPCA,
				     aGenPositiveTau.leadingTk(), aGenPositiveTau.p4()-aGenPositiveTau.leadingTk() );
    yTau=2.*aGenPositiveTau.leadingTk().Pt()/aGenPositiveTau.pt() - 1.;
  }
  else if(aGenNegativeTau.decayMode()!=tauDecay1ChargedPion0PiZero && isOneProng( aGenNegativeTau.decayMode() ) ){
    anglesIPRho = angleBetweenPlanes(aGenNegativeTau.leadingTk(), aGenNegativeTau.p4()-aGenNegativeTau.leadingTk(),
				     aGenPositiveTau.leadingTk(), positive_nPCA);
    yTau=2.*aGenNegativeTau.leadingTk().Pt()/aGenNegativeTau.pt() - 1.;
  }
  else return;
  
  float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;
  
  myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);
  myHistos_->fill1DHistogram("h1DPhi_nVecIP_"+hNameSuffix,shiftedIPrho,eventWeight);

  if(yTau>0) myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
  else myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<float,float> HTTAnalyzer::angleBetweenPlanes(const TLorentzVector &tau1, 
						       const TLorentzVector &tau1Daughter,
						       const TLorentzVector &tau2, 
						       const TLorentzVector &tau2Daughter,
						       bool sgn){
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

  if(sgn){
    ///phase
    float calO = direction * ( n1.Cross(n2) );
    if(calO<0)
      phi = 2*TMath::Pi() - phi;
  }

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
  
  //for(auto it: decayNamesGen) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodGenDecayMode = true;
  //for(auto it: decayNamesReco) if(it.find("Lepton1Prong0Pi0")!=std::string::npos) goodRecoDecayMode = true;
  
  for(auto it: decayNamesGen) if(it.find("Lepton1Prong")!=std::string::npos) goodGenDecayMode = true;
  for(auto it: decayNamesReco) if(it.find("Lepton1Prong")!=std::string::npos) goodRecoDecayMode = true;

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

  if(goodGenDecayMode){
    fillGenDecayPlaneAngle(sampleName+"GenNoOfflineSel", eventWeight);
  }

  ///This stands for core selection, that is common to all regions.
  bool tauKinematics = aTau.pt()>30 && fabs(aTau.eta())<2.3;
  bool tauID = aTau.tauID(byMediumCombinedIsolationDeltaBetaCorr3Hits);  
  bool muonKinematics = aMuon.pt()>22 && fabs(aMuon.eta())<2.1;
  bool trigger = aPair.trigger(HLT_IsoMu17_eta2p1);
  if(sampleName=="Data") trigger = aPair.trigger(HLT_IsoMu18) || aPair.trigger(HLT_IsoMu20);

  bool cpMuonSelection = aMuon.nPCARefitvx().Perp()>0.003;    
  bool cpTauSelection = (aTau.decayMode()==tauDecay1ChargedPion0PiZero && aTau.nPCARefitvx().Perp()>0.003) ||
                        (aTau.decayMode()!=tauDecay1ChargedPion0PiZero && isOneProng(aTau.decayMode()));    

  bool cpSelection = aEvent.nTracksInRefit()>=2 && cpMuonSelection && cpTauSelection;

  bool extraRequirements = true;
  //extraRequirements &= nJets30==0;
  //extraRequirements &= aTau.decayMode()!=5 && aTau.decayMode()!=6 && nJets30==0;
  //extraRequirements &= (aPair.m_vis()>50 && aPair.m_vis()<100);
  
  if(!myEventProxy.wpair->size()) return true;
  if(!tauKinematics || !tauID || !muonKinematics || !trigger) return true;
  if(!extraRequirements) return true;
  if(!cpSelection) return true;
    
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
    fillControlHistos(hNameSuffix, eventWeight);
    if(goodGenDecayMode || sampleName.find("WJets")!=std::string::npos){
      fillDecayPlaneAngle(hNameSuffix+"RefitPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"AODPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"GenPV", eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"Gen", eventWeight);
    }
  }

  ///Histograms for the QCD control region
  if(qcdSelectionSS){
    hNameSuffix = sampleName+"qcdselSS";
    ///SS ans OS isolation histograms are filled only for mT<40 to remove possible contamination
    //from TT in high mT region.
    if(aMuon.mt()<40) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
    ///Fill SS histos in signal mu isolation region. Those histograms
    ///provide shapes for QCD estimate in signal region and in various control regions.
    ///If control region has OS we still use SS QCD estimate.
    if(aMuon.mt()<40 && aMuon.iso()<0.1) fillControlHistos(hNameSuffix, eventWeight);
    if(aMuon.mt()>60 && aMuon.iso()<0.1){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS",aMuon.mt(),eventWeight);    
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS",aMuon.mt(),eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"wselOS",eventWeight);
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
    if(aPair.diq()==-1){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aMuon.mt(),eventWeight);
      fillDecayPlaneAngle(hNameSuffix+"OS",eventWeight);
    }
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
