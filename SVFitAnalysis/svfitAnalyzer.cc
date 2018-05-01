#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "TF1.h"
#include "Math/LorentzVector.h"

#include "svfitAnalyzer.h"
#include "svfitHistograms.h"
#include "MuTauSpecifics.h"
#include "TauTauSpecifics.h"
#include "MuMuSpecifics.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
svfitAnalyzer::svfitAnalyzer(const std::string & aName, const std::string & aDecayMode) : Analyzer(aName){

#pragma omp critical
  {
    if(aDecayMode=="MuTau") myChannelSpecifics = new MuTauSpecifics(this);
    else if (aDecayMode=="TauTau") myChannelSpecifics = new TauTauSpecifics(this);
    else if (aDecayMode=="MuMu") myChannelSpecifics = new MuMuSpecifics(this);
    myNumberOfCategories = myChannelSpecifics->getCategoryRejester().size();
    categoryDecisions.resize(myNumberOfCategories);

    ntupleFile_ = 0;

    svFitAlgo.addLogM_fixed(true, 4.0);
    svFitAlgo.setLikelihoodFileName("");
    svFitAlgo.setMaxObjFunctionCalls(100000);
    svFitAlgo.setVerbosity(1);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
svfitAnalyzer::~svfitAnalyzer(){

  if(myHistos_) delete myHistos_;
  if(myChannelSpecifics) delete myChannelSpecifics;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* svfitAnalyzer::clone() const {

  std::string myDecayMode = myChannelSpecifics->getDecayModeName();
  svfitAnalyzer* clone = new svfitAnalyzer(name(),myDecayMode);
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::initialize(TDirectory* aDir,
			       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  myHistos_ = new svfitHistograms(aDir, selectionFlavours_, myChannelSpecifics->getDecayModeName());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::finalize(){

  myHistos_->finalizeHistograms(myChannelSpecifics->getCategoryRejester());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  aEvent = *myEventProxy.event;
  aPair = (*myEventProxy.pairs)[0];

  TLorentzVector met4v(aPair.getMET().X(),
		       aPair.getMET().Y(),
		       0,
		       aPair.getMET().Mod());

  aMET = HTTParticle();
  aMET.setP4(met4v);
  myChannelSpecifics->setAnalysisObjects(myEventProxy);

  aSeparatedJets = getSeparatedJets(myEventProxy, 0.5);
  aJet1 = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
  aJet2 = aSeparatedJets.size()>1 ? aSeparatedJets[1] : HTTParticle();
  aBJet1 = HTTParticle();
  for(auto itJet: aSeparatedJets) {
    if(std::abs(itJet.getP4().Eta())<2.4 &&
       itJet.getProperty(PropertyEnum::bCSVscore)>0.8484){
      aBJet1 = itJet;
      break;
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> svfitAnalyzer::getSeparatedJets(const EventProxyHTT & myEventProxy,
							 float deltaR){

  std::vector<HTTParticle> separatedJets;

  for(auto aJet : *myEventProxy.jets) {
    float dRLeg2 = aJet.getP4().DeltaR(aLeg2.getP4());
    float dRLeg1 = aJet.getP4().DeltaR(aLeg1.getP4());
    bool loosePFJetID = aJet.getProperty(PropertyEnum::PFjetID)>=1;
    bool jetEtaCut = std::abs(aJet.getP4().Eta())<4.7;
    if(dRLeg1>deltaR && dRLeg2>deltaR && loosePFJetID && jetEtaCut) separatedJets.push_back(aJet);
  }
  return separatedJets;
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::addBranch(TTree *tree){ /*tree->Branch("muonPt",&muonPt);*/}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::computeMTT(const std::string & algoName){

  //Legs
  double mass1;
  int decay1 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type1;
  if(std::abs(aLeg1.getPDGid())==11) {
    mass1 = 0.51100e-3; //electron mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(aLeg1.getPDGid())==13) {
    mass1 = 0.10566; //muon mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay1 = aLeg1.getProperty(PropertyEnum::decayMode);
    mass1 = aLeg1.getP4().M();
    if(decay1==0)
      mass1 = 0.13957; //pi+/- mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }

  double mass2;
  int decay2 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type2;
  if(std::abs(aLeg2.getPDGid())==11) {
    mass2 = 0.51100e-3; //electron mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(aLeg2.getPDGid())==13) {
    mass2 = 0.10566; //muon mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay2 = aLeg2.getProperty(PropertyEnum::decayMode);
    mass2 = aLeg2.getP4().M();
    if(decay2==0) mass2 = 0.13957; //pi+/- mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }
  ///TEST 
  //mass1 = 0.51100e-3;
  //type1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;

  //decay2 = 0;
  //mass2 = 0.51100e-3;
  //type2 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;

  //decay1 = 0;
  //mass1 = 0.13957;
  //type1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  ///////
  

  //Leptons for SVFit
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type1, aLeg1.getP4().Pt(), aLeg1.getP4().Eta(),
								aLeg1.getP4().Phi(), mass1, decay1) );
  measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type2, aLeg2.getP4().Pt(), aLeg2.getP4().Eta(),
								aLeg2.getP4().Phi(), mass2, decay2) );
  //MET
  TVector2 aMET = aPair.getMET();
  TMatrixD covMET(2, 2);
  covMET[0][0] = aPair.getMETMatrix().at(0);
  covMET[0][1] = aPair.getMETMatrix().at(1);
  covMET[1][0] = aPair.getMETMatrix().at(2);
  covMET[1][1] = aPair.getMETMatrix().at(3);

  if(covMET[0][0]==0 && covMET[1][0]==0 && covMET[0][1]==0 && covMET[1][1]==0) return TLorentzVector(); //singular covariance matrix

  TLorentzVector aResult;
  if(algoName=="svfit") aResult = runSVFitAlgo(measuredTauLeptons, aMET, covMET);
  if(algoName=="fastMTT") aResult = runFastMTTAlgo(measuredTauLeptons, aMET, covMET);
  
  return aResult;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::runFastMTTAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
                                           const TVector2 &aMET, const TMatrixD &covMET){

  fastMTTAlgo.run(measuredTauLeptons, aMET.X(), aMET.Y(), covMET);
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > aP4 = fastMTTAlgo.getBestP4();
  
  TLorentzVector p4SVFit(aP4.X(), aP4.Y(), aP4.Z(), aP4.E());
   
  return p4SVFit;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
                                           const TVector2 &aMET, const TMatrixD &covMET){

  svFitAlgo.setVerbosity(0);

  svFitAlgo.addLogM_fixed(true, 4.0);
  if(myChannelSpecifics->getDecayModeName()=="MuTau") svFitAlgo.addLogM_fixed(true, 4.0);
  else if(myChannelSpecifics->getDecayModeName()=="TauTau") svFitAlgo.addLogM_fixed(true, 5.0);
  else if(myChannelSpecifics->getDecayModeName()=="MuMu") svFitAlgo.addLogM_fixed(true, 3.0);
  
  svFitAlgo.setMaxObjFunctionCalls(100000);
  //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  //svFitAlgo.setTreeFileName("markovChainTree.root");
  TVector3 recPV = aEvent.getRefittedPV();
  TVector3 svLeg2 = aLeg2.getSV();
  svLeg2 -=recPV;

  TVector3 svLeg1 = aGenLeg1.getSV();
  svLeg1 -=recPV;

  //svFitAlgo.addSVData(svLeg1, svLeg2);

  svFitAlgo.integrate(measuredTauLeptons, aMET.X(), aMET.Y(), covMET);
  classic_svFit::DiTauSystemHistogramAdapter* diTauAdapter = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter());
  //float cpuTime = svFitAlgo.getComputingTime_cpu();
  float mcMass = diTauAdapter->getMass();

  TLorentzVector p4SVFit;
  if(svFitAlgo.isValidSolution() )
    {//Get solution
      p4SVFit.SetPtEtaPhiM(diTauAdapter->getPt(),
			   diTauAdapter->getEta(),
			   diTauAdapter->getPhi(),
			   mcMass);
    }
  return p4SVFit;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::tuple<double, double> svfitAnalyzer::getTauMomentum(const TLorentzVector & visP4, const double &cosGJ){

  double mTau = 1.77685;
  double mTau2 = std::pow(mTau,2);
  
  double mVis =  visP4.M();
  double mVis2 = std::pow(mVis,2);
  double pVis =  visP4.P();
  double pVis2 = std::pow(pVis,2);
  double cosGJ2 = std::pow(cosGJ,2);
  double sinGJ2 = 1.0 - cosGJ2;

  double b2 = (mVis2 + mTau2)*pVis*cosGJ;

  double delta = (mVis2 + pVis2)*(std::pow(mVis2 - mTau2, 2) - 4.0*mTau2*pVis2*sinGJ2);
  if(delta<0){
    //std::cout<<"mVis2: "<<mVis2<<std::endl;
    return  std::tuple<double, double>(0, 0);
  }

  double twoA = 2.0*(mVis2 + pVis2*sinGJ2);
  
  double solution1 = (b2 - sqrt(delta))/twoA;

  double solution2 = (b2 + sqrt(delta))/twoA;
  /*
  std::cout<<" decayMode: "<<aGenLeg2.getProperty(PropertyEnum::decayMode)
           <<" solution1: "<<solution1
	   <<" solution2: "<<solution2
	   <<" true value: "<<aGenLeg2.getP4().P()
	   <<std::endl;
  */

  return std::tuple<double, double>(solution1, solution2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::fillControlHistos(const std::string & hNameSuffix){

  const TLorentzVector & aVisSum = aLeg1.getP4() + aLeg2.getP4();
  float visMass = aVisSum.M();

  myHistos_->fill1DHistogram("h1DMassSVStandalone"+hNameSuffix,aPair.getP4().M());
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix, visMass);
  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aPair.getMTMuon());

  TVector3 genPV = aEvent.getGenPV();
  TVector3 recPV = aEvent.getRefittedPV();

  ////
  //MET
  TMatrixD covMET(2, 2);
  covMET[0][0] = aPair.getMETMatrix().at(0);
  covMET[0][1] = aPair.getMETMatrix().at(1);
  covMET[1][0] = aPair.getMETMatrix().at(2);
  covMET[1][1] = aPair.getMETMatrix().at(3);


  TLorentzVector nunuGen = aGenLeg1.getP4() + aGenLeg2.getP4() -
    aGenLeg1.getChargedP4() - aGenLeg2.getChargedP4();
  TLorentzVector tautauGen = aGenLeg1.getP4() + aGenLeg2.getP4();

  /*
  TMatrixD A(2,2);
  A[0][0] = sin(aGenLeg1.getChargedP4().Theta())*
    cos(aGenLeg1.getChargedP4().Phi());

  A[1][1] = sin(aGenLeg2.getChargedP4().Theta())*
    sin(aGenLeg2.getChargedP4().Phi());

  A[0][1] = sin(aGenLeg2.getChargedP4().Theta())*
    cos(aGenLeg2.getChargedP4().Phi());

  A[1][0] = sin(aGenLeg1.getChargedP4().Theta())*
    sin(aGenLeg1.getChargedP4().Phi());

  TMatrixD invA = A.Invert();

  double e1 = invA[0][0]*nunuGen.Px() +
    invA[0][1]*nunuGen.Py();

  double e2 = invA[1][0]*nunuGen.Px() +
    invA[1][1]*nunuGen.Py();

  //std::cout<<"E1: "<<e1<<" nu1 e:"<< (aGenLeg1.getP4() -  aGenLeg1.getChargedP4()).E()<<std::endl;
  //std::cout<<"E2: "<<e2<<" nu2 e:"<< (aGenLeg2.getP4() -  aGenLeg2.getChargedP4()).E()<<std::endl;

  double x1 =  aGenLeg1.getChargedP4().E()/(aGenLeg1.getChargedP4().E() + e1);
  double x2 =  aGenLeg2.getChargedP4().E()/(aGenLeg2.getChargedP4().E() + e2);
  double caMassGen = (aGenLeg1.getChargedP4() + aGenLeg2.getChargedP4()).M()/std::sqrt(x1*x2);  
  myHistos_->fill1DHistogram("h1DMassSVCA"+hNameSuffix, caMassGen);
  */
  
  /*
  std::cout<<"mVis: "<<(aLeg1.getP4() + aLeg2.getP4()).M()<<std::endl;
  std::cout<<"mVisLeg1: "<<aLeg1.getP4().M()<<std::endl;
  std::cout<<"mVisLeg2: "<<aLeg2.getP4().M()<<std::endl;
  */
  TLorentzVector svFitP4 = computeMTT("fastMTT");
  
  myHistos_->fill1DHistogram("h1DMassSVFast"+hNameSuffix,svFitP4.M());
  myHistos_->fill1DHistogram("h1DCpuTimeFast"+hNameSuffix,fastMTTAlgo.getCpuTime("scan"));

  double delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaFast"+hNameSuffix,delta);

  delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiFast"+hNameSuffix,delta);

  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp(); 
  myHistos_->fill1DHistogram("h1DDeltaPtFast"+hNameSuffix,delta);
    
  //svFitP4 = computeMTT("svfit");
  myHistos_->fill1DHistogram("h1DMassSVClassic"+hNameSuffix,svFitP4.M());
  myHistos_->fill1DHistogram("h1DCpuTimeClassic"+hNameSuffix,svFitAlgo.getComputingTime_cpu());
  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtClassic"+hNameSuffix,delta);
  
  svFitP4 = aPair.getP4();
  delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaStandalone"+hNameSuffix,delta);

  delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiStandalone"+hNameSuffix,delta);

  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtStandalone"+hNameSuffix,delta);
  /////



  TVector3 genSV = aGenLeg2.getSV();
  TLorentzVector a1P4 = aGenLeg2.getChargedP4() + aGenLeg2.getNeutralP4();
  double cosGJGen = (genSV - genPV).Unit()*a1P4.Vect().Unit();

  TVector3 recSV = aLeg2.getSV();  
  double cosGJReco = (recSV - recPV).Unit()*aLeg2.getP4().Vect().Unit();

  bool isThreeProng = aLeg2.getProperty(PropertyEnum::decayMode)>=HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion0PiZero &&
                      aLeg2.getProperty(PropertyEnum::decayMode)<HTTAnalysis::hadronicTauDecayModes::tauDecayMuon;

  if(isThreeProng){
 
    std::tuple<double, double> gen = getTauMomentum(a1P4, cosGJGen);
    std::tuple<double, double> reco = getTauMomentum(aLeg2.getP4(), cosGJReco);

    if(std::get<0>(reco)<1E-3){
      delta = (aLeg2.getP4().M() - a1P4.M())/a1P4.M();
      myHistos_->fill1DHistogram("h1DDeltaMA1"+hNameSuffix,delta);

      delta = (aLeg2.getP4().P() - a1P4.P())/a1P4.P();
      myHistos_->fill1DHistogram("h1DDeltaPA1"+hNameSuffix,delta);

      delta = (std::pow(cosGJReco,2) - std::pow(cosGJGen,2))/(1 - std::pow(cosGJGen,2));
      delta *= 0.5;    
      myHistos_->fill1DHistogram("h1DDeltaCosGJ"+hNameSuffix,delta);  
    }


    double genTmp = std::max(std::get<0>(gen), std::get<1>(gen));
    double recoTmp = std::max(std::get<0>(reco), std::get<1>(reco));  
    
    delta = (recoTmp - genTmp)/genTmp;   
    myHistos_->fill1DHistogram("h1DDeltaSolution1"+hNameSuffix,delta);

    genTmp = std::min(std::get<0>(gen), std::get<1>(gen));
    recoTmp = std::min(std::get<0>(reco), std::get<1>(reco));  
    
    delta = (recoTmp - genTmp)/genTmp;   
    myHistos_->fill1DHistogram("h1DDeltaSolution2"+hNameSuffix,delta);
    
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool svfitAnalyzer::passCategory(unsigned int iCategory){

  if(categoryDecisions.size()==0) return false;
  else return categoryDecisions[iCategory];

  return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool svfitAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  sampleName = getSampleName(myEventProxy);

  std::string hNameSuffix = sampleName;

  if(!myEventProxy.pairs->size()) return true;
  setAnalysisObjects(myEventProxy);

  bool isGoodReco = aGenLeg1.getP4().DeltaR(aLeg1.getP4())<0.4 &&
		    aGenLeg2.getP4().DeltaR(aLeg2.getP4())<0.4;


  isGoodReco |= aGenLeg2.getP4().DeltaR(aLeg1.getP4())<0.4 &&
		aGenLeg1.getP4().DeltaR(aLeg2.getP4())<0.4;
  
  bool goodGenTau = aGenLeg1.getP4().E()>1.0 && aGenLeg2.getP4().E()>1.0;

 
  //double delta = aLeg2.getP4().E() - aGenLeg2.getChargedP4().E();
  //delta /= aGenLeg2.getChargedP4().E();
  //isGoodReco &= std::abs(delta)<0.1;

  if(sampleName=="WAllJets"){
    goodGenTau = true;
    isGoodReco = true;
  }

  if(isGoodReco && goodGenTau){
    fillControlHistos(hNameSuffix);
  }

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
