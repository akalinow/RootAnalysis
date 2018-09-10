#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "TF1.h"
#include "Math/LorentzVector.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinder.h"
#include "Math/Boost.h"
#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include "svfitAnalyzer.h"
#include "svfitHistograms.h"
#include "MLObjectMessenger.h"
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
		//std::cout<<aLeg1.getProperty(PropertyEnum::charge)<<std::endl;
		//std::cout<<aLeg1.getProperty(PropertyEnum::PDGId)<<std::endl;
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
 
  //Leptons for SVFit
  TVector3 recPV = aEvent.getRefittedPV();
  TVector3 recSV = aLeg1.getSV();  
  double cosGJReco = (recSV - recPV).Unit()*aLeg1.getP4().Vect().Unit();
  classic_svFit::MeasuredTauLepton aLepton1(type1, aLeg1.getP4().Pt(), aLeg1.getP4().Eta(),
					    aLeg1.getP4().Phi(), mass1, decay1);

  TVector3 genPV = aEvent.getGenPV();
  TVector3 genSV = aGenLeg1.getSV();
  TLorentzVector aP4 = aGenLeg1.getChargedP4();  
  double ip3D = aLeg1.getPCA().Mag();
 
  aLepton1.setCosGJ(cosGJReco);
  aLepton1.setIP3D(ip3D);

  classic_svFit::MeasuredTauLepton aLepton2(type2, aLeg2.getP4().Pt(), aLeg2.getP4().Eta(),
					    aLeg2.getP4().Phi(), mass2, decay2);
  recSV = aLeg2.getSV();  
  cosGJReco = (recSV - recPV).Unit()*aLeg2.getP4().Vect().Unit();
  ip3D = aLeg2.getPCA().Mag();

  genSV = aGenLeg2.getSV();
  aP4 = aGenLeg2.getChargedP4();
   
  aLepton2.setCosGJ(cosGJReco);
  aLepton2.setIP3D(ip3D);

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(aLepton1);
  measuredTauLeptons.push_back(aLepton2);
  
  //MET
	TVector2 aMET = aPair.getMET();
  TMatrixD covMET(2, 2);
  covMET[0][0] = aPair.getMETMatrix().at(0);
  covMET[0][1] = aPair.getMETMatrix().at(1);
  covMET[1][0] = aPair.getMETMatrix().at(2);
  covMET[1][1] = aPair.getMETMatrix().at(3); 
  if(covMET[0][0]==0 && covMET[1][0]==0 && covMET[0][1]==0 && covMET[1][1]==0) return TLorentzVector(); //singular covariance matrix
	else std::cout<<"non zero covMET"<<std::endl;

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

  aP4 = fastMTTAlgo.getTau1P4();
  svFitLeg1P4 = TLorentzVector(aP4.X(), aP4.Y(), aP4.Z(), aP4.E());
  
  aP4 = fastMTTAlgo.getTau2P4();
  svFitLeg2P4 = TLorentzVector(aP4.X(), aP4.Y(), aP4.Z(), aP4.E());
  
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

  svFitAlgo.integrate(measuredTauLeptons, aMET.X(), aMET.Y(), covMET);
  classic_svFit::DiTauSystemHistogramAdapter* diTauAdapter = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo.getHistogramAdapter());
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
std::tuple<double, double> svfitAnalyzer::getTauMomentum(const TLorentzVector & visP4, double cosGJ){

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
    return std::tuple<double, double>(0, 0);
  }

  double twoA = 2.0*(mVis2 + pVis2*sinGJ2);
  double solution1 = (b2 - sqrt(delta))/twoA;
  double solution2 = (b2 + sqrt(delta))/twoA;

  double tauEnergy1 = sqrt(pow(solution1,2) + mTau2);
  double tauEnergy2 = sqrt(pow(solution2,2) + mTau2);

  solution1 = tauEnergy1;
  solution2 = tauEnergy2;

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
    aGenLeg1.getChargedP4() - aGenLeg2.getChargedP4() -
    aGenLeg1.getNeutralP4() - aGenLeg2.getNeutralP4();
  TLorentzVector tautauGen = aGenLeg1.getP4() + aGenLeg2.getP4();

  myHistos_->fill1DHistogram("h1DMassGen"+hNameSuffix,tautauGen.M());

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
  
 
  //std::cout<<"mVisLeg1: "<<aLeg1.getP4().M()<<std::endl;
  //std::cout<<"mVisLeg2: "<<aLeg2.getP4().M()<<std::endl;

  TLorentzVector svFitP4 = aPair.getP4();
  double delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaStandalone"+hNameSuffix,delta);

  delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiStandalone"+hNameSuffix,delta);

  delta = (svFitP4.P() - tautauGen.P())/tautauGen.P();
  myHistos_->fill1DHistogram("h1DDeltaPStandalone"+hNameSuffix,delta);
  
  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtStandalone"+hNameSuffix,delta);

  delta = (svFitP4.Px() - tautauGen.Px())/tautauGen.Px(); 
  myHistos_->fill1DHistogram("h1DDeltaPxStandalone"+hNameSuffix,delta);

  delta = (svFitP4.Py() - tautauGen.Py())/tautauGen.Py(); 
  myHistos_->fill1DHistogram("h1DDeltaPyStandalone"+hNameSuffix,delta);

  delta = (svFitP4.Pz() - tautauGen.Pz())/tautauGen.Pz(); 
  myHistos_->fill1DHistogram("h1DDeltaPzStandalone"+hNameSuffix,delta);
  
  delta = (svFitP4.E() - tautauGen.E())/tautauGen.E();
  myHistos_->fill1DHistogram("h1DDeltaEStandalone"+hNameSuffix,delta);

  //svFitP4 = computeMTT("svfit");
  myHistos_->fill1DHistogram("h1DMassSVClassic"+hNameSuffix,svFitP4.M());
  myHistos_->fill1DHistogram("h1DCpuTimeClassic"+hNameSuffix,svFitAlgo.getComputingTime_cpu());
  delta = (svFitP4.Py() - tautauGen.Py())/tautauGen.Py();
  myHistos_->fill1DHistogram("h1DDeltaPtClassic"+hNameSuffix,delta);

  delta = (svFitP4.E() - tautauGen.E())/tautauGen.E();
  myHistos_->fill1DHistogram("h1DDeltaEClassic"+hNameSuffix,delta);
  
  svFitP4 = computeMTT("fastMTT");  
  myHistos_->fill1DHistogram("h1DMassSVFast"+hNameSuffix,svFitP4.M());
  myHistos_->fill1DHistogram("h1DCpuTimeFast"+hNameSuffix,fastMTTAlgo.getCpuTime("scan"));
  
  double mVis = (aLeg1.getP4() + aLeg2.getP4()).M();
  double superFast = mVis*exp(1/3.0)*0.96;
  myHistos_->fill1DHistogram("h1DMassSVSuperFast"+hNameSuffix, superFast);
  
  delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaFast"+hNameSuffix,delta);

   delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiFast"+hNameSuffix,delta);

  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtFast"+hNameSuffix,delta);

  delta = (svFitP4.P() - tautauGen.P())/tautauGen.P(); 
  myHistos_->fill1DHistogram("h1DDeltaPFast"+hNameSuffix,delta);

  delta = (svFitP4.Px() - tautauGen.Px())/tautauGen.Px();
  myHistos_->fill1DHistogram("h1DDeltaPxFast"+hNameSuffix,delta);

  delta = (svFitP4.Py() - tautauGen.Py())/tautauGen.Py(); 
  myHistos_->fill1DHistogram("h1DDeltaPyFast"+hNameSuffix,delta);

  delta = (svFitP4.Pz() - tautauGen.Pz())/tautauGen.Pz(); 
  myHistos_->fill1DHistogram("h1DDeltaPzFast"+hNameSuffix,delta);

  delta = (svFitP4.E() - tautauGen.E())/tautauGen.E(); 
  myHistos_->fill1DHistogram("h1DDeltaEFast"+hNameSuffix,delta);

  double recoX1 = 0.0, recoX2 = 0.0;
  std::tie(recoX1, recoX2) = fastMTTAlgo.getBestX();
  double bestLH = fastMTTAlgo.getBestLikelihood();
  
  delta = (svFitP4.M() - tautauGen.M())/svFitP4.M();
  myHistos_->fill2DHistogram("h2DDeltaM"+hNameSuffix, delta, -bestLH*10);

  myHistos_->fill2DHistogram("h2DMVisRatio"+hNameSuffix, mVis, svFitP4.M()/tautauGen.M());
    
  myHistos_->fill2DHistogram("h2DMHvsMVis"+hNameSuffix,  mVis, tautauGen.M());

  delta = aLeg1.getP4().E()/(aGenLeg1.getChargedP4() + aGenLeg1.getNeutralP4()).E();
  myHistos_->fill1DHistogram("h1DDeltaELeg1"+hNameSuffix, delta);
  
  delta = aLeg2.getP4().E()/(aGenLeg2.getChargedP4() + aGenLeg2.getNeutralP4()).E();
  myHistos_->fill1DHistogram("h1DDeltaELeg2"+hNameSuffix, delta);

  delta = (aMET.getP4().Perp() - nunuGen.Perp())/ nunuGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaMET"+hNameSuffix, delta);

  /////
  bool isThreeProng = aLeg2.getProperty(PropertyEnum::decayMode)>=HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion0PiZero &&
                      aLeg2.getProperty(PropertyEnum::decayMode)<HTTAnalysis::hadronicTauDecayModes::tauDecayMuon;

  TVector3 genSVLeg1 = aGenLeg1.getSV();
  delta = (genSVLeg1 - genPV).Mag();
  myHistos_->fill1DHistogram("h1DDelta_FlightPathLeg1_Gen"+hNameSuffix,delta*10);
  
  delta /= (aGenLeg1.getP4().Gamma()*aGenLeg1.getP4().Beta());
  myHistos_->fill1DHistogram("h1DDelta_FlightPathLeg1_Gen_CMS"+hNameSuffix,delta*10);
  /*
  TVector3 d = (genSVLeg1 - genPV);
  TVector3 pT = aGenLeg1.getChargedP4().Vect();

  double cosGJ = d.Unit().Dot(pT.Unit());
  double sinGJ = sqrt(1 - cosGJ*cosGJ);
  double pca1 = aGenLeg1.getPCA().Mag();
  double pca2 = d.Mag()*sinGJ;
  std::cout<<"d: "<<d.Mag()
	   <<" sinGJ: "<<sinGJ
	   <<" gamma*beta: "<<aGenLeg1.getP4().Gamma()*aGenLeg1.getP4().Beta()
	   <<std::endl;
  
  std::cout<<"pca1: "<<pca1<<" pca2: "<<pca2<<std::endl;
  */
  if(false && isThreeProng){

    TVector3 genSV = aGenLeg2.getSV();
    TLorentzVector a1P4 = aGenLeg2.getChargedP4() + aGenLeg2.getNeutralP4();
    double cosGJGen = (genSV - genPV).Unit()*a1P4.Vect().Unit();
    
    TVector3 recSV = aLeg2.getSV();  
    double cosGJReco = (recSV - recPV).Unit()*aLeg2.getP4().Vect().Unit();

    double delta = (recSV.X() - genSV.X());
    myHistos_->fill1DHistogram("h1DDeltaSV_X"+hNameSuffix,delta*10);

    delta = (recSV.Y() - genSV.Y());
    myHistos_->fill1DHistogram("h1DDeltaSV_Y"+hNameSuffix,delta*10);
    
    delta = (recSV.Z() - genSV.Z());
    myHistos_->fill1DHistogram("h1DDeltaSV_Z"+hNameSuffix,delta*10);

    delta = (recSV - recPV).Mag();
    myHistos_->fill1DHistogram("h1DDelta_FlightPathLeg2_Reco"+hNameSuffix,delta*10);

    delta = (genSV - genPV).Mag();
    myHistos_->fill1DHistogram("h1DDelta_FlightPathLeg2_Gen"+hNameSuffix,delta*10);
 
    std::tuple<double, double> gen = getTauMomentum(a1P4, cosGJGen);
    std::tuple<double, double> reco = getTauMomentum(aLeg2.getP4(), cosGJReco);

    myHistos_->fill1DHistogram("h1DCosGJGen"+hNameSuffix, 1E3*(cosGJGen-1));
    myHistos_->fill1DHistogram("h1DCosGJReco"+hNameSuffix, 1E3*(cosGJReco-1));
    
    delta = (aLeg2.getP4().M() - a1P4.M())/a1P4.M();
    myHistos_->fill1DHistogram("h1DDeltaMA1"+hNameSuffix,delta);

    //if(std::abs(delta)>0.1) return;//TEST

    delta = (aLeg2.getP4().P() - a1P4.P())/a1P4.P();
    myHistos_->fill1DHistogram("h1DDeltaPA1"+hNameSuffix,delta);

    delta = (1 - std::pow(cosGJReco,2) - (1 - std::pow(cosGJGen,2)))/(1 - std::pow(cosGJGen,2));
    myHistos_->fill1DHistogram("h1DDeltaSinGJ"+hNameSuffix,delta);

    double mTau = 1.77685;
    double mTau2 = std::pow(mTau,2);

    TLorentzVector visP4 = aLeg2.getP4();
    double cosGJ = cosGJReco;
      
    long double mVis =  visP4.M();
    long double mVis2 = std::pow(mVis,2);
    long double pVis =  visP4.P();
    long double pVis2 = std::pow(pVis,2);
    long double cosGJ2 = std::pow(cosGJ,2);
    long double sinGJ2 = 1.0 - cosGJ2;
    long double deltaReco = (mVis2 + pVis2)*(std::pow(mVis2 - mTau2, 2) - 4.0*mTau2*pVis2*sinGJ2);
    getTauMomentum(aLeg2.getP4(), cosGJReco);
         
    visP4 = aGenLeg2.getChargedP4() + aGenLeg2.getNeutralP4();
    cosGJ = cosGJGen;
    
    mVis =  visP4.M();
    mVis2 = std::pow(mVis,2);
    pVis =  visP4.P();
    pVis2 = std::pow(pVis,2);
    cosGJ2 = std::pow(cosGJ,2);
    sinGJ2 = 1.0 - cosGJ2;
    long double deltaGen = (mVis2 + pVis2)*(std::pow(mVis2 - mTau2, 2) - 4.0*mTau2*pVis2*sinGJ2);

    delta = (deltaReco - deltaGen)/deltaGen;
    myHistos_->fill1DHistogram("h1DDeltaDelta"+hNameSuffix, delta);

    myHistos_->fill1DHistogram("h1DDeltaReco"+hNameSuffix, deltaReco/1E5);
    myHistos_->fill1DHistogram("h1DDeltaGen"+hNameSuffix, deltaGen/1E5);

    double sinGJmax = (mTau2 - mVis2)/(2*mTau*pVis);
    delta = sqrt(sinGJ2)/sinGJmax;
    myHistos_->fill1DHistogram("h1DDeltaRatioGJMax"+hNameSuffix, delta);
        
    double recoLow = std::get<0>(reco);
    double recoHigh = std::get<1>(reco);

    double genLow = std::get<0>(gen);
    double genHigh = std::get<1>(gen);
        
    delta = (recoLow - genLow)/genLow;
    myHistos_->fill1DHistogram("h1DDeltaSolution1"+hNameSuffix,delta);

    delta = (recoHigh - genHigh)/genHigh;
    myHistos_->fill1DHistogram("h1DDeltaSolution2"+hNameSuffix,delta);
      
    delta = (svFitP4.M() - (aGenLeg1.getP4() + aGenLeg2.getP4()).M())/svFitP4.M();
    myHistos_->fill2DHistogram("h2DDeltaM"+hNameSuffix, svFitP4.M(), delta);

    visP4 = aGenLeg2.getChargedP4() + aGenLeg2.getNeutralP4();
    bool isGenSolution1 = (-aGenLeg2.getP4().Beta() + visP4.Beta()*cosGJGen)>0;

    //if(isGenSolution1) return;//TEST

    visP4 = aLeg2.getP4();
    double tauE = fastMTTAlgo.getTau2P4().E();
    double tauP = sqrt(tauE*tauE - classic_svFit::tauLeptonMass2); 
    double betaTau = tauP/tauE; 
    bool isRecoSolution1 = (-betaTau + visP4.Beta()*cosGJ)>0;

    if(isGenSolution1){
       myHistos_->fill1DHistogram("h1DDeltaCorrectSolution1",isGenSolution1&&isRecoSolution1);
       myHistos_->fill1DHistogram("h1DDeltaMETRatioS1",(1- betaTau)/(1-aGenLeg2.getP4().Beta()));
    }
    if(!isGenSolution1){
      myHistos_->fill1DHistogram("h1DDeltaCorrectSolution2",(!isGenSolution1)&&(!isRecoSolution1));
      myHistos_->fill1DHistogram("h1DDeltaMETRatioS2",(1- betaTau)/(1-aGenLeg2.getP4().Beta()));
    }
     
    delta = (fastMTTAlgo.getTau1P4().E() - aGenLeg1.getP4().E())/aGenLeg1.getP4().E();   
    myHistos_->fill1DHistogram("h1DDeltaLeg1"+hNameSuffix, delta);
    delta = (fastMTTAlgo.getTau2P4().E() - aGenLeg2.getP4().E())/aGenLeg2.getP4().E();
    myHistos_->fill1DHistogram("h1DDeltaLeg2"+hNameSuffix, delta);

    delta = (svFitP4.M() - tautauGen.M())/tautauGen.M();
    myHistos_->fill1DHistogram("h1DDeltaMassResolutionFast"+hNameSuffix, delta); 

    myHistos_->fill1DHistogram("h1DMassSV3ProngFast"+hNameSuffix, svFitP4.M());
    
    fastMTTAlgo.enableComponent(fastMTT::ENERGY);
    //fastMTTAlgo.disableComponent(fastMTT::MET);
    //svFitP4 = computeMTT("fastMTT");
    fastMTTAlgo.disableComponent(fastMTT::ENERGY);
    fastMTTAlgo.enableComponent(fastMTT::MET);

    myHistos_->fill1DHistogram("h1DMassSV3ProngFastE"+hNameSuffix, svFitP4.M());
    myHistos_->fill1DHistogram("h1DMassSV3ProngStandalone"+hNameSuffix, aPair.getP4().M()); 
  
    delta = (fastMTTAlgo.getTau1P4().E() - aGenLeg1.getP4().E())/aGenLeg1.getP4().E();   
    myHistos_->fill1DHistogram("h1DDeltaLeg1E"+hNameSuffix, delta);
    delta = (fastMTTAlgo.getTau2P4().E() - aGenLeg2.getP4().E())/aGenLeg2.getP4().E();
    myHistos_->fill1DHistogram("h1DDeltaLeg2E"+hNameSuffix, delta);

    delta = (aPair.getP4().M() - tautauGen.M())/tautauGen.M();
    myHistos_->fill1DHistogram("h1DDeltaMassResolutionStandalone"+hNameSuffix,delta);
    delta = (svFitP4.M() - tautauGen.M())/tautauGen.M();
    myHistos_->fill1DHistogram("h1DDeltaMassResolutionFastE"+hNameSuffix,delta);     

    double trueX1 = aLeg1.getP4().E()/aGenLeg1.getP4().E();
    double trueX2 = aLeg2.getP4().E()/aGenLeg2.getP4().E();
    double deltaX1, deltaX2, lhTest;
    double testX[2];
    int nGridPoints = 200;
    
    testX[0] = recoX1;
    testX[1] = recoX2;
 
    //fastMTTAlgo.enableComponent(fastMTT::ENERGY);
    //fastMTTAlgo.disableComponent(fastMTT::MASS);
    //fastMTTAlgo.disableComponent(fastMTT::MET);
    //fastMTTAlgo.disableComponent(fastMTT::PX);
    //fastMTTAlgo.disableComponent(fastMTT::PY);
    for(int iX2 = 1; iX2<nGridPoints;++iX2){
      testX[1] = 1.0*(double)iX2/nGridPoints;
      for(int iX1 = 1; iX1<nGridPoints;++iX1){
	testX[0] = 1.0*(double)iX1/nGridPoints;
	
	deltaX1 = testX[0] - trueX1;
	deltaX2 = testX[1] - trueX2;

	myHistos_->fill2DHistogram("h2DLikelihoodMapNorm"+hNameSuffix, deltaX1, deltaX2, 1.0);

	fastMTTAlgo.enableComponent(fastMTT::ENERGY);//TEST
	lhTest = fastMTTAlgo.getLikelihoodForX(testX);
	myHistos_->fill2DHistogram("h2DLikelihoodMap"+hNameSuffix, deltaX1, deltaX2, lhTest);

	fastMTTAlgo.enableComponent(fastMTT::ENERGY);
	fastMTTAlgo.disableComponent(fastMTT::MET);
	fastMTTAlgo.disableComponent(fastMTT::MASS);
	fastMTTAlgo.disableComponent(fastMTT::PX);
	fastMTTAlgo.disableComponent(fastMTT::PY);
	lhTest = fastMTTAlgo.getLikelihoodForX(testX);
	myHistos_->fill2DHistogram("h2DLikelihoodMapE"+hNameSuffix, deltaX1, deltaX2, lhTest);
	fastMTTAlgo.disableComponent(fastMTT::ENERGY);
	fastMTTAlgo.enableComponent(fastMTT::MET);
	fastMTTAlgo.enableComponent(fastMTT::MASS);
	fastMTTAlgo.enableComponent(fastMTT::PX);
	fastMTTAlgo.enableComponent(fastMTT::PY);
      }
    }
    fastMTTAlgo.disableComponent(fastMTT::ENERGY);
    fastMTTAlgo.enableComponent(fastMTT::MASS);
    fastMTTAlgo.enableComponent(fastMTT::MET);
    fastMTTAlgo.enableComponent(fastMTT::PX);
    fastMTTAlgo.enableComponent(fastMTT::PY);
    //////////////////////
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
bool svfitAnalyzer::analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){

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

  isGoodReco = aGenLeg2.getP4().DeltaR(aLeg2.getP4())<0.4;
  goodGenTau = aGenLeg2.getP4().E()>1.0;

  if(sampleName=="WAllJets"){
    goodGenTau = true;
    isGoodReco = true;
  }
  
  if(sampleName=="DYAllJetsMatchL"){
    goodGenTau = true;
    isGoodReco = true;
  }
  
  if(isGoodReco && goodGenTau){
    fillControlHistos(hNameSuffix);
  }
	aVisSumM = (aLeg1.getP4() + aLeg2.getP4()).M();
	aGenSumM = (aGenLeg1.getP4() + aGenLeg2.getP4()).M();
	if(aMessenger and std::string("MLObjectMessenger").compare((aMessenger->name()).substr(0,17))==0 ) // if NULL it will do nothing
	{
		// Putting data to MLObjectMessenger
		try
		{
			MLObjectMessenger* mess = (MLObjectMessenger*)aMessenger;
			mess->putObjectVector(const_cast <const HTTParticle*>(&aLeg1), "legs");
			mess->putObjectVector(const_cast <const HTTParticle*>(&aLeg2), "legs");
			mess->putObject(const_cast <const HTTParticle*>(&aMET), "amet");
			mess->putObject(&aPair.getMETMatrix().at(0), "covMET00");
			mess->putObject(&aPair.getMETMatrix().at(1), "covMET01");
			mess->putObject(&aPair.getMETMatrix().at(2), "covMET10");
			mess->putObject(&aPair.getMETMatrix().at(3), "covMET11");
			//float beta_score = aBJet1.getProperty(PropertyEnum::bCSVscore);
			//mess->putObject(beta_score, "beta_score");
			mess->putObject(&aVisSumM, "visible_mass");
			mess->putObject(&aGenSumM, "gen_visible_mass");
			//mess->putObject(&aPair.getMTMuon(), "higgs_mass_trans");
		}
		catch(const std::exception& e)
		{
			std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN "+std::string( __func__ )+ " WHEN PUTTING DATA TO MLObjectMessenger!"));
		}                    
	}

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
