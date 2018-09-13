#include "SVfitAnalyzer.h"

#include "Tools.h"
#include "Math/LorentzVector.h"
#include "MLObjectMessenger.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
SVfitAnalyzer::SVfitAnalyzer(
    const std::string & aName,
    const std::string & aDecayMode)
  : Analyzer(aName), decayMode_(aDecayMode), covMET_(2,2){

#pragma omp critical
    {
      SVfitAlgo_.addLogM_fixed(true, 4.0);
      SVfitAlgo_.setLikelihoodFileName("");
      SVfitAlgo_.setMaxObjFunctionCalls(100000);
      SVfitAlgo_.setVerbosity(1);
    }
  }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
SVfitAnalyzer::~SVfitAnalyzer(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* SVfitAnalyzer::clone() const {

  SVfitAnalyzer* clone = new SVfitAnalyzer(name(),decayMode_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void SVfitAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  event_ = *myEventProxy.event;
  pair_ = (*myEventProxy.pairs)[0];

  //TLorentzVector met4v(pair_.getMET().X(),
  //    pair_.getMET().Y(),
  //    0,
  //    pair_.getMET().Mod());
  //MET_.setP4(met4v);

  //TODO use genMET() when it contains actual value instead of nunuGen
  TLorentzVector nunuGen = genLeg1_.getP4() + genLeg2_.getP4() -
    genLeg1_.getChargedP4() - genLeg2_.getChargedP4() -
    genLeg1_.getNeutralP4() - genLeg2_.getNeutralP4();
  MET_.setP4(nunuGen);

  covMET_[0][0] = pair_.getMETMatrix().at(0);
  covMET_[0][1] = pair_.getMETMatrix().at(1);
  covMET_[1][0] = pair_.getMETMatrix().at(2);
  covMET_[1][1] = pair_.getMETMatrix().at(3);

  separatedJets_ = HTTAnalysis::getSeparatedJets(myEventProxy, leg1_, leg2_, 0.5);
  jet1_ = separatedJets_.size() ? separatedJets_[0] : HTTParticle();
  jet2_ = separatedJets_.size()>1 ? separatedJets_[1] : HTTParticle();

  genSumM_ = (genLeg1_.getP4() + genLeg2_.getP4()).M();
  higgsMassTrans_ = pair_.getMTMuon();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TLorentzVector SVfitAnalyzer::computeMTT(const std::string & algoName){

  //Legs
  double mass1;
  int decay1 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type1;
  if(std::abs(leg1_.getPDGid())==11) {
    mass1 = 0.51100e-3; //electron mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(leg1_.getPDGid())==13) {
    mass1 = 0.10566; //muon mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay1 = leg1_.getProperty(PropertyEnum::decayMode);
     mass1 = leg1_.getP4().M();
    if(decay1==0)
      mass1 = 0.13957; //pi+/- mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }

  double mass2;
  int decay2 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type2;
  if(std::abs(leg2_.getPDGid())==11) {
    mass2 = 0.51100e-3; //electron mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(leg2_.getPDGid())==13) {
    mass2 = 0.10566; //muon mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay2 = leg2_.getProperty(PropertyEnum::decayMode);
    mass2 = leg2_.getP4().M();
    if(decay2==0) mass2 = 0.13957; //pi+/- mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }

  //Leptons for SVfit
  TVector3 recPV = event_.getRefittedPV();
  TVector3 recSV = leg1_.getSV();  
  double cosGJReco = (recSV - recPV).Unit()*leg1_.getP4().Vect().Unit();
  classic_svFit::MeasuredTauLepton aLepton1(type1, leg1_.getP4().Pt(), leg1_.getP4().Eta(),
      leg1_.getP4().Phi(), mass1, decay1);

  TVector3 genPV = event_.getGenPV();
  TVector3 genSV = genLeg1_.getSV();
  TLorentzVector aP4 = genLeg1_.getChargedP4();  
  double ip3D = leg1_.getPCA().Mag();

  aLepton1.setCosGJ(cosGJReco);
  aLepton1.setIP3D(ip3D);

  classic_svFit::MeasuredTauLepton aLepton2(type2, leg2_.getP4().Pt(), leg2_.getP4().Eta(),
      leg2_.getP4().Phi(), mass2, decay2);
  recSV = leg2_.getSV();  
  cosGJReco = (recSV - recPV).Unit()*leg2_.getP4().Vect().Unit();
  ip3D = leg2_.getPCA().Mag();

  genSV = genLeg2_.getSV();
  aP4 = genLeg2_.getChargedP4();

  aLepton2.setCosGJ(cosGJReco);
  aLepton2.setIP3D(ip3D);

  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(aLepton1);
  measuredTauLeptons.push_back(aLepton2);

  //MET
  if(covMET_[0][0]==0 && covMET_[1][0]==0 && covMET_[0][1]==0 && covMET_[1][1]==0) return TLorentzVector(); //singular covariance matrix

  TLorentzVector aResult;
  if(algoName=="SVfit") aResult = runSVfitAlgo(measuredTauLeptons);
  if(algoName=="fastMTT") aResult = runFastMTTAlgo(measuredTauLeptons);

  return aResult;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector SVfitAnalyzer::runFastMTTAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons){

  fastMTTAlgo_.run(measuredTauLeptons, MET_.getP4().X(), MET_.getP4().Y(), covMET_);
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > aP4 = fastMTTAlgo_.getBestP4();

  TLorentzVector p4SVfit(aP4.X(), aP4.Y(), aP4.Z(), aP4.E());

  return p4SVfit;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector SVfitAnalyzer::runSVfitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons){

  SVfitAlgo_.setVerbosity(0);

  SVfitAlgo_.addLogM_fixed(true, 4.0);
  if     (decayMode_=="MuTau" ) SVfitAlgo_.addLogM_fixed(true, 4.0);
  else if(decayMode_=="TauTau") SVfitAlgo_.addLogM_fixed(true, 5.0);
  else if(decayMode_=="MuMu"  ) SVfitAlgo_.addLogM_fixed(true, 3.0);

  SVfitAlgo_.setMaxObjFunctionCalls(100000);
  //SVfitAlgo_.setLikelihoodFileName("testClassicSVfit.root");
  //SVfitAlgo_.setTreeFileName("markovChainTree.root");
  TVector3 recPV = event_.getRefittedPV();
  TVector3 svLeg2 = leg2_.getSV();
  svLeg2 -=recPV;

  TVector3 svLeg1 = genLeg1_.getSV();
  svLeg1 -=recPV;

  SVfitAlgo_.integrate(measuredTauLeptons, MET_.getP4().X(), MET_.getP4().Y(), covMET_);
  classic_svFit::DiTauSystemHistogramAdapter* diTauAdapter = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(SVfitAlgo_.getHistogramAdapter());
  float mcMass = diTauAdapter->getMass();

  TLorentzVector p4SVfit;
  if(SVfitAlgo_.isValidSolution() )
  {//Get solution
    p4SVfit.SetPtEtaPhiM(diTauAdapter->getPt(),
        diTauAdapter->getEta(),
        diTauAdapter->getPhi(),
        mcMass);
  }
  return p4SVfit;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool SVfitAnalyzer::analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  sampleName_ = HTTAnalysis::getSampleName(myEventProxy);

  if(!myEventProxy.pairs->size()) return true;
  setAnalysisObjects(myEventProxy);

  bool isGoodReco = genLeg1_.getP4().DeltaR(leg1_.getP4())<0.4 &&
    genLeg2_.getP4().DeltaR(leg2_.getP4())<0.4;

  isGoodReco |= genLeg2_.getP4().DeltaR(leg1_.getP4())<0.4 &&
    genLeg1_.getP4().DeltaR(leg2_.getP4())<0.4;

  bool goodGenTau = genLeg1_.getP4().E()>1.0 && genLeg2_.getP4().E()>1.0;							   

  isGoodReco = genLeg2_.getP4().DeltaR(leg2_.getP4())<0.4;
  goodGenTau = genLeg2_.getP4().E()>1.0;

  if(sampleName_=="WAllJets"){
    goodGenTau = true;
    isGoodReco = true;
  }

  if(sampleName_=="DYAllJetsMatchL"){
    goodGenTau = true;
    isGoodReco = true;
  }

  if(isGoodReco && goodGenTau){
    computeMTT("SVfit");
    computeMTT("fastMTT");  
  }


  if(aMessenger and std::string("MLObjectMessenger").compare((aMessenger->name()).substr(0,17))==0 ) // if NULL it will do nothing
  {
    // Putting data to MLObjectMessenger
    try
    {
      MLObjectMessenger* mess = (MLObjectMessenger*)aMessenger;
      mess->putObjectVector(const_cast <const HTTParticle*>(&leg1_), "legs");
      mess->putObjectVector(const_cast <const HTTParticle*>(&leg2_), "legs");
      mess->putObject(const_cast <const HTTParticle*>(&MET_), "met");
      mess->putObject(&covMET_[0][0], "covMET00");
      mess->putObject(&covMET_[0][1], "covMET01");
      mess->putObject(&covMET_[1][0], "covMET10");
      mess->putObject(&covMET_[1][1], "covMET11");
      mess->putObject(&genSumM_, "gen_mass");
      mess->putObject(higgsMassTrans_, "higgs_mass_trans");
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
