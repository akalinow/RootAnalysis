#include <bitset>

#include "Math/LorentzVector.h"

#include "svfitAnalyzer.h"
#include "MLObjectMessenger.h"
#include "MuTauSpecifics.h"
#include "TauTauSpecifics.h"
#include "MuMuSpecifics.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
svfitAnalyzer::svfitAnalyzer(
	const std::string & aName,
	const std::string & aDecayMode)
 : Analyzer(aName), decayMode(aDecayMode), aCovMET(2,2){

#pragma omp critical
  {
    if(aDecayMode=="MuTau") myChannelSpecifics = new MuTauSpecifics(this);
    else if (aDecayMode=="TauTau") myChannelSpecifics = new TauTauSpecifics(this);
    else if (aDecayMode=="MuMu") myChannelSpecifics = new MuMuSpecifics(this);
    myNumberOfCategories = myChannelSpecifics->getCategoryRejester().size();
    categoryDecisions.resize(myNumberOfCategories);


    svFitAlgo.addLogM_fixed(true, 4.0);
    svFitAlgo.setLikelihoodFileName("");
    svFitAlgo.setMaxObjFunctionCalls(100000);
    svFitAlgo.setVerbosity(1);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
svfitAnalyzer::~svfitAnalyzer(){

  if(myChannelSpecifics) delete myChannelSpecifics;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* svfitAnalyzer::clone() const {

  svfitAnalyzer* clone = new svfitAnalyzer(name(),decayMode);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::initialize(TDirectory* aDir,
			       pat::strbitset *aSelections){
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void svfitAnalyzer::finalize(){

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
	aMET.setP4(met4v);
  //TLorentzVector nunuGen = aGenLeg1.getP4() + aGenLeg2.getP4() -
  //  aGenLeg1.getChargedP4() - aGenLeg2.getChargedP4() -
  //  aGenLeg1.getNeutralP4() - aGenLeg2.getNeutralP4();
	//aMET.setP4(nunuGen);
	
  aCovMET[0][0] = aPair.getMETMatrix().at(0);
  aCovMET[0][1] = aPair.getMETMatrix().at(1);
  aCovMET[1][0] = aPair.getMETMatrix().at(2);
  aCovMET[1][1] = aPair.getMETMatrix().at(3);

  myChannelSpecifics->setAnalysisObjects(myEventProxy);

  aSeparatedJets = getSeparatedJets(myEventProxy, 0.5);
  aJet1 = aSeparatedJets.size() ? aSeparatedJets[0] : HTTParticle();
  aJet2 = aSeparatedJets.size()>1 ? aSeparatedJets[1] : HTTParticle();
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
  if(aCovMET[0][0]==0 && aCovMET[1][0]==0 && aCovMET[0][1]==0 && aCovMET[1][1]==0) return TLorentzVector(); //singular covariance matrix

  TLorentzVector aResult;
  if(algoName=="svfit") aResult = runSVFitAlgo(measuredTauLeptons);
  if(algoName=="fastMTT") aResult = runFastMTTAlgo(measuredTauLeptons);
  
  return aResult;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::runFastMTTAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons){

  fastMTTAlgo.run(measuredTauLeptons, aMET.getP4().X(), aMET.getP4().Y(), aCovMET);
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
TLorentzVector svfitAnalyzer::runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons){

  svFitAlgo.setVerbosity(0);

  svFitAlgo.addLogM_fixed(true, 4.0);
  if     (decayMode=="MuTau" ) svFitAlgo.addLogM_fixed(true, 4.0);
  else if(decayMode=="TauTau") svFitAlgo.addLogM_fixed(true, 5.0);
  else if(decayMode=="MuMu"  ) svFitAlgo.addLogM_fixed(true, 3.0);
  
  svFitAlgo.setMaxObjFunctionCalls(100000);
  //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  //svFitAlgo.setTreeFileName("markovChainTree.root");
  TVector3 recPV = aEvent.getRefittedPV();
  TVector3 svLeg2 = aLeg2.getSV();
  svLeg2 -=recPV;

  TVector3 svLeg1 = aGenLeg1.getSV();
  svLeg1 -=recPV;

  svFitAlgo.integrate(measuredTauLeptons, aMET.getP4().X(), aMET.getP4().Y(), aCovMET);
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
    computeMTT("svfit");
		computeMTT("fastMTT");  
  }

	aGenSumM = (aGenLeg1.getP4() + aGenLeg2.getP4()).M();
	higgs_mass_trans = aPair.getMTMuon();

	if(aMessenger and std::string("MLObjectMessenger").compare((aMessenger->name()).substr(0,17))==0 ) // if NULL it will do nothing
	{
		// Putting data to MLObjectMessenger
		try
		{
			MLObjectMessenger* mess = (MLObjectMessenger*)aMessenger;
			mess->putObjectVector(const_cast <const HTTParticle*>(&aLeg1), "legs");
			mess->putObjectVector(const_cast <const HTTParticle*>(&aLeg2), "legs");
			mess->putObject(const_cast <const HTTParticle*>(&aMET), "amet");
			mess->putObject(&aCovMET[0][0], "covMET00");
			mess->putObject(&aCovMET[0][1], "covMET01");
			mess->putObject(&aCovMET[1][0], "covMET10");
			mess->putObject(&aCovMET[1][1], "covMET11");
			mess->putObject(&aGenSumM, "gen_visible_mass");
			mess->putObject(higgs_mass_trans, "higgs_mass_trans");
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
