#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "TF1.h"

#include "svfitAnalyzer.h"
#include "svfitHistograms.h"
#include "MuTauSpecifics.h"
#include "TauTauSpecifics.h"
#include "Tools.h"

#include "FastMTT.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Double_t likelihoodFunc(Double_t *x, Double_t *par)
{
  Float_t mH =x[0];
  Double_t mVis = par[0];
  Double_t mVisLeg1 = par[1];
  Double_t mVisLeg2 = par[2];
  Double_t coeff1 = par[3];
  Double_t scale = par[4];

  Double_t mTau = 1.77685;

  Double_t x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  Double_t x2Min = std::max(std::pow(mVisLeg2/mTau,2), std::pow(mVis/mH,2));
  Double_t x2Max = std::min(1.0, std::pow(mVis/mH,2)/x1Min);

  Double_t value = 2.0*std::pow(mVis,2)*std::pow(mH,-coeff1)*(log(x2Max)-log(x2Min) + std::pow(mVis/mH,2)*(1 - std::pow(x2Min,-1)));
  if(mH<mVis) return 0.0;
  return scale*value;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
svfitAnalyzer::svfitAnalyzer(const std::string & aName, const std::string & aDecayMode) : Analyzer(aName){

#pragma omp critical
  {
    if(aDecayMode=="MuTau") myChannelSpecifics = new MuTauSpecifics(this);
    else if (aDecayMode=="TauTau") myChannelSpecifics = new TauTauSpecifics(this);
    myNumberOfCategories = myChannelSpecifics->getCategoryRejester().size();
    categoryDecisions.resize(myNumberOfCategories);

    ntupleFile_ = 0;

    svFitAlgo.addLogM_fixed(true, 4.0);
    svFitAlgo.setLikelihoodFileName("");
    svFitAlgo.setMaxObjFunctionCalls(100000);
    svFitAlgo.setVerbosity(1);
  }

  fLikelihood = new TF1("likelihood",likelihoodFunc,0,5000,5);
  fLikelihood->SetParNames("mVis","mVisLeg1","mVisLeg2", "coefficient","scale");
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
TLorentzVector svfitAnalyzer::runFastSVFitAlgo(const TLorentzVector & leg1P4,
                                               const TLorentzVector & leg2P4,
                                               const TLorentzVector &metP4,
                                               const TMatrixD &covMET){

  double x1 = 0.0;
  double x2 = 0.0;

  double x1Max = 0.0;
  double x2Max = 0.0;

  TLorentzVector tau1P4, tau2P4;
  TLorentzVector nuP4;
  TLorentzVector met1;
  TLorentzVector p4SVFit;
  double maxLLH = 0;
  double llh = 0.0;
  double mH = 0.0;
  double mVis = (leg1P4 + leg2P4).M();
  double mVisLeg1 = leg1P4.M();
  double mVisLeg2 = leg2P4.M();
  if(mVisLeg2>1.5) mVisLeg2 = 0.6;
  mVisLeg2 = 0.3;
  double x2Min = 0.0;

  TLorentzVector nuP4Max;
  TLorentzVector nunuGen = aGenLeg1.getP4() + aGenLeg2.getP4() -
    leg1P4  - leg2P4;

  double x1True =  leg1P4.E()/aGenLeg1.getP4().E();
  double x2True =  leg2P4.E()/aGenLeg2.getP4().E();

  fLikelihood->FixParameter(0,mVis);
  fLikelihood->FixParameter(1,mVisLeg1);
  fLikelihood->FixParameter(2,mVisLeg2);
  fLikelihood->FixParameter(3,6);
  fLikelihood->FixParameter(4,1.0);

  int nGridPoints = 100;
  int nCalls = 0;
  for(int iX2 = 1; iX2<nGridPoints;++iX2){
    //for(int iX2 = 1; iX2<2;++iX2){

    x2 = (double)iX2/nGridPoints;
    //x2 = x2True;

    x2Min = std::pow(mVisLeg2/1.77685,2);
    if(x2<-x2Min) continue; //TEST "-"

    tau2P4 = leg2P4*(1.0/x2);
    tau2P4.SetVectM(tau2P4.Vect(),1.77685);

    for(int iX1 = 1; iX1<nGridPoints;++iX1){
      //for(int iX1 = 1; iX1<2;++iX1){

      x1 = (double)iX1/nGridPoints;

      //x1 = std::pow(mVis,2)/std::pow(123.0,2)/x2;
      //x1 = x1True;

      tau1P4 = leg1P4*(1.0/x1);
      tau1P4.SetVectM(tau1P4.Vect(),1.77685);

      TLorentzVector nu1P4 = tau1P4 - leg1P4;
      TLorentzVector nu2P4 = tau2P4 - leg2P4;
      nuP4 = nu1P4 + nu2P4;

      mH = (tau1P4+tau2P4).M();
      x2Min = std::max(std::pow(mVisLeg2/1.77685,2), std::pow(mVis/mH,2));
      //if(x2<x2Min) continue; TEST

      llh = EvalMET_TF(metP4, nuP4, covMET);
      //llh *= fLikelihood->Eval(mH/1.17);
      ++nCalls;

      if(llh>maxLLH){
	maxLLH = llh;
	p4SVFit = (tau1P4 + tau2P4);
	nuP4Max = nuP4;
	x1Max = x1;
	x2Max = x2;
      }
    }
  }
  /*
    std::cout<<" nCalls: "<<nCalls
    <<" x1Max: "<<x1Max
    <<" x2Max: "<<x2Max
    <<" max LLH: "<<maxLLH
    <<std::endl;
  */

  myHistos_->fill1DHistogram("h1DLLH_1",x1Max/x1True);
  myHistos_->fill1DHistogram("h1DLLH_2",x2Max/x2True);

  double delta = leg2P4.E() - aGenLeg2.getChargedP4().E();
  delta /= aGenLeg2.getChargedP4().E();
  myHistos_->fill1DHistogram("h1DLLH_3",delta);

  delta = leg1P4.E() - aGenLeg1.getChargedP4().E();
  delta /= aGenLeg1.getChargedP4().E();
  myHistos_->fill1DHistogram("h1DLLH_4",delta);


  myHistos_->fill2DHistogram("h2DDelta_1",
			     (metP4.Vect().DeltaPhi(nunuGen.Vect())),
			     (metP4.Vect().Mag() - nunuGen.Vect().Mag())/nunuGen.Vect().Mag());

  myHistos_->fill2DHistogram("h2DDelta_2",
			     (x1Max - x1True)/x1True,
			     (x2Max - x2True)/x2True);


  if(delta>0.05){

    double deltaR1 = nunuGen.DeltaPhi(metP4);
    deltaR1 = nunuGen.DeltaPhi(nuP4Max);

    //deltaR1 = nunuGen.E()/nuP4Max.E();
    deltaR1 = nunuGen.E()/metP4.E();

    myHistos_->fill1DHistogram("h1DDeltaR_1",deltaR1);
  }

  ////
  return p4SVFit;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double svfitAnalyzer::EvalMET_TF(const TLorentzVector & metP4,
                                 const TLorentzVector & nuP4,
                                 const TMatrixD& covMET){

  double  aMETx = metP4.X();
  double  aMETy = metP4.Y();
  // determine transfer matrix for MET
  double invCovMETxx = covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy = covMET(0,0);
  double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet)<1E-10){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
	      <<"METx: "<<aMETy<<" METy: "<<aMETy
	      << std::endl;
    //errorCode_ |= MatrixInversion; //FIXME violates const
    return 0;
  }
  double const_MET = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));

  // evaluate transfer function for MET/hadronic recoil
  double residualX = aMETx - (nuP4.X());
  double residualY = aMETy - (nuP4.Y());

  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
    residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2/=covDet;

  return const_MET*TMath::Exp(-0.5*pull2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::computeSvFit(){

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
  //decay1 = aLeg2.getProperty(PropertyEnum::decayMode);
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

  return runSVFitAlgo(measuredTauLeptons, aMET, covMET);

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector svfitAnalyzer::runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
                                           const TVector2 &aMET, const TMatrixD &covMET){

  svFitAlgo.setVerbosity(0);
  svFitAlgo.addLogM_fixed(true, 4.0);
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

  //std::cout<<"SVFit mass: "<<mcMass<<std::endl;

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
    std::cout<<"mVisLeg2: "<<aLeg2.getP4().M()<<std::endl;
    TLorentzVector test = computeSvFit();
    std::cout<<"Mass: "<<test.M()<<std::endl;
  */

  std::vector<double> shapeParams = {6, 1/1.17};
  FastMTT testMinimizer;
  testMinimizer.setLikelihoodParams(shapeParams);
  testMinimizer.scan(aLeg1.getP4(), aLeg2.getP4(), aMET.getP4(), covMET, 1, 1);
  //testMinimizer.scan(aLeg1.getP4(), aLeg2.getP4(), aMET.getP4(), covMET, 1, 1);
  //testMinimizer.scan(aGenLeg1.getChargedP4(), aGenLeg2.getChargedP4(), nunuGen, covMET, 1, 1);
  //testMinimizer.scan(aLeg1.getP4(), aLeg2.getP4(), nunuGen, covMET, 1, 1);
  TLorentzVector svFitP4 = testMinimizer.getBestP4();
  //TLorentzVector svFitP4 = runFastSVFitAlgo(aLeg1.getP4(), aLeg2.getP4(), aMET.getP4(), covMET);
  //TLorentzVector svFitP4 = runFastSVFitAlgo(aGenLeg1.getChargedP4(), aLeg2.getP4(), nunuGen, covMET);
  //TLorentzVector svFitP4 = runFastSVFitAlgo(aLeg1.getP4(), aLeg2.getP4(), aMET.getP4(), covMET);
  myHistos_->fill1DHistogram("h1DMassSVFast"+hNameSuffix,svFitP4.M());
  myHistos_->fill1DHistogram("h1DCpuTimeFast"+hNameSuffix,testMinimizer.getCpuTime("scan"));

  double delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaFast"+hNameSuffix,delta);

  delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiFast"+hNameSuffix,delta);

  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtFast"+hNameSuffix,delta);
  ////
  /*
    svFitP4 = computeSvFit();
    myHistos_->fill1DHistogram("h1DMassSVClassic"+hNameSuffix,svFitP4.M());
    myHistos_->fill1DHistogram("h1DCpuTimeClassic"+hNameSuffix,svFitAlgo.getComputingTime_cpu());
    delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
    myHistos_->fill1DHistogram("h1DDeltaPtClassic"+hNameSuffix,delta);
  */

  svFitP4 = aPair.getP4();
  delta = svFitP4.Eta() - tautauGen.Eta();
  myHistos_->fill1DHistogram("h1DDeltaEtaStandalone"+hNameSuffix,delta);

  delta = svFitP4.Vect().DeltaPhi(tautauGen.Vect());
  myHistos_->fill1DHistogram("h1DDeltaPhiStandalone"+hNameSuffix,delta);

  delta = (svFitP4.Perp() - tautauGen.Perp())/tautauGen.Perp();
  myHistos_->fill1DHistogram("h1DDeltaPtStandalone"+hNameSuffix,delta);
  /////


  double sinThetaGenReco = (aGenLeg1.getSV() - genPV).Unit()*aLeg1.getChargedP4().Vect().Unit();
  sinThetaGenReco = sqrt(1.0 - std::pow(sinThetaGenReco, 2));
  double flightPathRec = aLeg1.getPCARefitPV().Mag()/sinThetaGenReco;
  flightPathRec /= aGenLeg1.getP4().Gamma();
  flightPathRec /= aGenLeg1.getP4().Beta();

  if(aGenLeg1.getP4().DeltaR(aLeg1.getP4())<0.4 && aGenLeg1.getP4().E()>1.0){
    myHistos_->fill1DHistogram("h1DFlightPathPCARecLeg1"+hNameSuffix,flightPathRec);
  }

  TVector3 genSV = aGenLeg2.getSV();
  TVector3 genPCA = aGenLeg2.getPCA();

  TVector3 recSV = aLeg2.getSV();
  TVector3 recPCA = aLeg2.getPCARefitPV();

  //std::cout<<hNameSuffix<<std::endl;
  bool isOneProng = HTTAnalysis::isOneProng(aLeg2.getProperty(PropertyEnum::decayMode));
  //bool isSeparated = aGenLeg2.getP4().DeltaR(aGenLeg2.getP4())>0.4;
  //bool isSeparated = aGenLeg2.getP4().DeltaR(aLeg2.getP4())>0.4;

  double sinThetaReco = (aGenLeg2.getSV() - genPV).Unit()*aLeg2.getChargedP4().Vect().Unit();
  sinThetaReco = sqrt(1.0 - std::pow(sinThetaReco, 2));
  double sinThetaGen = (genSV - genPV).Unit()*aGenLeg2.getChargedP4().Vect().Unit();
  sinThetaGen = sqrt(1.0 - std::pow(sinThetaGen, 2));

  //isGoodReco &= recPCA.Mag()>0.05;
  //isGoodReco &= std::abs(recPCA.Mag() - genPCA.Mag())/genPCA.Mag() < 0.1;
  //isGoodReco &= std::abs(sinThetaReco - sinThetaGen)/sinThetaGen<0.1;

  //svFitP4 = computeSvFit();
  //myHistos_->fill1DHistogram("h1DMassSVRecalculated"+hNameSuffix,svFitP4.M());
  //isOneProng = false;
  if(!isOneProng){

    flightPathRec = (recSV - recPV).Mag();
    flightPathRec /= aGenLeg2.getP4().Gamma();
    flightPathRec /= aGenLeg2.getP4().Beta();
    //std::cout<<" flightPathRec: "<<flightPathRec<<std::endl;
    myHistos_->fill1DHistogram("h1DFlightPathRec"+hNameSuffix,flightPathRec);

    flightPathRec = recPCA.Mag()/sinThetaReco;
    flightPathRec /= aGenLeg2.getP4().Gamma();
    flightPathRec /= aGenLeg2.getP4().Beta();
    myHistos_->fill1DHistogram("h1DFlightPathPCARec"+hNameSuffix,flightPathRec);
    //std::cout<<" flightPathRec PCA: "<<flightPathRec;

    double flightPathGen = (genSV - genPV).Mag();
    flightPathGen /= aGenLeg2.getP4().Gamma();
    flightPathGen /= aGenLeg2.getP4().Beta();
    myHistos_->fill1DHistogram("h1DFlightPathGen"+hNameSuffix,flightPathGen);
    //std::cout<<" flightPathGen PCA: "<<flightPathGen<<std::endl;

    double deltaR = (recPCA.Mag() - genPCA.Mag())/genPCA.Mag();
    myHistos_->fill2DHistogram("h2DFlightPathVsDeltaRGen"+hNameSuffix, deltaR, flightPathGen);

    double flightPathGenPCA = genPCA.Mag()/sinThetaGen;
    flightPathGenPCA /= aGenLeg2.getP4().Gamma();
    flightPathGenPCA /= aGenLeg2.getP4().Beta();
    myHistos_->fill1DHistogram("h1DFlightPathPCAGen"+hNameSuffix,flightPathGenPCA);
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

  bool isGoodReco = aGenLeg2.getP4().DeltaR(aLeg2.getP4())<0.4;
  bool goodGenTau = aGenLeg2.getP4().E()>1.0;
  double delta = aLeg2.getP4().E() - aGenLeg2.getChargedP4().E();
  delta /= aGenLeg2.getChargedP4().E();
  //isGoodReco &= std::abs(delta)<0.1;

  if(isGoodReco && goodGenTau){
    fillControlHistos(hNameSuffix);
  }

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
