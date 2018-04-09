#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "svfitAnalyzer.h"
#include "svfitHistograms.h"
#include "MuTauSpecifics.h"
#include "TauTauSpecifics.h"
#include "Tools.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
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

          //TEST svFitAlgo.addLogM_fixed(true, 4.0);
          svFitAlgo.addLogM_fixed(false, 4.0);
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
TLorentzVector svfitAnalyzer::runFastSVFitAlgo(const TLorentzVector & leg1P4,
                                               const TLorentzVector & leg2P4,
                                               const TLorentzVector &metP4,
                                               const TMatrixD &covMET){

  double x1 = 0.0;
  double x2 = 0.0;
  TLorentzVector tau1P4, tau2P4;
  TLorentzVector nuP4;
  TLorentzVector met1;
  TLorentzVector p4SVFit;
  double maxLLH = 0;
  double llh = 0.0;
  double mH = 0.0;
  double mVis = (leg1P4 + leg2P4).M();
  double mVisLeg2 = leg2P4.M();
  double x2Min = 0.0;

    for(int iX2 = 1; iX2<100;++iX2){

      x2 = iX2/100.0;
      if(mVisLeg2>1.5) mVisLeg2 = 1.5;
      x2Min = std::pow(mVisLeg2/1.77685,2);
      if(x2<x2Min) continue;
      tau2P4 = leg2P4*(1.0/x2);

      for(int iX1 = 1; iX1<100;++iX1){

      x1 = iX1/100.0;

    tau1P4 = leg1P4*(1.0/x1);

    nuP4 = tau1P4 - aLeg1.getP4();
    nuP4 += tau2P4 - aLeg2.getP4();

    mH = (tau1P4+tau2P4).M();
    x2Min = std::max(std::pow(mVisLeg2/1.77685,2), std::pow(mVis/mH,2));
    if(x2<x2Min) continue;

    llh = EvalMET_TF(metP4, nuP4, covMET);
    //llh *= std::pow(mVis,2)*std::pow(mH,-3)*(2*log(mH/mVis) + std::pow(mVis/mH,2)*(1 - std::pow(mH/mVis,2)));

    llh *= std::pow(mVis,2)*std::pow(mH,-3)*(-log(x2Min) + std::pow(mVis/mH,2)*(1 - std::pow(x2Min,-1)));
    //std::cout<<"llh: "<<llh<<std::endl;
    if(llh>maxLLH){
      maxLLH = llh;
      p4SVFit = tau1P4 + tau2P4;
    }
  }
  }
  if(p4SVFit.M()>20000){
    std::cout<<"fast Mass: "<<p4SVFit.M()
             <<" classic mass: "<<computeSvFit().M()
             <<std::endl;
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
    svFitAlgo.setMaxObjFunctionCalls(10000);
    //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
    //svFitAlgo.setTreeFileName("markovChainTree.root");
    TVector3 recPV = aEvent.getRefittedPV();
    TVector3 svLeg2 = aLeg2.getSV();
    svLeg2 -=recPV;

    TVector3 svLeg1 = aGenLeg1.getSV();
    svLeg1 -=recPV;

    svFitAlgo.addSVData(svLeg1, svLeg2);

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

        myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.getP4().M());
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

        TLorentzVector svFitP4 = runFastSVFitAlgo(aLeg1.getP4(), aLeg2.getP4(), aMET.getP4(), covMET);
        myHistos_->fill1DHistogram("h1DMassSVRecalculated"+hNameSuffix,svFitP4.M());
        ////


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

        bool isGoodReco = aGenLeg2.getP4().DeltaR(aLeg2.getP4())<0.4;
        bool goodGenTau = aGenLeg2.getP4().E()>1.0;

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
        if(!isOneProng && goodGenTau && isGoodReco){

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

        fillControlHistos(hNameSuffix);

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
