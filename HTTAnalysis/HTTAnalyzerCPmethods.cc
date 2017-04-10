#include <sstream>
#include <bitset>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "Tools.h"

#include "TF1.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight,
                                      const HTTAnalysis::sysEffects & aSystEffect){

  ///Method from http://arxiv.org/abs/1108.0670 (S. Berge)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles;

  ///Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  ///Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  std::pair<float,float>  anglesIPRho;

  const TLorentzVector & muonTk = aLeg1.getChargedP4();
  const TLorentzVector & tauLeadingTk = aLeg2.getChargedP4();

  TLorentzVector muonPCA(aLeg1.getPCA(),0);
  TLorentzVector tauPCA(aLeg2.getPCA(),0);

  if(hNameSuffix.find("AODPV")!=std::string::npos){
    muonPCA.SetVect(aLeg1.getPCA());
    tauPCA.SetVect(aLeg2.getPCA());
  }
  if(hNameSuffix.find("RefitPV")!=std::string::npos){
    muonPCA.SetVect(aLeg1.getPCARefitPV());
    tauPCA.SetVect(aLeg2.getPCARefitPV());
  }
  if(hNameSuffix.find("GenPV")!=std::string::npos){
    muonPCA.SetVect(aLeg1.getPCAGenPV());
    tauPCA.SetVect(aLeg2.getPCAGenPV());
  }

  if(aLeg1.getCharge()>0) {
    angles = angleBetweenPlanes(muonTk, muonPCA, tauLeadingTk, tauPCA);
    anglesIPRho = angleBetweenPlanes(muonTk, muonPCA, tauLeadingTk, aLeg2.getNeutralP4());
    //angles = angleBetweenPlanes(muonTk, muonPCA, aGenLeg2.getP4(), aGenLeg2.getChargedP4());//TEST
  }
  else{
    angles = angleBetweenPlanes(tauLeadingTk, tauPCA, muonTk, muonPCA);
    anglesIPRho = angleBetweenPlanes(tauLeadingTk, aLeg2.getNeutralP4(), muonTk, muonPCA);
    //angles = angleBetweenPlanes(aGenLeg2.getP4(), aGenLeg2.getChargedP4(), muonTk, muonPCA);//TEST
  }

  float yTau =  2.*tauLeadingTk.Pt()/aLeg2.getP4().Pt() - 1.;
  float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;

  if(aLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero &&
    HTTAnalysis::isOneProng(aLeg2.getProperty(PropertyEnum::decayMode))){
     myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);

     myHistos_->fill1DHistogram("h1DPhi-nVecIP"+hNameSuffix,shiftedIPrho,eventWeight);
     myHistos_->fill2DUnrolledHistogram("h1DUnRollMassSVYCP"+hNameSuffix, shiftedIPrho, aPair.getP4(aSystEffect).M(),eventWeight);

     if(yTau>0){
       myHistos_->fill1DHistogram("h1DPhi-nVecIP-yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
     else{
       myHistos_->fill1DHistogram("h1DPhi-nVecIP-yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
  }

  if(aLeg2.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero){
    float cosPhiNN =  tauPCA.Vect().Unit().Dot(muonPCA.Vect().Unit());


    float x =  muonPCA.Vect().Mag()/0.004;
    float eventWeightTmp = eventWeight;
    //eventWeightTmp /= exp(-0.656064-0.343095*x);
    //eventWeightTmp /= exp(-0.656064-0.343095*x);
    //if(x>5) eventWeightTmp = 0;

    myHistos_->fill1DHistogram("h1DPhi-nVectors"+hNameSuffix,angles.first,eventWeightTmp);
    myHistos_->fill2DUnrolledHistogram("h1DUnRollMassSVPhiCP"+hNameSuffix,angles.first, aPair.getP4(aSystEffect).M(),eventWeight);
    myHistos_->fill1DHistogram("h1DCosPhiNN"+hNameSuffix,cosPhiNN);

    float cosMuon =  muonPCA.Vect().Unit()*aGenLeg1.getPCA().Unit();
    float cosTau = tauPCA.Vect().Unit()*aGenLeg2.getPCA().Unit();

    myHistos_->fillProfile("hProfPhiVsMag"+hNameSuffix,aGenLeg1.getPCA().Mag(),cosMuon);
    myHistos_->fillProfile("hProfPhiVsMag"+hNameSuffix,aGenLeg2.getPCA().Mag(),cosTau);

    myHistos_->fillProfile("hProfRecoVsMagGen"+hNameSuffix,aGenLeg1.getPCA().Mag(),muonPCA.Vect().Mag());
    myHistos_->fillProfile("hProfRecoVsMagGen"+hNameSuffix,aGenLeg2.getPCA().Mag(),tauPCA.Vect().Mag());

//////////////////TEST
    float deltaPhiLegs = std::abs(aLeg1.getP4().Vect().DeltaPhi(aLeg2.getP4().Vect()))/M_PI;
    float deltaPhiMET = std::abs(aLeg1.getP4().Vect().DeltaPhi(aMET.getP4().Vect()))/M_PI;

    float dxyLeg1 = aLeg1.getProperty(PropertyEnum::dxy);
    float dxyLeg2 = aLeg2.getProperty(PropertyEnum::dxy);

    float dzLeg1 = aLeg1.getProperty(PropertyEnum::dz);
    float dzLeg2 = aLeg2.getProperty(PropertyEnum::dz);

    float delta_dxy = std::abs(dxyLeg1 + dxyLeg2);
    float delta_dz = std::abs(dzLeg1 + dzLeg2);
    float abseta = std::abs(aLeg1.getP4().Eta());

    TLorentzVector leg1SVFitP4 = aPair.getLeg1P4();
    TLorentzVector leg2SVFitP4 = aPair.getLeg2P4();


    TLorentzVector getTauPCA(aGenLeg2.getPCA(),0);
    angles = angleBetweenPlanes(muonTk, muonPCA, aGenLeg2.getChargedP4(), getTauPCA);
    std::pair<float,float> anglesSVFit1 = angleBetweenPlanes(muonTk, leg1SVFitP4, aGenLeg2.getChargedP4(), aGenLeg2.getP4());
    std::pair<float,float> anglesSVFit = angleBetweenPlanes(leg1SVFitP4, muonTk, aGenLeg2.getP4(), aGenLeg2.getChargedP4());
    std::pair<float,float> anglesGen = angleBetweenPlanes(aGenLeg1.getP4(), aGenLeg1.getChargedP4(), aGenLeg2.getP4(), aGenLeg2.getChargedP4());

    float delta1 = leg1SVFitP4.Vect().Angle(aGenLeg1.getP4().Vect())/M_PI;
    //delta1 = leg1SVFitP4.DeltaR(aGenLeg1.getP4());
    //delta1 = leg1SVFitP4.E()/aGenLeg1.getP4().E();
    //delta1 = (leg1SVFitP4.Vect().X()-aGenLeg1.getP4().Vect().X())/aGenLeg1.getP4().Vect().X();

    float delta2 = leg2SVFitP4.Vect().Angle(aGenLeg2.getP4().Vect())/M_PI;
    //delta2 = leg2SVFitP4.DeltaR(aGenLeg2.getP4());
    //delta2 = leg2SVFitP4.E()/aGenLeg2.getP4().E();
    //delta2 = (leg1SVFitP4.Vect().Z()-aGenLeg1.getP4().Vect().Z())/aGenLeg1.getP4().Vect().Z();

    float delta3 = muonTk.Vect().Angle(aGenLeg1.getChargedP4().Vect());
    float delta4 = leg1SVFitP4.Vect().Unit().Angle(muonPCA.Vect().Unit());
    //delta4 = aGenLeg1.getP4().Vect().Unit().Angle(muonPCA.Vect().Unit());
    delta4-=M_PI/2;
    delta4*=-1;
    delta4/=M_PI;

    float delta5 = muonPCA.Vect().Angle(aGenLeg1.getPCA())/M_PI;

    float delta6 = leg1SVFitP4.Vect().Angle(muonTk.Vect())/M_PI;
    float delta7 = aGenLeg1.getP4().Vect().Angle(aGenLeg1.getChargedP4().Vect())/M_PI;
    float delta8 = muonPCA.Vect().Mag()/0.004;

    float delta9 = aGenLeg1.getP4().Vect().Angle(aGenLeg1.getPCA());
    delta9-=M_PI/2;
    delta9*=-1;
    delta9/=M_PI;

    float delta10 = aGenLeg1.getPCA().Mag();
    float delta11 = aGenLeg1.getPCA().Angle(muonPCA.Vect())/M_PI;
          //delta11 = muonPCA.Vect().Mag()/aGenLeg1.getPCA().Mag();
    float delta12 = (aGenLeg1.getP4() + aGenLeg2.getP4()).Perp();
          //delta12 = (aLeg1.getP4() + aLeg2.getP4()).Phi();
          delta12 = aGenLeg1.getChargedP4().Vect().Eta() - aGenLeg1.getP4().Vect().Eta();
          //delta12 = aLeg1.getChargedP4().Vect().Eta() - leg1SVFitP4.Vect().Eta();

    float delta13 = anglesSVFit.first - anglesSVFit1.first;
    float delta14 = aEvent.getGenPV().X() - aEvent.getRefittedPV().X();


    float gamma = leg1SVFitP4.Gamma();
    //gamma = aGenLeg1.getP4().Gamma();

    if(delta7>0.5){
    std::cout << "hNameSuffix: "<<hNameSuffix
    <<" delta4: "<<delta4 <<" delta7: "<<delta7
    <<" delta4/delta7: "<<delta4/delta7
    <<" delta9: "<<delta9
    <<" Tau pt: "<<aGenLeg1.getP4().Perp()<<" mupn pt: "<<muonTk.Perp()
    <<" Higgs pt: "<<(aGenLeg1.getP4() + aGenLeg2.getP4()).Perp()
              << '\n';
  }

    //delta4*=gamma;

    float deltaSVFit = std::abs(anglesSVFit.first - anglesGen.first);
    float deltaNPCA = angles.first - anglesGen.first;
    float delta = std::min(deltaSVFit, deltaNPCA);

    if(delta>M_PI) delta-=2*M_PI;
    delta = std::abs(delta/M_PI);

    if(deltaNPCA>M_PI) deltaNPCA-=2*M_PI;
    else if (deltaNPCA<-M_PI) deltaNPCA+=2*M_PI;
    deltaNPCA = deltaNPCA/M_PI;

    if(deltaSVFit>M_PI) deltaSVFit-=2*M_PI;
    deltaSVFit = std::abs(deltaSVFit/M_PI);

    eventWeight = 1.0;
    eventWeight /= exp(-0.656064-0.343095*delta8);
    eventWeight /= exp(-0.656064-0.343095*delta8);
    if(delta8<5) myHistos_->fill2DHistogram("h2DTestHisto"+hNameSuffix,delta8,angles.first, eventWeight);
    myHistos_->fillProfile("hProfTest1"+hNameSuffix,delta8,deltaNPCA);

//if(std::abs(delta4)<10.005) myHistos_->fill2DHistogram("h2DTestHisto"+hNameSuffix,delta11,deltaNPCA);
///////////////////////////
    }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight){

  ///Method from http://arxiv.org/abs/1108.0670 (Berger)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles;

  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  std::pair<float,float>  anglesIPRho;

  TLorentzVector muonTk = aGenLeg1.getChargedP4();
  TLorentzVector tauLeadingTk = aGenLeg2.getChargedP4();
  TLorentzVector muonPCA(aGenLeg1.getPCA(),0);
  TLorentzVector tauPCA(aGenLeg2.getPCA(),0);

  if(aGenLeg1.getCharge()>0) {
    angles = angleBetweenPlanes(muonTk, muonPCA, tauLeadingTk, tauPCA);
    anglesIPRho = angleBetweenPlanes(muonTk, muonPCA, tauLeadingTk, aGenLeg2.getNeutralP4());
  }
  else{
    angles = angleBetweenPlanes(tauLeadingTk, tauPCA, muonTk, muonPCA);
    anglesIPRho = angleBetweenPlanes(tauLeadingTk, aGenLeg2.getNeutralP4(), muonTk, muonPCA);
  }

  ///////////////////////////////////////////////////////////
  HTTAnalysis::sysEffects sysType = HTTAnalysis::NOMINAL;

  if(aGenLeg2.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero){
    float cosPhiNN =  muonPCA.Vect().Unit().Dot(tauPCA.Vect().Unit());
    myHistos_->fill1DHistogram("h1DCosPhiNN"+hNameSuffix,cosPhiNN);
    myHistos_->fill1DHistogram("h1DPhi-nVectors"+hNameSuffix,angles.first,eventWeight);
    myHistos_->fill2DUnrolledHistogram("h1DUnRollMassSVPhiCP"+hNameSuffix,angles.first, aPair.getP4(sysType).M(),eventWeight);
  }
  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  if(aGenLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero && HTTAnalysis::isOneProng(aGenLeg2.getProperty(PropertyEnum::decayMode))){

    float yTau =  2.*tauLeadingTk.Pt()/(aGenLeg2.getChargedP4()+aGenLeg2.getNeutralP4()).Pt() - 1.;
    float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;

    myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);
    myHistos_->fill1DHistogram("h1DPhi-nVecIP"+hNameSuffix,shiftedIPrho,eventWeight);
    myHistos_->fill2DUnrolledHistogram("h1DUnRollMassSVYCP"+hNameSuffix, shiftedIPrho, aPair.getP4(sysType).M(),eventWeight);

    if(yTau>0) myHistos_->fill1DHistogram("h1DPhi-nVecIP-yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
    else myHistos_->fill1DHistogram("h1DPhi-nVecIP-yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
  }
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
  //TVector3 boost = aPair.getP4().BoostVector();
  //TVector3 boost = (aGenLeg1.getP4()+aGenLeg2.getP4()).BoostVector();

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
