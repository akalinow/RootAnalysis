#include <sstream>
#include <bitset>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "Tools.h"

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
  }
  else{
    angles = angleBetweenPlanes(tauLeadingTk, tauPCA, muonTk, muonPCA);
    anglesIPRho = angleBetweenPlanes(tauLeadingTk, aLeg2.getNeutralP4(), muonTk, muonPCA);
  }

  myHistos_->fillProfile("hProfRecoVsMagGen_"+hNameSuffix,
			 aGenLeg2.getPCA().Mag(),
			 tauPCA.Vect().Mag(),
			 eventWeight);

  float yTau =  2.*tauLeadingTk.Pt()/aLeg2.getP4().Pt() - 1.;
  float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;

  if(aLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero && HTTAnalysis::isOneProng(aLeg2.getProperty(PropertyEnum::decayMode))){
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


    myHistos_->fill1DHistogram("h1DPhi-nVectors"+hNameSuffix,angles.first,eventWeight);
    myHistos_->fill2DUnrolledHistogram("h1DUnRollMassSVPhiCP"+hNameSuffix,angles.first, aPair.getP4(aSystEffect).M(),eventWeight);
    myHistos_->fill1DHistogram("h1DCosPhiNN"+hNameSuffix,cosPhiNN);

    float cosMuon =  muonPCA.Vect().Unit()*aGenLeg1.getPCA().Unit();
    float cosTau = tauPCA.Vect().Unit()*aGenLeg2.getPCA().Unit();

    myHistos_->fillProfile("hProfPhiVsMag"+hNameSuffix,aGenLeg1.getPCA().Mag(),cosMuon);
    myHistos_->fillProfile("hProfPhiVsMag"+hNameSuffix,aGenLeg2.getPCA().Mag(),cosTau);

    myHistos_->fillProfile("hProfRecoVsMagGen"+hNameSuffix,aGenLeg1.getPCA().Mag(),muonPCA.Vect().Mag());
    myHistos_->fillProfile("hProfRecoVsMagGen"+hNameSuffix,aGenLeg2.getPCA().Mag(),tauPCA.Vect().Mag());

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
  HTTAnalysis::sysEffects sysType = HTTAnalysis::NOMINAL_SVFIT;

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
