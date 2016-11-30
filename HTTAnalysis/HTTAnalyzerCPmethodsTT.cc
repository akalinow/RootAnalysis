#include <sstream>
#include <bitset>

#include "HTTAnalyzerTT.h"
#include "HTTHistogramsTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight){

  ///Method from http://arxiv.org/abs/1108.0670 (S. Berge)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles;

  ///Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  ///Angle between rho decay planes in the rho-rho.
  std::pair<float,float>  anglesRhoRho;  
  ///Here we take rho on one side, pi+/- on the other
  std::pair<float,float>  anglesIPRho;  

  TLorentzVector tau1LeadingTk = aTau1.getChargedP4();
  TLorentzVector tau2LeadingTk = aTau2.getChargedP4();
  TLorentzVector tau1PCA(aTau1.getPCA(),0);
  TLorentzVector tau2PCA(aTau2.getPCA(),0);

  if(hNameSuffix.find("AODPV")!=std::string::npos){
    tau1PCA = TLorentzVector(aTau1.getPCA(),0);
    tau2PCA = TLorentzVector(aTau2.getPCA(),0);
  }    
  if(hNameSuffix.find("RefitPV")!=std::string::npos){
    tau1PCA = TLorentzVector(aTau1.getPCARefitPV(),0);
    tau2PCA = TLorentzVector(aTau2.getPCARefitPV(),0);
  }
  if(hNameSuffix.find("GenPV")!=std::string::npos){
    tau1PCA = TLorentzVector(aTau1.getPCAGenPV(),0);
    tau2PCA = TLorentzVector(aTau2.getPCAGenPV(),0);
  }

  if(aTau1.getCharge()>0) {
    angles = angleBetweenPlanes(tau1LeadingTk, tau1PCA, tau2LeadingTk, tau2PCA);
    anglesIPRho = angleBetweenPlanes(tau1LeadingTk, tau1PCA, tau2LeadingTk, aTau2.getNeutralP4());
    anglesRhoRho = angleBetweenPlanes(tau1LeadingTk, aTau1.getNeutralP4(), tau2LeadingTk, aTau2.getNeutralP4());
  }
  else{
    angles = angleBetweenPlanes(tau2LeadingTk, tau2PCA, tau1LeadingTk, tau1PCA);
    anglesIPRho = angleBetweenPlanes(tau2LeadingTk, tau2PCA, tau1LeadingTk, aTau1.getNeutralP4());
    anglesRhoRho = angleBetweenPlanes(tau2LeadingTk, aTau2.getNeutralP4(), tau1LeadingTk, aTau1.getNeutralP4());
  }

  myHistos_->fillProfile("hProfRecoVsMagGen_"+hNameSuffix,
			 aGenTau1.getPCA().Mag(),
			 tau1PCA.Vect().Mag(),
			 eventWeight);
  myHistos_->fillProfile("hProfRecoVsMagGen_"+hNameSuffix,
			 aGenTau2.getPCA().Mag(),
			 tau2PCA.Vect().Mag(),
			 eventWeight);
 
  ///////////////////////////////////////////////////////////
  float yTau1 =  2.*tau1LeadingTk.Pt()/aTau1.getP4().Pt() - 1.;
  float yTau2 =  2.*tau2LeadingTk.Pt()/aTau2.getP4().Pt() - 1.;
  float yTau = aTau1.getCharge()>0 ? yTau2 : yTau1; //FIXME: correct?
  float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;
 
  if(aTau2.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero && isOneProng(aTau2.getProperty(PropertyEnum::decayMode))){
     myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);
     myHistos_->fill1DHistogram("h1DPhi_nVecIP_"+hNameSuffix,shiftedIPrho,eventWeight);
     if(yTau>0){
       myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
     else{
       myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
     }
  }

  if(aTau2.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero){
    float cosPhiNN =  tau1PCA.Vect().Unit().Dot(tau2PCA.Vect().Unit());

    myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);
    myHistos_->fill1DHistogram("h1DCosPhiNN_"+hNameSuffix,cosPhiNN);
    
    float cosTau1 = tau1PCA.Vect().Unit()*aGenTau1.getPCA().Unit();
    float cosTau2 = tau2PCA.Vect().Unit()*aGenTau2.getPCA().Unit();
    
    myHistos_->fillProfile("hProfPhiVsMag_"+hNameSuffix,aGenTau1.getPCA().Mag(),cosTau1);
    myHistos_->fillProfile("hProfPhiVsMag_"+hNameSuffix,aGenTau2.getPCA().Mag(),cosTau2);       
    }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzerTT::fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight){
  
  ///Method from http://arxiv.org/abs/1108.0670 (Berger)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame  
  std::pair<float,float>  angles;

  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  std::pair<float,float>  anglesRhoRho;  
  ///Here we take rho on one side, pi+- on the other
  std::pair<float,float>  anglesIPRho;

  TLorentzVector tau1LeadingTk = aGenTau1.getChargedP4();
  TLorentzVector tau2LeadingTk = aGenTau2.getChargedP4();
  TLorentzVector tau1PCA(aGenTau1.getPCA(),0);
  TLorentzVector tau2PCA(aGenTau2.getPCA(),0);

  if(aTau1.getCharge()>0) {
    angles = angleBetweenPlanes(tau1LeadingTk, tau1PCA, tau2LeadingTk, tau2PCA);
    anglesIPRho = angleBetweenPlanes(tau1LeadingTk, tau1PCA, tau2LeadingTk, aGenTau2.getNeutralP4());
    anglesRhoRho = angleBetweenPlanes(tau1LeadingTk, aGenTau1.getNeutralP4(), tau2LeadingTk, aGenTau2.getNeutralP4());
  }
  else{
    angles = angleBetweenPlanes(tau2LeadingTk, tau2PCA, tau1LeadingTk, tau1PCA);
    anglesIPRho = angleBetweenPlanes(tau2LeadingTk, tau2PCA, tau1LeadingTk, aGenTau1.getNeutralP4());
    anglesRhoRho = angleBetweenPlanes(tau2LeadingTk, aGenTau2.getNeutralP4(), tau1LeadingTk, aGenTau1.getNeutralP4());
  }

  if(aGenTau2.getProperty(PropertyEnum::decayMode)==tauDecay1ChargedPion0PiZero){
    float cosPhiNN =  tau1PCA.Vect().Unit().Dot(tau2PCA.Vect().Unit());
    myHistos_->fill1DHistogram("h1DCosPhiNN_"+hNameSuffix,cosPhiNN);
    myHistos_->fill1DHistogram("h1DPhi_nVectors"+hNameSuffix,angles.first,eventWeight);
  }
  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho.
  ///Here we take rho on one side, mu on the other
  if(aGenTau2.getProperty(PropertyEnum::decayMode)!=tauDecay1ChargedPion0PiZero && isOneProng(aTau2.getProperty(PropertyEnum::decayMode))){
    
    float yTau =  2.*tau2LeadingTk.Pt()/(aGenTau2.getChargedP4()+aGenTau2.getNeutralP4()).Pt() - 1.;
    float shiftedIPrho = anglesIPRho.first +(yTau<0)*(1-2*(anglesIPRho.first>M_PI))*M_PI;
   
    myHistos_->fill1DHistogram("h1DyTau"+hNameSuffix,yTau,eventWeight);
    myHistos_->fill1DHistogram("h1DPhi_nVecIP_"+hNameSuffix,shiftedIPrho,eventWeight);
    if(yTau>0) myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauPos"+hNameSuffix,anglesIPRho.first,eventWeight);
    else myHistos_->fill1DHistogram("h1DPhi_nVecIP_yTauNeg"+hNameSuffix,anglesIPRho.first,eventWeight);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<float,float> HTTAnalyzerTT::angleBetweenPlanes(const TLorentzVector &tau1, 
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
