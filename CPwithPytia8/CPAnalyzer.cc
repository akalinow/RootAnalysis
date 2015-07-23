#include <sstream>

#include "CPAnalyzer.h"
#include "EventProxyCPNtuple.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
CPAnalyzer::CPAnalyzer(const std::string & aName):Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
CPAnalyzer::~CPAnalyzer(){

  if(myHistos_) delete myHistos_;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void CPAnalyzer::initialize(TFileDirectory& aDir,
			    pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new CPHistograms(&aDir, selectionFlavours_);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void CPAnalyzer::finalize(){ 

 myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void CPAnalyzer::fillAngles(const EventProxyCPNtuple & myEvent,
			    const std::string & sysType){

  ///Angles between decay planes, and decay products
  std::pair<float,float> angles = angleBetweenPlanes(*myEvent.tauMinus, *myEvent.piMinus,
						     *myEvent.tauPlus, *myEvent.piPlus);
  
  myHistos_->fill1DHistogram("h1DPhi"+sysType,angles.first);
  myHistos_->fill1DHistogram("h1DRho"+sysType,angles.second);

  TVector3 pv = *myEvent.thePV;
  TVector3 svMinus = *myEvent.svMinus;
  TVector3 svPlus = *myEvent.svPlus;

  ///Smear pv and sv
  double pullX = 0, pullZ = 0;
  
  if(sysType.find("smearPV")!=std::string::npos){
    double sigmaX = 12E-3, sigmaY = 12E-3, sigmaZ = 30E-3;//[mm]
    pullX = aRndm.Gaus(0.0,sigmaX);
    pullZ = aRndm.Gaus(0.0,sigmaZ);
    pv.SetXYZ(pv.X()+pullX,
	      pv.Y()+aRndm.Gaus(0.0,sigmaY),
	      pv.Z()+pullZ);
  }
  //////////////////
  myHistos_->fill1DHistogram("h1DVxPullX"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullZ"+sysType,pullZ);
  
  TVector3 nMinus = impactParameter(pv,svMinus,*myEvent.piMinus);
  TVector3 nPlus  = impactParameter(pv,svPlus,*myEvent.piPlus);
  
  if(sysType.find("smear")!=std::string::npos &&
     sysType.find("PCA")!=std::string::npos){
    double sigmaX = 10E-3, sigmaY = 10E-3, sigmaZ = 30E-3; //[mm]
    /*
    nMinus.SetPerp(nMinus.Perp()+aRndm.Gaus(0.0,sigmaX));
    nMinus.SetZ(nMinus.Z()+aRndm.Gaus(0.0,sigmaZ));

    nPlus.SetPerp(nPlus.Perp()+aRndm.Gaus(0.0,sigmaX));
    nPlus.SetZ(nPlus.Z()+aRndm.Gaus(0.0,sigmaZ));
    */
    
    nMinus.SetXYZ(nMinus.X()+aRndm.Gaus(0.0,sigmaX),
		  nMinus.Y()+aRndm.Gaus(0.0,sigmaY),
		  nMinus.Z()+aRndm.Gaus(0.0,sigmaZ));
    
    nPlus.SetXYZ(nPlus.X()+aRndm.Gaus(0.0,sigmaX),
		 nPlus.Y()+aRndm.Gaus(0.0,sigmaY),
		 nPlus.Z()+aRndm.Gaus(0.0,sigmaZ));    
  }
  ////////////////////////////

  
  float cosPhiTauMinusPi = myEvent.tauMinus->Vect().Unit().Dot(myEvent.piMinus->Vect().Unit());
  float cosPhiTauPlusPi = myEvent.tauPlus->Vect().Unit().Dot(myEvent.piPlus->Vect().Unit());
  
  myHistos_->fill1DHistogram("h1DIPMinus"+sysType,nMinus.Mag());
  myHistos_->fill1DHistogram("h1DIPPlus"+sysType,nPlus.Mag());

  ///Method from http://arxiv.org/abs/1108.0670
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles2 = angleBetweenPlanes(*myEvent.piMinus,TLorentzVector(nMinus,0),
						       *myEvent.piPlus,TLorentzVector(nPlus,0));

  myHistos_->fill1DHistogram("h1DPhi2"+sysType,angles2.first);
  myHistos_->fill1DHistogram("h1DRho2"+sysType,angles2.second);


  //take visible tau instead full tau momentum.
  //calculate angles in pi+ - pi- rest frame
  std::pair<float,float> anglesRho = angleBetweenPlanes(*myEvent.visTauMinus, *myEvent.piMinus,
							*myEvent.visTauPlus,  *myEvent.piPlus);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyCPNtuple & myEvent = static_cast<const EventProxyCPNtuple&>(iEvent);


  if( !(isOneProng(myEvent.decModeMinus) || isLepton(myEvent.decModeMinus) ) ||
      !(isOneProng(myEvent.decModePlus)  || isLepton(myEvent.decModePlus) ) ) return true;

  //if( (*myEvent.piMinus+*myEvent.piPlus).Perp()>100) return true;

  //if(myEvent.piMinus->Perp()<20 || fabs(myEvent.piMinus->Eta())>2.1) return true;
  //if(myEvent.piPlus->Perp()<20 || fabs(myEvent.piPlus->Eta())>2.1) return true;

  //if(myEvent.tauMinus->Perp()<20 || fabs(myEvent.tauMinus->Eta())>2.1) return true;
  //if(myEvent.tauPlus->Perp()<20 || fabs(myEvent.tauPlus->Eta())>2.1) return true;

  
  std::string smearType = "ideal";
  std::string decayName = getDecayName(myEvent.decModeMinus, myEvent.decModePlus);
  std::string motherName = getMotherName(myEvent.bosonId);
  std::string name = motherName+"_"+decayName+"_"+smearType;
  
  fillAngles(myEvent, name);

  


  /*
  sysType = "smearPV";
  fillAngles(myEvent, sysType);
  sysType = "smearPV_PCA";
  fillAngles(myEvent, sysType);
  */


  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::pair<float,float> CPAnalyzer::angleBetweenPlanes(const TLorentzVector &tau1, 
					  const TLorentzVector &tau1Daughter,
					  const TLorentzVector &tau2, 
					  const TLorentzVector &tau2Daughter){
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
  ///angle between tau1 and tau2 daughter momentums
  float rho=TMath::ACos( (tau1DaughterStar.Vect().Unit() )*(tau2DaughterStar.Vect().Unit() ) );

  return std::make_pair(phi,rho);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TVector3 CPAnalyzer::impactParameter(const TVector3& pv, 
				     const TVector3& sv, 
				     const TLorentzVector& p4){
  
  TVector3 dir = (p4.Vect()).Unit();
  TVector3 n = (sv-pv) - ((sv-pv)*dir)*dir;

  return n;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string CPAnalyzer::getDecayName(int decModeMinus, int decModePlus){


  if(decModeMinus==kOneProng0pi0 && decModePlus==kOneProng0pi0) return "PiPi0Pi0";

  if(isOneProng(decModeMinus) && isOneProng(decModePlus) ) return "1Prong1Prong";

  if( (decModeMinus==kOneProng0pi0 && isLepton(decModePlus) ) ||
      (isLepton(decModeMinus) && decModePlus==kOneProng0pi0)) return "Lepton1Prong0Pi0";
    
  if( (isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
      ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ) return "Lepton1Prong";

  if(decModeMinus==kOneProng1pi0 && decModePlus==kOneProng1pi0 ) return "PiPlusPiMinus2Pi0";

  if( isOneProng(decModeMinus) && decModeMinus!=kOneProng0pi0 && 
      isOneProng(decModePlus) && decModePlus!=kOneProng0pi0 ) return "1Prong1ProngXPi0";

  if(isLepton(decModePlus) && isLepton(decModeMinus))  return "LeptonLepton";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string CPAnalyzer::getMotherName(int bosonId){

  if(abs(bosonId) == 0) return "Unknown";
  if(abs(bosonId) == 23) return "Z0";
  if(abs(bosonId) == 25) return "h0";
  if(abs(bosonId) == 35) return "H0";
  if(abs(bosonId) == 36) return "A0";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::isOneProng(int decMode){
  if(decMode==kOneProng0pi0 ||
     decMode==kOneProng1pi0 ||
     decMode==kOneProng2pi0 ||
     decMode==kOneProng3pi0 ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::isLepton(int decMode){
  if(decMode==kElectron || decMode==kMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////







float deltaPhi(float phi1, float phi2) { 
  float result = phi1 - phi2;
  while( result > TMath::Pi() ) result -= 2.*TMath::Pi();
  while( result <= -TMath::Pi() ) result += 2.*TMath::Pi();
  return result;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float deltaR2(float phi1, float eta1, float phi2, float eta2){
  return (eta1-eta2)*(eta1-eta2)+deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float deltaR2(const TLorentzVector &p41, const TLorentzVector &p42){
  return deltaR2(p41.Phi(),p41.Eta(),p42.Phi(),p42.Eta());
 }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
