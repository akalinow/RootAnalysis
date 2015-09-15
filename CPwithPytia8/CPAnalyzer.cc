#include <sstream>

#include "CPAnalyzer.h"
#include "EventProxyCPNtuple.h"
#include "HTTEvent.h"

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
void CPAnalyzer::fillAngles(const DiTauData* aEvent,
			    const std::string & sysType){

  ///Angles between decay planes, and decay products
  ///Angle betweeny decay products in tau-tau rest frame: http://arxiv.org/abs/hep-ph/0503172, pp. 81
  std::pair<float,float> angles = angleBetweenPlanes(aEvent->tauMinus_, aEvent->piMinus_,
						     aEvent->tauPlus_, aEvent->piPlus_);
  
  myHistos_->fill1DHistogram("h1DPhi"+sysType,angles.first);
  myHistos_->fill1DHistogram("h1DRho"+sysType,angles.second);

  TVector3 pv = aEvent->thePV_;
  TVector3 svMinus = aEvent->svMinus_;
  TVector3 svPlus = aEvent->svPlus_;

  ///Smear pv and sv
  double pullX = 0, pullZ = 0;
  
  if(sysType.find("smearPV")!=std::string::npos){
    double sigmaX = 10E-4, sigmaY = 10E-4, sigmaZ = 30E-4;//[cm]?
    pullX = aRndm.Gaus(0.0,sigmaX);
    pullZ = aRndm.Gaus(0.0,sigmaZ);
    pv.SetXYZ(pv.X()+pullX,
	      pv.Y()+aRndm.Gaus(0.0,sigmaY),
	      pv.Z()+pullZ);
  }
  //////////////////
  myHistos_->fill1DHistogram("h1DVxPullX"+sysType,pullX);
  myHistos_->fill1DHistogram("h1DVxPullZ"+sysType,pullZ);
  
  //TVector3 nMinus = impactParameter(pv,svMinus,aEvent->piMinus_);
  //TVector3 nPlus  = impactParameter(pv,svPlus,aEvent->piPlus_);
  TVector3 nMinus = aEvent->nPiMinus_;
  TVector3 nPlus  = aEvent->nPiPlus_;
  
  if(sysType.find("smear")!=std::string::npos &&
     sysType.find("PCA")!=std::string::npos){
    double sigmaX = 20E-4, sigmaY = 20E-4, sigmaZ = 20E-4; //[cm]?
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

  
  float cosPhiTauMinusPi = (svMinus-pv).Unit().Dot(aEvent->piMinus_.Vect().Unit());
  float cosPhiTauPlusPi = (svPlus-pv).Unit().Dot(aEvent->piPlus_.Vect().Unit());

  myHistos_->fill1DHistogram("h1DIP_PCA"+sysType,nMinus.Mag());
  myHistos_->fill1DHistogram("h1DIP_3DIP"+sysType,(svMinus-pv).Mag());

  myHistos_->fill1DHistogram("h1DCosPhi_collinearMinus"+sysType,cosPhiTauMinusPi);
  myHistos_->fill1DHistogram("h1DCosPhi_collinearPlus"+sysType,cosPhiTauPlusPi);

  ///Method from http://arxiv.org/abs/1108.0670 (Berger)
  ///take impact parameters instead of tau momentum.
  ///calculate angles in pi+ - pi- rest frame
  std::pair<float,float>  angles2 = angleBetweenPlanes(aEvent->piMinus_,TLorentzVector(nMinus,0),
  						       aEvent->piPlus_,TLorentzVector(nPlus,0));

  //if(cosPhiTauMinusPi>0.9999 && cosPhiTauPlusPi>0.9999){
  myHistos_->fill1DHistogram("h1DPhi_nVectors"+sysType,angles2.first);
  myHistos_->fill1DHistogram("h1DRho_nVectors"+sysType,angles2.second);
  //}

  /*
  //Method from http://arxiv.org/abs/hep-ph/0204292 (Was)
  //Angle between rho decay planes in the rho-rho 
  //system.
  float rho = angleBetweenPlanes(aEvent->visTauMinus_, aEvent->piMinus_,
				 aEvent->visTauPlus_,  aEvent->piPlus_).first;

  ///Take Monte Carlo level y1, y2 (calculated in tau rest frames)
  if(aEvent->yPlus_*aEvent.yMinus_>0) myHistos_->fill1DHistogram("h1DPhi_Rho_y1y2Plus"+sysType,rho);
  //else myHistos_->fill1DHistogram("h1DPhi_Rho_y1y2Minus"+sysType,rho);

  myHistos_->fill1DHistogram("h1DPhi_Rho_y1y2Minus"+sysType,rho*aEvent->yPlus_*aEvent.yMinus_);

  ///Take reconstruction level y1, y2 (calculated in LAB)
  if(aEvent->yPlusLab2_*aEvent->yMinusLab2_>0) myHistos_->fill1DHistogram("h1DPhi_Rho_y1y2PlusLAB"+sysType,rho);
  else myHistos_->fill1DHistogram("h1DPhi_Rho_y1y2MinusLAB"+sysType,rho);
  */

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyCPNtuple & myEvent = static_cast<const EventProxyCPNtuple&>(iEvent);
  DiTauData* aEventGen = &(myEvent.event->genEvent_);
  DiTauData* aEventReco = &(myEvent.event->recoEvent_);
  DiTauData* aEvent =  aEventReco;

  
  //Skip 3-prong and unknown decays.
  if( !(isOneProng(aEventGen->decModeMinus_) || isLepton(aEventGen->decModeMinus_) ) ||
      !(isOneProng(aEventGen->decModePlus_)  || isLepton(aEventGen->decModePlus_) ) ) return true;
  ///Skip lepton decays
  if(isLepton(aEventGen->decModeMinus_) || isLepton(aEventGen->decModePlus_)) return true;
  
  bool accepted = analysisSelection(aEvent);
  accepted = true;//TEST
  ///
  std::vector<std::string> decayNames = getDecayName(aEventGen->decModeMinus_, aEventGen->decModePlus_);
  
  std::string motherName = getMotherName(myEvent.event->bosonId_);
  std::string smearType = "ideal";
  std::string name;
  ///
  for(auto decayName:decayNames){
    smearType = "ideal";
    name = "_"+motherName+"_"+decayName+"_"+smearType;
    fillAngles(aEvent,name);
    if(accepted) name+="_selected";
    fillAngles(aEvent, name);
    smearType = "smearPV";
    name = "_"+motherName+"_"+decayName+"_"+smearType;
    fillAngles(aEvent, name);
    smearType = "smearPV_PCA";
    name = "_"+motherName+"_"+decayName+"_"+smearType;
    fillAngles(aEvent, name);
    if(accepted) name+="_selected";
    fillAngles(aEvent, name);
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::analysisSelection(const DiTauData * aEvent){
  // Minimal kinematic selection
  // lepton+lepton
  if(isLepton(aEvent->decModeMinus_) && !isLepton(aEvent->decModePlus_) ){
    if(aEvent->visTauMinus_.Pt()>20 &&
       std::abs(aEvent->visTauMinus_.Eta())<2.4 &&
       aEvent->visTauPlus_.Pt()>20 &&
       std::abs(aEvent->visTauPlus_.Eta())<2.4 ) return true;
    else return false;
  }
  // lepton+tau
  if( isLepton(aEvent->decModeMinus_) && !isLepton(aEvent->decModePlus_) ){
    if(aEvent->visTauMinus_.Pt()>20 &&
       std::abs(aEvent->visTauMinus_.Eta())<2.1 &&
       aEvent->visTauPlus_.Pt()>25 &&
       std::abs(aEvent->visTauPlus_.Eta())<2.3 ) return true;
    else return false;
  }
  if( !isLepton(aEvent->decModeMinus_) && isLepton(aEvent->decModePlus_) ){
    if(aEvent->visTauMinus_.Pt()>25 &&
       std::abs(aEvent->visTauMinus_.Eta())<2.3 &&
       aEvent->visTauPlus_.Pt()>20 &&
       std::abs(aEvent->visTauPlus_.Eta())<2.1 ) return true;
    else return false;
  }
  // tau+tau
  if( !isLepton(aEvent->decModeMinus_) && !isLepton(aEvent->decModePlus_) ){
    if(aEvent->visTauMinus_.Pt()>40 &&
       std::abs(aEvent->visTauMinus_.Eta())<2.3 &&
       aEvent->visTauPlus_.Pt()>40 &&
       std::abs(aEvent->visTauPlus_.Eta())<2.3 ) return true;
    else return false;
  }
  //undefined
  return false;
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
std::vector<std::string> CPAnalyzer::getDecayName(int decModeMinus, int decModePlus){

  std::vector<std::string> types;

  if(decModeMinus==tauDecay1ChargedPion0PiZero && decModePlus==tauDecay1ChargedPion0PiZero) types.push_back("PiPi0Pi0");

  if(isOneProng(decModeMinus) && isOneProng(decModePlus) ) types.push_back("1Prong1Prong");

  if( (decModeMinus==tauDecay1ChargedPion0PiZero && isLepton(decModePlus) ) ||
      (isLepton(decModeMinus) && decModePlus==tauDecay1ChargedPion0PiZero)) types.push_back("Lepton1Prong0Pi0");
    
  if( (isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
      ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ) types.push_back("Lepton1Prong");

  if(decModeMinus==tauDecay1ChargedPion1PiZero && decModePlus==tauDecay1ChargedPion1PiZero ) types.push_back("PiPlusPiMinus2Pi0");


  if( isOneProng(decModeMinus) && decModeMinus!=tauDecay1ChargedPion0PiZero && 
      isOneProng(decModePlus) && decModePlus!=tauDecay1ChargedPion0PiZero )   types.push_back("1Prong1ProngXPi0");

  if(isLepton(decModePlus) && isLepton(decModeMinus)) types.push_back("LeptonLepton");

  return types;
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
  if(decMode==tauDecay1ChargedPion0PiZero ||
     decMode==tauDecay1ChargedPion1PiZero ||
     decMode==tauDecay1ChargedPion2PiZero ||
     decMode==tauDecay1ChargedPion3PiZero ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool CPAnalyzer::isLepton(int decMode){
  if(decMode==tauDecaysElectron || decMode==tauDecayMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




