#include <sstream>

#include "CPAnalyzer.h"
#include "EventProxyCPNtuple.h"



//declare helper functins
TVector3 impactParameter(const TVector3&, 
			 const TVector3&, 
			 const TLorentzVector&);
std::pair<float,float> angleBetweenPlanes(const TLorentzVector&, const TLorentzVector&,
					  const TLorentzVector&, const TLorentzVector&);
float deltaPhi(float, float);
float deltaR2(float, float, float, float);
float deltaR2(const TLorentzVector&, const TLorentzVector&);
enum tauDecayModes {kElectron, kMuon, 
		    kOneProng0pi0, kOneProng1pi0, kOneProng2pi0, kOneProng3pi0,
		    kThreeProng0pi0, kThreeProng1pi0,
		    kOther, kUndefined};
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
bool CPAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyCPNtuple & myEvent = static_cast<const EventProxyCPNtuple&>(iEvent);
  
  std::pair<float,float> angles, angles2;
  angles = angleBetweenPlanes(*myEvent.tauMinus,
			      *myEvent.piMinus,
			      *myEvent.tauPlus,
			      *myEvent.piPlus);


  myHistos_->fill1DHistogram("h1DPhi",angles.first);
  myHistos_->fill1DHistogram("h1DRho",angles.second);
  
  
  /*
  TVector3 nMinus = impactParameter(*thePV,*svMinus,*piMinus);
  TVector3 nPlus  = impactParameter(*thePV,*svPlus,*piPlus);
  angles2 = angleBetweenPlanes(*piMinus,TLorentzVector(nMinus,0),
			       *piPlus,TLorentzVector(nPlus,0));
  
  hPhi->Fill(angles.first); hRho->Fill(angles.second);      
  hDPhi->Fill(deltaPhi(angles.first,phi));
  hPhi2->Fill(angles2.first); hRho2->Fill(angles2.second);
  hDPhi2->Fill(deltaPhi(angles2.first,phi2));
  */
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TVector3 impactParameter(const TVector3& pv, 
			 const TVector3& sv, 
			 const TLorentzVector& p4){
  
  TVector3 dir = (p4.Vect()).Unit();
  TVector3 n = (sv-pv) - ((sv-pv)*dir)*dir;

  return n;
}

std::pair<float,float> angleBetweenPlanes(const TLorentzVector &prod1, 
					  const TLorentzVector &prod12,
					  const TLorentzVector &prod2, 
					  const TLorentzVector &prod22){
  //Boost to correct frame
  TVector3 boost = (prod1+prod2).BoostVector();
  TLorentzVector prod1Star = prod1; prod1Star.Boost(-boost);
  TLorentzVector prod12Star = prod12; prod12Star.Boost(-boost);
  TLorentzVector prod2Star = prod2; prod2Star.Boost(-boost);
  TLorentzVector prod22Star = prod22; prod22Star.Boost(-boost);

  //define common direction and normal vectors to decay planes
  TVector3 direction = prod1Star.Vect().Unit();
  TVector3 n1 = ( direction.Cross( prod12Star.Vect() ) ).Unit(); 
  TVector3 n2 = ( direction.Cross( prod22Star.Vect() ) ).Unit(); 

  float phi=TMath::ACos(n1*n2);
  float rho=TMath::ACos( (prod12Star.Vect().Unit() )*(prod22Star.Vect().Unit() ) );

  return std::make_pair(phi,rho);
}

float deltaPhi(float phi1, float phi2) { 
  float result = phi1 - phi2;
  while( result > TMath::Pi() ) result -= 2.*TMath::Pi();
  while( result <= -TMath::Pi() ) result += 2.*TMath::Pi();
  return result;
}
float deltaR2(float phi1, float eta1, float phi2, float eta2){
  return (eta1-eta2)*(eta1-eta2)+deltaPhi(phi1,phi2)*deltaPhi(phi1,phi2);
}

float deltaR2(const TLorentzVector &p41, const TLorentzVector &p42){
  return deltaR2(p41.Phi(),p41.Eta(),p42.Phi(),p42.Eta());
 }
