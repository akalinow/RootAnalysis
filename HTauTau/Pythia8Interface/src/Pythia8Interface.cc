#include <sstream>
#include <bitset>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TParticlePDG.h"

#include "Pythia8Interface.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Pythia8Interface::Pythia8Interface(const std::string & aName) : Analyzer(aName), covMET(2,2){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Pythia8Interface::~Pythia8Interface(){

  if(myFile){
    myFile->Write();
    delete myFile;
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* Pythia8Interface::clone() const {

  Pythia8Interface* clone = new Pythia8Interface(name());
  return clone;
};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::initialize(TDirectory* aDir,
                             pat::strbitset *aSelections){

  //initializePythia(125, HTTAnalysis::hadronicTauDecayModes::tauDecayMuon);

  myParticles.SetClass("TParticle", 2000);

  std::string fileName = "Pythia8_SVFitData.root";
  myFile = new TFile(fileName.c_str(),"RECREATE");
  httEvent = new HTTEvent();
  myTree = new TTree("HTauTauTree","");
  myTree->SetDirectory(myFile);
  myTree->Branch("HTTEvent.",&httEvent);
  myTree->Branch("HTTPairCollection",&httPairCollection);
  myTree->Branch("HTTJetCollection",&httJetCollection);
  myTree->Branch("HTTLeptonCollection",&httLeptonCollection);
  myTree->Branch("HTTGenLeptonCollection",&httGenLeptonCollection);
  hStats = new TH1F("hStats","Bookkeeping histogram",11,-0.5,10.5);
  hStats->SetDirectory(myFile);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::finalize(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::initializePythia(double mH, int decayMode){

  //pythia8.ReadString("Init:showChangedSettings = off");
  //pythia8.ReadString("Init:showChangedParticleData = off");
  pythia8.ReadString("Init:showProcesses = off");
  pythia8.ReadString("Init:showMultipartonInteractions = off");
  /*
    pythia8.ReadString("Higgs:useBSM  = on");
    pythia8.ReadString("HiggsBSM:gg2A3 = on"); //Higgs production by gluon-gluon fusion
    pythia8.ReadString("36:m0 = 125");       //Higgs mass
    pythia8.ReadString("36:onMode = no");    //switch off all Higgs decay channels
    pythia8.ReadString("36:onIfMatch =  15 15"); //switch back on Higgs -> tau tau
  */

  std::string tmpString = "25:m0 = " + std::to_string(mH);
  pythia8.ReadString(tmpString.c_str());       //Higgs mass
  tmpString = "25:mMin = " + std::to_string(mH-5.0);
  pythia8.ReadString(tmpString.c_str());       //Higgs Breit-Wigner range lower edge
  tmpString = "25:mMax = " + std::to_string(mH+5.0);
  pythia8.ReadString(tmpString.c_str());       //Higgs Breit-Wigner range upper edge

  pythia8.ReadString("HiggsSM:gg2H = on"); //Higgs production by gluon-gluon fusion  
  pythia8.ReadString("25:onMode = no");    //switch off all Higgs decay channels
  pythia8.ReadString("25:onIfMatch = 15 -15"); //switch back on Higgs -> tau tau

  //pythia8.ReadString("15:onMode = no");    //switch off all tau decay channels
  //pythia8.ReadString("15:onIfMatch =  211 16"); //switch back on tau -> pi nu

  pythia8.ReadString("15:offIfAny = 11");//switch off tau -> e nu

  if(decayMode!=HTTAnalysis::hadronicTauDecayModes::tauDecayMuon){
    pythia8.ReadString("15:offIfAny = 13");//switch off tau -> mu nu
  }   
  //pythia8.ReadString("15:offIfAny = 211"); //switch back on tau -> pi nu
    
  pythia8.Initialize(2212 /* p */, 2212 /* p */, 13000. /* TeV */);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int Pythia8Interface::getDetailedTauDecayMode(const TParticle & aTau) const{

  int tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayOther;
  
  int numElectrons      = 0;
  int numMuons          = 0;
  int numChargedPions   = 0;
  int numNeutralPions   = 0;
  int numPhotons        = 0;
  int numNeutrinos      = 0;
  int numOtherParticles = 0;

  TParticle* daughter = 0;  
  for(int iDaughter=aTau.GetFirstDaughter();
      iDaughter<=aTau.GetLastDaughter(); ++iDaughter){

    daughter = (TParticle*) myParticles.At(iDaughter);    
    int pdg_id = std::abs(daughter->GetPdgCode());    
    if(pdg_id == 11) numElectrons++;
    else if(pdg_id == 13) numMuons++;
    else if(pdg_id == 211 || pdg_id == 321 ) numChargedPions++; //Count both pi+ and K+
    else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 ) numNeutralPions++; //Count both pi0 and K0_L/S
    else if(pdg_id == 12 || 
	    pdg_id == 14 || 
	    pdg_id == 16) {
      numNeutrinos++;
    }
    else if(pdg_id == 22) numPhotons++;
    else {
      numOtherParticles++;
    }
  }

  ///Pythia does not decay pi0.
  numPhotons += 2*numNeutralPions; 
 
  if(numElectrons>1){//sometimes there are gamma->ee conversions 
    numPhotons += numElectrons/2;
    numElectrons -= 2*(numElectrons/2);
  }
  
  if( numOtherParticles == 0 ){
    if( numElectrons == 1 ){
      //--- tau decays into electrons
      tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecaysElectron;
    } else if( numMuons == 1 ){
      //--- tau decays into muons
      tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayMuon;
    } else {
      //--- hadronic tau decays
      switch ( numChargedPions ){
      case 1 :
	if( numNeutralPions != 0 ){
	  ///Pythia level does not decay pi0
	  //tauDecayMode =  HTTAnalysis::hadronicTauDecayModes::tauDecayOther;
	  //break;
	}
	switch ( numPhotons ){
	case 0:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay1ChargedPion0PiZero;
	  break;
	case 2:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay1ChargedPion1PiZero;
	  break;
	case 4:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay1ChargedPion2PiZero;
	  break;
	case 6:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay1ChargedPion3PiZero;
	  break;
	case 8:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay1ChargedPion4PiZero;
	  break;
	default:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayOther;
	  break;
	}
	break;
      case 3 : 
	if( numNeutralPions != 0 ){
	  //tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayOther;
	  //break;
	}
	switch ( numPhotons ){
	case 0 : 
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion0PiZero;
	  break;
	case 2 : 
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion1PiZero;
	  break;
	case 4 : 
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion2PiZero;
	  break;
	case 6 : 
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion3PiZero;
	  break;
	case 8 : 
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecay3ChargedPion4PiZero;
	  break;
	default:
	  tauDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayOther;
	  break;
	}
	break;
      }
    }
  } 
  return tauDecayMode;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::getGenTaus(){

  int numberOfParticles = myParticles.GetEntriesFast();
  
  TParticle* aPart = 0;
  for (int iParticle=0;iParticle<numberOfParticles;++iParticle){
    aPart = (TParticle*) myParticles.At(iParticle);
    int pdgCode = aPart->GetPdgCode();
    if(pdgCode==15) myTauMinus = *aPart;
    if(pdgCode==-15) myTauPlus = *aPart;
  }      
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Pythia8Interface::getChargedComponent(const TParticle & aTau) const{

  TLorentzVector aP4;
  TParticle* daughter = 0;
  int charge = 0;
  
  for(int iDaughter=aTau.GetFirstDaughter();
      iDaughter<=aTau.GetLastDaughter(); ++iDaughter){
    daughter = (TParticle*) myParticles.At(iDaughter);
    charge = daughter->GetPDG()->Charge();
    if(charge) aP4+=TLorentzVector(daughter->Px(),
				   daughter->Py(),
				   daughter->Pz(),
				   daughter->Energy());
  }
  return aP4;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TLorentzVector Pythia8Interface::getNeutralComponent(const TParticle & aTau) const{

  TLorentzVector aP4;
  TParticle* daughter = 0;
  int pdgID = 0;
  int charge = 1;
  for(int iDaughter=aTau.GetFirstDaughter();
      iDaughter<=aTau.GetLastDaughter(); ++iDaughter){
    daughter = (TParticle*) myParticles.At(iDaughter);
    pdgID = daughter->GetPdgCode();
    charge = daughter->GetPDG()->Charge();
    if(pdgID != 12 && pdgID != 14 && pdgID != 16 && charge==0)  aP4+=TLorentzVector(daughter->Px(),
										    daughter->Py(),
										    daughter->Pz(),
										    daughter->Energy());
  }
  return aP4;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TVector3 Pythia8Interface::get3DImpactPoint(const TLorentzVector & aTauP4,
					    const TLorentzVector & aChargedP4,
					    const TVector3 & sv) const{

  double cosGJ = aTauP4.Vect().Unit().Dot(aChargedP4.Vect().Unit());
  double ipDistance = sv.Mag()*sqrt(1.0 - cosGJ*cosGJ);
  TVector3 ip = sv - aChargedP4.Vect().Unit()*ipDistance;

  return ip;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::fillHTTEvent(unsigned long int eventNumber){

  int runNumber = 1;
  httEvent->setRun(runNumber);
  httEvent->setEvent(eventNumber);
  httEvent->setNPV(1.0);
  //httEvent->setGenPV(pvGen);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool Pythia8Interface::checkDecayMode(const TParticle & aTau1,
				      const TParticle & aTau2,
				      int pairDecayMode){

 int decayMode1 = getDetailedTauDecayMode(aTau1);
 int decayMode2 = getDetailedTauDecayMode(aTau2);
 
 if(decayMode1!=pairDecayMode &&
    decayMode2!=pairDecayMode) return false;

 if(decayMode1==pairDecayMode &&
    decayMode2==pairDecayMode) return false;

 return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTParticle Pythia8Interface::makeTau(const TParticle & aTau){

  int daughterIndex = aTau.GetFirstDaughter();
  TParticle* aDaughter = (TParticle*) myParticles.At(daughterIndex);

  TLorentzVector p4(aTau.Px(), aTau.Py(), aTau.Pz(), aTau.Energy());
  TLorentzVector p4Charged = getChargedComponent(aTau);
  TLorentzVector p4Neutral = getNeutralComponent(aTau);
  TVector3 pv(aTau.Vx(), aTau.Vy(), aTau.Vz());
  pv *= 0.1;//convert mm [Pythia] to cm [CMSSW]
  TVector3 sv(aDaughter->Vx(), aDaughter->Vy(), aDaughter->Vz());
  sv *= 0.1;//convert mm [Pythia] to cm [CMSSW]
  TVector3 pca = get3DImpactPoint(p4, p4Charged, sv);
  
  HTTParticle aLepton;
		    
  aLepton.setP4(p4);
  aLepton.setChargedP4(p4Charged);
  aLepton.setNeutralP4(p4Neutral);
  aLepton.setPCA(pca);
  aLepton.setSV(sv);

  std::vector<Double_t> aProperties;
  int decayMode = getDetailedTauDecayMode(aTau);

  int pdgId = 15;
  if(decayMode == HTTAnalysis::hadronicTauDecayModes::tauDecaysElectron) pdgId = 11;
  else if(decayMode == HTTAnalysis::hadronicTauDecayModes::tauDecayMuon) pdgId = 16;
  else pdgId = 15;
  
  aProperties.push_back(pdgId);
  aProperties.push_back(aTau.GetPDG()->Charge());
  aProperties.push_back(decayMode);

  aLepton.setProperties(aProperties);
  
  return aLepton;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::makeMET(const HTTParticle & aTau1, const HTTParticle & aTau2){
  
  TLorentzVector nunuGen = aTau1.getP4() + aTau2.getP4() -
    aTau1.getChargedP4() - aTau2.getChargedP4() -
    aTau1.getNeutralP4() - aTau2.getNeutralP4();

  myMET.SetX(nunuGen.X());
  myMET.SetY(nunuGen.Y());

  covMET[0][0] = 1.0;
  covMET[0][1] = 0.0;
  covMET[1][0] = 0.0;
  covMET[1][1] = 1.0;
      
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void Pythia8Interface::makeRecoTaus(const TParticle & aTau1, const TParticle & aTau2){

  int decayModeTau1 = getDetailedTauDecayMode(aTau1);
  int decayModeTau2 = getDetailedTauDecayMode(aTau2);

  if(decayModeTau1==HTTAnalysis::hadronicTauDecayModes::tauDecayMuon &&
     decayModeTau2!=HTTAnalysis::hadronicTauDecayModes::tauDecayMuon){
    mySmearedLeg1 = aTau1;
    mySmearedLeg2 = aTau2;
  }
  else if(decayModeTau1!=HTTAnalysis::hadronicTauDecayModes::tauDecayMuon &&
	  decayModeTau2==HTTAnalysis::hadronicTauDecayModes::tauDecayMuon){
    mySmearedLeg1 = aTau1;
    mySmearedLeg2 = aTau2;
  }
  else if(aTau1.Pt()>aTau2.Pt()){
    mySmearedLeg1 = aTau1;
    mySmearedLeg2 = aTau2;
  }
  else{
    mySmearedLeg1 = aTau2;
    mySmearedLeg2 = aTau1;
  }

  TLorentzVector p4Vis = getChargedComponent(mySmearedLeg1) + getNeutralComponent(mySmearedLeg1);
  mySmearedLeg1.SetMomentum(p4Vis);

  p4Vis = getChargedComponent(mySmearedLeg2) + getNeutralComponent(mySmearedLeg2);
  mySmearedLeg2.SetMomentum(p4Vis);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTPair Pythia8Interface::makePair(const HTTParticle & aTau1, const HTTParticle & aTau2,
				   const TVector2 & met, const TMatrixD & covMET){

  TLorentzVector p4 = aTau1.getP4() + aTau2.getP4();

  double mTLeg1 = 0.0;//FIXME
  double mTLeg2 = 0.0;//FIXME

  HTTPair aHTTpair;
  aHTTpair.setP4(p4);
  aHTTpair.setMET(met);
  aHTTpair.setMETMatrix(covMET[0][0], covMET[0][1], covMET[1][0], covMET[1][1]);
  aHTTpair.setMTLeg1(mTLeg1);
  aHTTpair.setMTLeg2(mTLeg2);
  aHTTpair.setLeg1(aTau1);
  aHTTpair.setLeg2(aTau2);

  return aHTTpair;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool Pythia8Interface::analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){

  double mH = 125;
  double tmp = 0;

  unsigned int long eventsToGenerate = *aMessenger->getObject(&tmp,"eventsToGenerate");  
  int pairDecayMode = HTTAnalysis::hadronicTauDecayModes::tauDecayMuon;

  for(int iMass=0;iMass<100;++iMass){
    mH = 50 + 2*iMass;
    std::cout<<"Generating mass: "<<mH<<std::endl;
    initializePythia(mH, pairDecayMode);
    for(unsigned int long eventNumber=0;eventNumber<eventsToGenerate;++eventNumber){
      
      pythia8.GenerateEvent();
      //pythia8.EventListing(); 
      pythia8.ImportParticles(&myParticles,"All");      
      getGenTaus();
      bool isCorrectPairDecayMode = checkDecayMode(myTauPlus, myTauMinus, pairDecayMode);
      if(!isCorrectPairDecayMode) continue;
      
      makeRecoTaus(myTauPlus, myTauMinus);
      
      fillHTTEvent(eventNumber);
      
      httGenLeptonCollection.clear();
      HTTParticle aGenTau1 = makeTau(myTauPlus);
      HTTParticle aGenTau2 = makeTau(myTauMinus);
      makeMET(aGenTau1, aGenTau2);
      httGenLeptonCollection.push_back(aGenTau1);
      httGenLeptonCollection.push_back(aGenTau2);

      httLeptonCollection.clear();
      HTTParticle aLeg1 = makeTau(mySmearedLeg1);
      HTTParticle aLeg2 = makeTau(mySmearedLeg2);
      httLeptonCollection.push_back(aLeg1);
      httLeptonCollection.push_back(aLeg2);

      httPairCollection.clear();
      httPairCollection.push_back(makePair(aLeg1, aLeg2, myMET, covMET));
      
      myTree->Fill();
      hStats->Fill(0);      
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
