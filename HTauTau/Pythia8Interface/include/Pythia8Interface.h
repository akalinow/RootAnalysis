#ifndef RootAnalysis_Pythia8Interface_H
#define RootAnalysis_Pythia8Interface_H

#include <string>
#include <vector>
#include <map>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "HTTEvent.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"
#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TRandom3.h"

#include "Analyzer.h"


class TLorentzVector;
class TTree;
class TFile;
class TH1F;

class Pythia8Interface: public Analyzer{

 public:

  Pythia8Interface(const std::string & aName);

  virtual ~Pythia8Interface();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent) {return analyze(iEvent, NULL);}

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger);

  virtual void finalize();

  virtual void clear(){;};

  Analyzer* clone() const;

 private:

  void initializePythia(double mH, int decayMode);

  void getGenTaus();

  int getDetailedTauDecayMode(const TParticle & aTau) const;

  TLorentzVector getNeutralComponent(const TParticle & aTau) const;
  
  TLorentzVector getChargedComponent(const TParticle & aTau) const;

  TVector3 get3DImpactPoint(const TLorentzVector & aTauP4,
			    const TLorentzVector & aChargedP4,
			    const TVector3 & sv) const;
  
  void fillHTTEvent(unsigned long int eventNumber);

  bool checkDecayMode(const TParticle & aTau1,
		      const TParticle & aTau2,
		      int pairDecayMode);

  HTTParticle makeTau(const TParticle & aTau);

  HTTPair makePair(const HTTParticle & aTau1, const HTTParticle & aTau2);

  void makeRecoTaus(const TParticle & aTau1, const TParticle & aTau2);

  TRandom3 myRandGenerator;
  TPythia8 pythia8;
  TClonesArray myParticles;

  TParticle myTauPlus, myTauMinus;
  TParticle mySmearedLeg1, mySmearedLeg2; 

  std::vector<HTTPair> httPairCollection;
  std::vector<HTTParticle> httJetCollection;
  std::vector<HTTParticle> httLeptonCollection;
  std::vector<HTTParticle> httGenLeptonCollection;
  TTree *myTree;
  TFile *myFile;
  HTTEvent *httEvent;
  TH1F* hStats;

};

#endif
