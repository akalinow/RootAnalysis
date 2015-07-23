#ifndef RootAnalysis_CPAnalyzer_H
#define RootAnalysis_CPAnalyzer_H

#include <string>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"

#include "strbitset.h"
#include "TFileDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"
#include "TRandom3.h"

#include "Analyzer.h"
#include "CPHistograms.h"

class EventProxyCPNtuple;
class TLorentzVector;
class TVector3;

class CPAnalyzer: public Analyzer{

 public:

  enum tauDecayModes {kElectron, kMuon, 
		      kOneProng0pi0, kOneProng1pi0, kOneProng2pi0, kOneProng3pi0,
		      kThreeProng0pi0, kThreeProng1pi0,
		      kOther, kUndefined};
  
  CPAnalyzer(const std::string & aName);

  virtual ~CPAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *tree){;};

  virtual void addCutHistos(TList *aList){;};

  const std::string & name(){return myName_;};
  bool filter() const{ return filterEvent_;};

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  ///Calculate variabmes for the pi+pi- decay mode
  ///sysType enumarates different reality scenarios
  //(w/o cuts, PV smearing etc.)
  void fillAngles(const EventProxyCPNtuple & myEvent,
		  const std::string & sysType);
  

  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter);

  //Calculate the transverse impace parameter using primary vertex (pv)
  //decay (secondary) vertex (sc), and decaying particle momentum.
  TVector3 impactParameter(const TVector3& pv, const TVector3& sv, 
			   const TLorentzVector& momentum);

  ///Return string encoding di-tau decay mode
  std::string getDecayName(int decModeMinus, int decModePlus);

  ///Return mother resonance name.
  std::string getMotherName(int bosonId);

  ///Check if the decMode points to single prong hadronic tau decay
  bool isOneProng(int decMode);

  ///Check if the decMode points to leptonic tau decay
  bool isLepton(int decMode);
 
  std::string myName_;

  ///Random numbeer gnerator
  TRandom3 aRndm;
  
  ///Histograms storage
  CPHistograms *myHistos_;
  
  //should this CPAnalyzer be able to filter events
  bool filterEvent_;
 
};

#endif
