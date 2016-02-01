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
class HTTEvent;
class DiTauData;

class CPAnalyzer: public Analyzer{

 public:
  
  ///Copy from DataFormats/TauReco/interface/PFTauDecayMode.h
  enum hadronicTauDecayModes 
  {
    tauDecay1ChargedPion0PiZero,
    tauDecay1ChargedPion1PiZero,  // rho (770 MeV) mediated)
    tauDecay1ChargedPion2PiZero,  // a1  (1.2 GeV) mediated
    tauDecay1ChargedPion3PiZero,  // contaminated or unmerged photo
    tauDecay1ChargedPion4PiZero,  // contaminated or unmerged photo
    tauDecay2ChargedPion0PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion1PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion2PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion3PiZero,  // extra track or un-recod track
    tauDecay2ChargedPion4PiZero,  // extra track or un-recod track
    tauDecay3ChargedPion0PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion1PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion2PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion3PiZero,  // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion4PiZero,  // a1  (1.2 GeV) mediated
    tauDecaysElectron,
    tauDecayMuon,
    tauDecayOther                 // catch-all
  };
  
  CPAnalyzer(const std::string & aName);

  virtual ~CPAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual Analyzer* clone() const;

  virtual void addBranch(TTree *tree){;};

  virtual void addCutHistos(TList *aList){;};

  bool filter() const{ return filterEvent_;};

  void setHistos(CPHistograms *histos) { myHistos_ = histos;};

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  //Check if the event passes analysis like selection cuts.
  bool analysisSelection(const DiTauData *aEvent);

  ///Calculate variabmes for the pi+pi- decay mode
  ///sysType enumarates different reality scenarios
  //(w/o cuts, PV smearing etc.)
  void fillAngles(const DiTauData *aEvent,
		  const std::string & sysType);

  ///Fill pulls between generator and various reco vertices.
  bool fillVertices(const DiTauData* aEventGen,
		    const DiTauData* aEventReco,
		    const std::string & sysType);
  

  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter);

  //Calculate the transverse impace parameter using primary vertex (pv)
  //decay (secondary) vertex (sc), and decaying particle momentum.
  TVector3 impactParameter(const TVector3& pv, const TVector3& sv, 
			   const TLorentzVector& momentum);

  ///Return string encoding di-tau decay mode.
  ///The event can belong to more than one category
  std::vector<std::string> getDecayName(int decModeMinus, int decModePlus);

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

  HTTEvent *myEvent;
  
  //should this CPAnalyzer be able to filter events
  bool filterEvent_;
 
};

#endif
