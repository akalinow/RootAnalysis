#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>
#include <vector>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TFileDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"

class HTTHistograms;

class TH1F;
class TLorentzVector;

class HTTAnalyzer: public Analyzer{

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
  
  HTTAnalyzer(const std::string & aName);

  virtual ~HTTAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *);

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  ///Check it tau decay modes (GEN and RECO) match selected (hardcoded)
  ///decay mode.
  std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it
  static std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return pileup reweighting weight.
  ///Weight is calculatedon fly using the ration of nPU 
  ///histograms for data and analyased sample.
  float getPUWeight(const EventProxyHTT & myEventProxy);

  ///Fill pulls between generator and various reco vertices.
  bool fillVertices(const std::string & sysType);
  
  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const std::string & hNameSuffix, float eventWeight);
  
  ///Fill histograms with cos(phi), where phi is the decay 
  ///between tau decay planes. Method used for reconstructed
  ///mu+tau_h mode
  void fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);

  ///Fill histograms with cos(phi), where phi is the decay 
  ///between tau decay planes. Method used for 
  ///generator level taus for all decay modes.
  void fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);
  
  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter);

  ///Return string encoding di-tau decay mode.
  ///The event can belong to more than one category
  std::vector<std::string> getTauDecayName(int decModeMinus, int decModePlus);

  ///Check if the decMode points to single prong hadronic tau decay
  bool isOneProng(int decMode);
  
  ///Check if the decMode points to leptonic tau decay
  bool isLepton(int decMode);

  ///Get jets separated by deltaR from tau an muon.
  std::vector<Wjet> getSeparatedJets(const EventProxyHTT & myEventProxy, const Wtau & aTau,
				     const Wmu & aMuon, float deltaR);
  
 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  void getPreselectionEff(const EventProxyHTT & myEventProxy);

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  HTTHistograms *myHistos_;

  ///ROOT file with PU histogram
  TFile *puDataFile_, *puMCFile_;

  ///ROOT file containing current TTree
  TFile *ntupleFile_;
 
  ///Vector of PU histograms for MC samples
  std::vector<TH1F*> hPUVec_;
 
  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;

  ///Reconstructed objects selected for given event.
  Wevent aEvent;
  Wpair aPair;
  Wtau aTau;
  Wtau aGenNegativeTau;
  Wtau aGenPositiveTau;
  Wmu aMuon;
  Wmet aMET;
  std::vector<Wjet> aSeparatedJets;
  Wjet aJet;
  int nJets30;

  ///Variables stored in a TTree
  Double_t muonPt;
 
};

#endif
