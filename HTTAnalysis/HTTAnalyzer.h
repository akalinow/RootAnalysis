#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>
#include <vector>
#include <map>

#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"
#include "ChannelSpecifics.h"
#include "AnalysisEnums.h"
#include "ObjectMessenger.h"
class HTTHistograms;

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;

class HTTAnalyzer: public Analyzer{

  friend class ChannelSpecifics;
  friend class MuTauSpecifics;
  friend class TauTauSpecifics;
  friend class MuMuSpecifics;

 public:

  HTTAnalyzer(const std::string & aName, const std::string & aDecayMode = "None");

  virtual ~HTTAnalyzer();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger);

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *);

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  ///Check it the event passes given category selections.
  bool passCategory(unsigned int iCategory);

  ///Return human readable sample name (Data, WJets, etc).
  std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it.
  ///Method used when sample coding in TTree is not present.
  ///In this case a ROOT file name is used to decode the sample type.
  std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy);

  ///Return sample name for DY. Name encoded jet bin, and decay mode.
  std::string getDYSampleName(const EventProxyHTT & myEventProxy);

  //Return name sample name suffix for different particles matched to reconstructed tau
  std::string getMatchingName(const EventProxyHTT & myEventProxy);

  ///Return pileup reweighting weight.
  ///Weight is calculatedon fly using the ration of nPU
  ///histograms for data and analyased sample.
  float getPUWeight(const EventProxyHTT & myEventProxy);

  ///Return event weight for systematic effects
  ///implemented by a global event weight.
  float getSystWeight(const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);

  ///Fill pulls between generator and various reco vertices.
  bool fillVertices(const std::string & sysType, float eventWeight);

  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const std::string & hNameSuffix, float eventWeight,
			 const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);


  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for reconstructed
  ///mu+tau_h mode
  void fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight,
  const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);

  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for
  ///generator level taus for all decay modes.
  void fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);

  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter,
					    bool sgn=true);

  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
					    float deltaR);

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  void getPreselectionEff(const EventProxyHTT & myEventProxy);

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Parts of code specific to give decay channel.
  ///In particular category and object selection.
  ChannelSpecifics *myChannelSpecifics;

  ///Histograms storage.
  HTTHistograms *myHistos_;

  ///ROOT file with PU histogram
  TFile *puDataFile_, *puMCFile_;

  ///ROOT file containing current TTree
  TFile *ntupleFile_;

  ///Histogram with event counts filled during preselection step.
  TH1F *hStatsFromFile;

  ///Vector of PU histograms for MC samples
  std::vector<TH1F*> hPUVec_;

  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;

  ///Map from file name to sample name.
  std::map<std::string, std::string> fileName2sampleName;

  ///Reconstructed objects selected for given event.
  HTTEvent aEvent;
  HTTPair aPair;
  std::string sampleName;

  HTTParticle aLeg2, aLeg1, aMET;
  HTTParticle aGenLeg1, aGenLeg2;
  HTTParticle aJet1, aJet2, aBJet1;
  std::vector<HTTParticle> aSeparatedJets;
  int nJets30;
  int nJetsInGap30;
  int nBJets;
  std::vector<bool> categoryDecisions;
  unsigned int myNumberOfCategories;

  //cut on nPCA
  float nPCAMin_;

};

#endif
