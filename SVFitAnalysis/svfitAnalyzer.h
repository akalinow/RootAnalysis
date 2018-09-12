
#ifndef RootAnalysis_svfitAnalyzer_H
#define RootAnalysis_svfitAnalyzer_H

#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"
#include "TRandom3.h"

#include "Analyzer.h"
#include "ChannelSpecifics.h"
#include "AnalysisEnums.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"

class TLorentzVector;

class svfitAnalyzer: public Analyzer{

  friend class ChannelSpecifics;
  friend class MuTauSpecifics;
  friend class TauTauSpecifics;
  friend class MuMuSpecifics;

 public:

  svfitAnalyzer(const std::string & aName, const std::string & aDecayMode = "None");

  virtual ~svfitAnalyzer();

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger) override;

  virtual bool analyze(const EventProxyBase& iEvent) override {return analyze(iEvent, nullptr); }

  Analyzer* clone() const override;

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

  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
    float deltaR);

 private:

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  TLorentzVector computeMTT(const std::string & algoName);
  
  TLorentzVector runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons);

  TLorentzVector runFastMTTAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons);

  std::tuple<double, double> getTauMomentum(const TLorentzVector & visP4, double cosGJ);

  ///Parts of code specific to give decay channel.
  ///In particular category and object selection.
  ChannelSpecifics *channelSpecifics_;

	//decayMode
	const std::string decayMode_;

  ///Reconstructed objects selected for given event.
  HTTEvent event_;
  HTTPair pair_;
  std::string sampleName_;

  HTTParticle leg2_;
	HTTParticle leg1_;

	HTTParticle MET_;
	TMatrixD covMET_;

  HTTParticle genLeg1_;
	HTTParticle genLeg2_;
  float genSumM_;
  float higgsMassTrans_;
  HTTParticle jet1_;
	HTTParticle jet2_;
  std::vector<HTTParticle> separatedJets_;
  int nJets30_; //TODO remove their usage in ChannelSpecifics
  int nJetsInGap30_;

  std::vector<bool> categoryDecisions_;
  unsigned int numberOfCategories_;

  ClassicSVfit svFitAlgo_;
  FastMTT fastMTTAlgo_;
};

#endif
