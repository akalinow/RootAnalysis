#ifndef RootAnalysis_SVfitAnalyzer_H
#define RootAnalysis_SVfitAnalyzer_H

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

#include "Analyzer.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"

class TLorentzVector;

class SVfitAnalyzer: public Analyzer{

  public:

  SVfitAnalyzer(const std::string & aName, const std::string & aDecayMode = "None");

  virtual ~SVfitAnalyzer();

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger) override;

  virtual bool analyze(const EventProxyBase& iEvent) override {return analyze(iEvent, nullptr); }

  Analyzer* clone() const override;

  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
      float deltaR);

  private:

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  TLorentzVector computeMTT(const std::string & algoName);

  TLorentzVector runSVfitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons);

  TLorentzVector runFastMTTAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons);

  std::tuple<double, double> getTauMomentum(const TLorentzVector & visP4, double cosGJ);

  //decayMode
  const std::string decayMode_;

  ///Reconstructed objects selected for given event.
  HTTEvent event_;
  HTTPair pair_;
  std::string sampleName_;

  HTTParticle leg2_;
  HTTParticle leg1_;
  HTTParticle genLeg1_;
  HTTParticle genLeg2_;

  HTTParticle jet1_;
  HTTParticle jet2_;
  std::vector<HTTParticle> separatedJets_;

  HTTParticle MET_;
  TMatrixD covMET_;

  float genSumM_;
  float higgsMassTrans_;

  ClassicSVfit SVfitAlgo_;
  FastMTT fastMTTAlgo_;

};

#endif
