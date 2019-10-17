#ifndef RootAnalysis_TauLFVAnalyzer_H
#define RootAnalysis_TauLFVAnalyzer_H

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
class TauLFVHistograms;

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;

class TauLFVAnalyzer: public Analyzer{

 public:

  TauLFVAnalyzer(const std::string & aName): Analyzer(aName){};

  virtual ~TauLFVAnalyzer();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual void finalize();

  virtual void clear();

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

  ///Find reconstructed objectd interesting for this analysis,
  ///for example the three muons.
  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const std::string & hNameSuffix, float eventWeight);

  ///Fill the bookkeeping histograms with number of events
  ///after the preselection.
  void getPreselectionEff(const EventProxyHTT & myEventProxy);
			 
  
 protected:

  pat::strbitset *mySelections_;

  virtual void setHistos(TauLFVHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  TauLFVHistograms *myHistos_;

  //should this TauLFVAnalyzer be able to filter events
  bool filterEvent_;

  ///Reconstructed objects selected for given event.
  HTTEvent * aEvent;
  std::string sampleName;

  HTTParticle aMuon1, aMuon2, aMuon3;

};

#endif
