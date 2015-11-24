#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>
#include <vector>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"

#include "strbitset.h"
#include "TFileDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"

class EventProxyHTT;
class HTTHistograms;
class TH1F;

class HTTAnalyzer: public Analyzer{

 public:
  
  HTTAnalyzer(const std::string & aName);

  virtual ~HTTAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  bool filter() const{ return filterEvent_;};

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it
  static std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return pileup reweighting weight.
  ///Weight is calculatedon fly using the ration of nPU 
  ///histograms for data and analyased sample.
  float getPUWeight(const EventProxyHTT & myEventProxy);

  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  ///Histograms storage.
  HTTHistograms *myHistos_;

  ///ROOT file with PU histogram
  TFile *puFile_;
 
  ///Vector of PU histograms for MC samples
  std::vector<TH1F*> hPUVec_;
 
  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;
 
};

#endif
