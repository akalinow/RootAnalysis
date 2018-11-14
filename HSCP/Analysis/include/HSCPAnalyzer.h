#ifndef RootAnalysis_HSCPAnalyzer_H
#define RootAnalysis_HSCPAnalyzer_H

#include <string>
#include <vector>
#include <map>

#include "EventProxyBase.h"
#include "EventProxyHSCP.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"
class HSCPHistograms;

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;

class HSCPAnalyzer: public Analyzer{

 public:

  HSCPAnalyzer(const std::string & aName);

  virtual ~HSCPAnalyzer();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger);

  virtual void finalize();

  virtual void clear(){;};

  Analyzer* clone() const;

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const HSCPEvent & myEvent);

 protected:

  void setHistos(HSCPHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  HSCPHistograms *myHistos_;

};

#endif
