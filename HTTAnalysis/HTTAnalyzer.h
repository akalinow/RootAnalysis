#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>

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

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  HTTHistograms *myHistos_;
  
  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;
 
};

#endif
