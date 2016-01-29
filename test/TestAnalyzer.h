#ifndef RootAnalysis_TestAnalyzer_H
#define RootAnalysis_TestAnalyzer_H

#include <string>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"

#include "strbitset.h"
#include "TFileDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"

class EventProxyTest;
class TestHistograms;

class TestAnalyzer: public Analyzer{

 public:
  
  TestAnalyzer(const std::string & aName);

  virtual ~TestAnalyzer();
  
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

  void setHistos(TestHistograms *histos) { myHistos_ = histos;};

  ///Histograms storage.
  TestHistograms *myHistos_;
  
  //should this TestAnalyzer be able to filter events
  bool filterEvent_;

  std::string tmpName;
 
};

#endif
