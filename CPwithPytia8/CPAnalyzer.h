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

#include "Analyzer.h"
#include "CPHistograms.h"


class CPAnalyzer: public Analyzer{

 public:

  CPAnalyzer(const std::string & aName);

  virtual ~CPAnalyzer();
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *tree){;};

  virtual void addCutHistos(TList *aList){;};

  const std::string & name(){return myName_;};
  bool filter() const{ return filterEvent_;};

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  std::string myName_;

  ///Histograms storage
  CPHistograms *myHistos_;
  
  //should this CPAnalyzer be able to filter events
  bool filterEvent_;
 
};

#endif
