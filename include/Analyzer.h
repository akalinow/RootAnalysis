#ifndef RootAnalysis_Analyzer_H
#define RootAnalysis_Analyzer_H

#include <string>

#include "ObjectMessenger.h"
#include "EventBase.h"

#include "strbitset.h"
#include "TFileDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"


class Analyzer{

 public:

  Analyzer(const std::string & aName);

  virtual ~Analyzer(){;}
  
  ///Initialize the analyzer
  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections){

    //default is to do whatever the analzer does
    filterEvent_ = false;
};
  
  
  virtual bool analyze(const EventBase& iEvent) = 0;

  virtual bool analyze(const EventBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize(){;};

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
  //should this Analyzer be able to filter events
  bool filterEvent_;
 
};

#endif
