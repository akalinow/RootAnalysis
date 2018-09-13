#ifndef RootAnalysis_Analyzer_H
#define RootAnalysis_Analyzer_H

#include <string>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"


class Analyzer{

 public:
  Analyzer(){}

  Analyzer(const std::string & aName);

  virtual ~Analyzer(){;}
  
  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections){

    //default is to do whatever the analyzer does
    filterEvent_ = false;
};
  

  virtual Analyzer * clone() const{ return 0;}
  
  virtual bool analyze(const EventProxyBase& iEvent) = 0;

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize(){;};

  virtual void clear(){;};

  virtual void addBranch(TTree *tree){;};

  virtual void addCutHistos(TList *aList){;};

  const std::string & name() const{return myName_;};
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
