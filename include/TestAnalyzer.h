#ifndef FWLiteTriggerAnalyzer_H
#define FWLiteTriggerAnalyzer_H

// -*- C++ -*-
//
//
// Original Author:  Artur Kalinowski
//         Created:  Fri Sep 18 10:37:10 CEST 2009
//
//

#include <string>
#include <map>

#include "Analyzer.h"

class TriggerHistograms;

class TestAnalyzer:public Analyzer{

 public:

  TestAnalyzer(const std::string & aName);

  virtual ~FWLiteTriggerAnalyzer();

  virtual void initialize(const TFileDirectory&,
			  pat::strbitset *aSelections);
  
  virtual bool analyze(const EventBase& iEvent);

  virtual void finalize();

  virtual void addBranch(TTree *tree);

 private:

  void registerCuts();
  bool checkSelections(const std::string & type);
  void clear();


  ///Variables stored in the TTree
  std::map<std::string,float> treeVariables_;

  //TriggerHistograms *myHistos_;
  float eventWeight_;
  
};

#endif
