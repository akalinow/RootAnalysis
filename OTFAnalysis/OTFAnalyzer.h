#ifndef OTFAnalyzer_H
#define OTFAnalyzer_H

#include <string>
#include <map>

#include "OTFHistograms.h"
#include "Analyzer.h"
#include "EventData.h"

class TriggerHistograms;

class OTFAnalyzer:public Analyzer{

 public:

  OTFAnalyzer(const std::string & aName);

  virtual ~OTFAnalyzer();

  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual void finalize();

  virtual void addBranch(TTree *tree);

 private:

  void fillTurnOnCurve(int & ptCut, const std::string & sysType);

  void registerCuts();
  bool checkSelections(const std::string & type);
  void clear();

  //Class stored in the TTree
  EventData       *theEvent;

  ///Histograms for this analysis
  OTFHistograms *myHistos_;

  ///Variables stored in the TTree
  std::map<std::string,float> treeVariables_;

  //TriggerHistograms *myHistos_;
  float eventWeight_;
  
};

#endif
