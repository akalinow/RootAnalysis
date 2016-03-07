
#ifndef HWvsEMULAnalyzer_H
#define HWvsEMULAnalyzer_H

#include <string>
#include <map>

#include "HWvsEMULHistograms.h"
#include "Analyzer.h"
#include "EventData.h"

class TriggerHistograms;

class HWvsEMULAnalyzer:public Analyzer{

 public:

  HWvsEMULAnalyzer(const std::string & aName);

  ~HWvsEMULAnalyzer();

  void initialize(TFileDirectory& aDir,
		  pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  void addBranch(TTree *tree);

  Analyzer* clone() const;

  void setHistos(HWvsEMULHistograms *histos) { myHistos_ = histos;};

 private:

  void registerCuts() {;};
  bool checkSelections(const std::string & type) {;};
  void clear();

  //Class stored in the TTree
  EventData       *theEvent;

  ///Histograms for this analysis
  HWvsEMULHistograms *myHistos_;

  ///Fill histograms with various distributions
  void fillCandidates();
  
  ///Variables stored in the TTree
  std::map<std::string,float> treeVariables_;

  std::vector<L1Obj> emulCandidates, hwCandidates;

  float eventWeight_;
  
};

#endif
