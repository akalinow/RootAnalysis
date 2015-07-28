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

  ~OTFAnalyzer();

  void initialize(TFileDirectory& aDir,
		  pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  void addBranch(TTree *tree);

  Analyzer* clone() const;

  void setHistos(OTFHistograms *histos) { myHistos_ = histos;};

 private:

  void fillTurnOnCurve(const int & ptCut, const std::string & sysType,
		               const std::string & selType);

  void fillRateHisto(const std::string & sysType,
  			         const std::string & selType);

  void fillGhostHisto(const std::string & sysType,
		      const std::string & selType);
  
  bool passQuality(std::vector<L1Obj> * myL1Coll,
		   const std::string & sysType, 
		   int iCand);

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

  std::map<int, float> tmpMap;
  
};

#endif
