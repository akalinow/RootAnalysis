#ifndef PACAnalyzer_H
#define PACAnalyzer_H

#include <string>
#include <map>

#include "PACHistograms.h"
#include "Analyzer.h"
#include "EventData.h"

class TriggerHistograms;

class PACAnalyzer:public Analyzer{

 public:

  PACAnalyzer(const std::string & aName);

  virtual ~PACAnalyzer();

  virtual void initialize(TFileDirectory& aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual void finalize();

  virtual void addBranch(TTree *tree);

 private:

  void fillTurnOnCurve(const int & ptCut, const std::string & sysType,
		               const std::string & selType);

  void fillRateHisto(const std::string & sysType,
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
  PACHistograms *myHistos_;

  ///Variables stored in the TTree
  std::map<std::string,float> treeVariables_;

  //TriggerHistograms *myHistos_;
  float eventWeight_;

  std::map<int, float> tmpMap;
  
};

#endif
