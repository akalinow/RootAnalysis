#ifndef OMTFAnalyzer_H
#define OMTFAnalyzer_H

#include <string>
#include <map>

#include "OMTFHistograms.h"
#include "Analyzer.h"

#include "L1ObjColl.h"

class TriggerHistograms;

class OMTFAnalyzer:public Analyzer{

 public:

  OMTFAnalyzer(const std::string & aName);

  ~OMTFAnalyzer();

  void initialize(TDirectory* aDir,
		  pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  void addBranch(TTree *tree);

  Analyzer* clone() const;

  void setHistos(OMTFHistograms *histos) { myHistos_ = histos;};

 private:

  void fillTurnOnCurve(const int & ptCut, const std::string & sysType,
		               const std::string & selType);

  void fillRateHisto(const std::string & sysType,
  			         const std::string & selType);

  void fillGhostHisto(const std::string & sysType,
		      const std::string & selType);
  
  bool passQuality(std::vector<L1Obj> * myL1Coll,
		   const std::string & sysType, 
		   unsigned int iCand);

  void registerCuts();
  bool checkSelections(const std::string & type);
  void clear();

  ///Histograms for this analysis
  OMTFHistograms *myHistos_;

  ///Variables stored in the TTree
  std::map<std::string,float> treeVariables_;

  //TriggerHistograms *myHistos_;
  float eventWeight_;

  std::map<int, float> tmpMap;
  
};

#endif
