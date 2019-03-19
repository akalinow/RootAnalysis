#ifndef OMTFAnalyzer_H
#define OMTFAnalyzer_H

#include <string>
#include <map>

#include "OMTFHistograms.h"
#include "Analyzer.h"

#include "EventObj.h"
#include "GenObjColl.h"
#include "L1ObjColl.h"

#include "TVector3.h"

class TriggerHistograms;

class OMTFAnalyzer:public Analyzer{

 public:

  OMTFAnalyzer(const std::string & aName);

  ~OMTFAnalyzer();

  void initialize(TDirectory* aDir,
		  pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  Analyzer* clone() const;

  void setHistos(OMTFHistograms *histos) { myHistos_ = histos;};

 private:

  void fillHistosForGenMuon();

  void fillTurnOnCurve(const int & ptCut, const std::string & sysType,
		               const std::string & selType);

  void fillRateHisto(const std::string & sysType,
		     const std::string & selType);

  bool passQuality(const L1Obj & aL1Cand,
		   const std::string & sysType);

  ///Histograms for this analysis
  OMTFHistograms *myHistos_;

  const EventObj  *myEventId;
  const GenObjColl *myGenObjColl;
  const L1ObjColl  *myL1ObjColl;

  TVector3 genMuMom;

  std::map<int, int> hitsPatterns;
  
};

#endif
