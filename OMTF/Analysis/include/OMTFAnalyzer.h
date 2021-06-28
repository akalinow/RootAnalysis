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
class OMTFHit;

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
		   const std::string & sysType,
		   const std::string & selType = "");

  void fixQualityHistos();

  void fillBendingHistos(const std::string & sysType);

  std::pair<double, double> getPtProfile();

  ///Histograms for this analysis
  OMTFHistograms *myHistos_;

  const EventObj  *myEventId;
  const GenObjColl *myGenObjColl;
  const L1ObjColl  *myL1ObjColl;
  GenObj myGenObj;
  std::vector<OMTFHit> myHits;

  TH3F *hGoldenPatterns;
  TH1D *hPtProfile;

  TVector3 genMuMom;

  std::map<unsigned long, int> quality_index_map;
  
};

#endif
