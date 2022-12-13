#ifndef GMTAnalyzer_H
#define GMTAnalyzer_H

#include <string>
#include <map>

#include "GMTHistograms.h"
#include "Analyzer.h"

#include "EventObj.h"
#include "GenObjColl.h"
#include "L1ObjColl.h"
#include "L1PhaseIIObjColl.h"
#include "L1PhaseIIObj.h"


#include "TVector3.h"

class TriggerHistograms;
class OMTFHit;

class GMTAnalyzer:public Analyzer{

 public:

  GMTAnalyzer(const std::string & aName);

  ~GMTAnalyzer();

  void initialize(TDirectory* aDir,
		  pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  Analyzer* clone() const;

  void setHistos(GMTHistograms *histos) { myHistos_ = histos;};

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

  double calibratedPt(const std::string & sysType, const double & ptRaw);

  TFile *calibrationFile{0};
  TH1F *hCalibationHisto{0};
  TF1 *calibrationFunc{0};
  ///Histograms for this analysis
  GMTHistograms *myHistos_;

  const EventObj  *myEventId;
  const GenObjColl *myGenObjColl;
  const L1ObjColl  *myL1ObjColl;
  const L1PhaseIIObjColl  *myL1PhaseIIObjColl;

  GenObj myGenObj;
  std::vector<OMTFHit> myHits;

  TH3F *hGoldenPatterns;
  TH1D *hPtProfile;

  TVector3 genMuMom;

  std::map<unsigned long, int> quality_index_map;
  
};

#endif
