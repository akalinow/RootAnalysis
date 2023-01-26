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

  void initialize(TDirectory* aDir, pat::strbitset *aSelections);

  bool analyze(const EventProxyBase& iEvent);

  void finalize();

  Analyzer* clone() const;

  void setHistos(GMTHistograms *histos) { myHistos_ = histos;};

 private:

  void fillHistosForGenMuon(const GenObj & aGenObj);

  void fillTurnOnCurve(const GenObj & aGenObj,
                      const int & ptCut, const std::string & sysType,
		                  const std::string & selType);

  void fillRateHisto(const GenObj & aGenObj,
                    const std::string & sysType,
		                const std::string & selType);

  bool passQuality(const L1PhaseIIObj & aL1Cand,
		              const std::string & sysType,
		              const std::string & selType = "");
 
  ///Histograms for this analysis
  GMTHistograms *myHistos_;

  const EventObj  *myEventId;
  const GenObjColl *myGenObjColl;
  const L1ObjColl  *myL1ObjColl;
  const L1PhaseIIObjColl  *myL1PhaseIIObjColl; 
};

#endif
