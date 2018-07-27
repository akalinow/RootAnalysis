#ifndef RootAnalysis_MuMuSpecifics_H
#define RootAnalysis_MuMuSpecifics_H

#include "ChannelSpecifics.h"

class svfitAnalyzer;
class EventProxyHTT;

class MuMuSpecifics: public ChannelSpecifics{

public:
  MuMuSpecifics(svfitAnalyzer *aAnalyzer);

  ~MuMuSpecifics() {}

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy);

  void testAllCategories(const HTTAnalysis::sysEffects & aSystEffect);

float getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect);
float getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect);

virtual std::string getDecayModeName() const {return decayModeName;}

private:

  std::string decayModeName = "MuMu";

};
#endif
