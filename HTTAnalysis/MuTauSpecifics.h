#ifndef RootAnalysis_MuTauSpecifics_H
#define RootAnalysis_MuTauSpecifics_H

#include "ChannelSpecifics.h"

class HTTAnalyzer;
class EventProxyHTT;

class MuTauSpecifics: public ChannelSpecifics{

public:
  MuTauSpecifics(HTTAnalyzer *aAnalyzer);

  ~MuTauSpecifics() {}

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy);

  void testAllCategories(const HTTAnalysis::sysEffects & aSystEffect);

float getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect);
float getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect);

virtual std::string getDecayModeName() const {return decayModeName;}

//float getTauIDSF(float eta, int mc_match);

private:

  std::string decayModeName = "MuTau";

};
#endif
