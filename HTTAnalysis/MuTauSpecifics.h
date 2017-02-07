#ifndef MuTauSpecifics_HTTAnalyzer_H
#define MuTauSpecifics_HTTAnalyzer_H

#include "ChannelSpecifics.h"

class HTTAnalyzer;
class EventProxyHTT;

class MuTauSpecifics: public ChannelSpecifics{

public:
  MuTauSpecifics(HTTAnalyzer *aAnalyzer);

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy);

  void testAllCategories(const sysEffects::sysEffectsEnum & aSystEffect);

  std::string decayModeName = "MuTau";

};
#endif
