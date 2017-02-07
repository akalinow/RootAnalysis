#ifndef ChannelSpecifics_HTTAnalyzer_H
#define ChannelSpecifics_HTTAnalyzer_H

#include <string>
#include "HTTEvent.h"

class EventProxyHTT;
class HTTAnalyzer;

class ChannelSpecifics {

public:
ChannelSpecifics(HTTAnalyzer *aAnalyzer){myAnalyzer = aAnalyzer;};


virtual void setAnalysisObjects(const EventProxyHTT & myEventProxy) = 0;

///Check it tau decay modes (GEN and RECO) match selected (hardcoded)
///decay mode.
virtual std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy) = 0;

///Check it the event passes category selections with given systematic effect.
virtual void testAllCategories(const sysEffects::sysEffectsEnum & aSystEffect) = 0;

std::string getDecayModeName() const {return decayModeName;}

protected:

HTTAnalyzer *myAnalyzer;

std::string decayModeName = "None";
};


#endif
