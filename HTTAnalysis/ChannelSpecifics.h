#ifndef RootAnalysis_ChannelSpecifics_H
#define RootAnalysis_ChannelSpecifics_H

#include <string>
#include "HTTEvent.h"
#include "AnalysisEnums.h"

class HTTAnalyzer;
class EventProxyHTT;
class TH2F;
class TH3F;

class ChannelSpecifics {

public:

ChannelSpecifics(HTTAnalyzer *aAnalyzer);

~ChannelSpecifics();

virtual void setAnalysisObjects(const EventProxyHTT & myEventProxy) = 0;

///Check it tau decay modes (GEN and RECO) match selected (hardcoded)
///decay mode.
virtual std::pair<bool, bool> checkTauDecayMode(const EventProxyHTT & myEventProxy) = 0;

///Check it the event passes category selections with given systematic effect.
virtual void testAllCategories(const HTTAnalysis::sysEffects & aSystEffect) = 0;

///Return cumulative MC corrections for the leg1
virtual float getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect) = 0;

///Return cumulative MC corrections for the leg2
virtual float getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect) = 0;

float getLeptonCorrection(float eta, float pt, HTTAnalysis::hadronicTauDecayModes tauDecayMode, bool useTauTrigger);

virtual std::string getDecayModeName() const {
        return decayModeName;
}

const std::vector<const HTTAnalysis::eventCategory*> & getCategoryRejester() const {return categoryRejester;}

protected:

virtual void defineCategories();

///Convert RooRealVar functions to histograms
void initializeCorrections();

HTTAnalyzer *myAnalyzer;

///Histograms with lepton corrections
TH2F *h2DMuonIdCorrections, *h2DMuonIsoCorrections, *h2DMuonTrgCorrections;
TH3F *h3DTauCorrections, *h3DTauTrgOSCorrections, *h3DTauTrgSSCorrections;

std::string decayModeName = "None";

std::vector<const HTTAnalysis::eventCategory*> categoryRejester;
HTTAnalysis::eventCategory *jet0, *boosted, *vbf;
HTTAnalysis::eventCategory *antiIso_jet0, *antiIso_boosted, *antiIso_vbf;
HTTAnalysis::eventCategory *mu_pi, *mu_rho, *pi_pi, *pi_rho, *rho_rho;

};


#endif
