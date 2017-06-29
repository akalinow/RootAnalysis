#ifndef RootAnalysis_ChannelSpecifics_H
#define RootAnalysis_ChannelSpecifics_H

#include <string>
#include "HTTEvent.h"
#include "AnalysisEnums.h"
#include "TH1F.h"

class HTTAnalyzer;
class EventProxyHTT;
class TH1F;
class TH2F;
class TH3F;
class TRandom3;
class TFile;

class BTagCalibration;
class BTagCalibrationReader;

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

float getLeptonCorrection(float eta, float pt, float iso, HTTAnalysis::hadronicTauDecayModes tauDecayMode, bool useTauTrigger, int mc_match, bool useXTrigger = false);

float getDYReweight(const std::string & categoryName, const HTTAnalysis::sysEffects & aSystEffect = HTTAnalysis::NOMINAL);

virtual bool promoteBJet(const HTTParticle &jet,
			 const HTTAnalysis::sysEffects &aSystEffect=HTTAnalysis::NOMINAL,
			 std::string correctionType="central");

virtual std::string getDecayModeName() const {
        return decayModeName;
}

const std::vector<const HTTAnalysis::eventCategory*> & getCategoryRejester() const {
        return categoryRejester;
}

float getTauIDSF(float eta, int mc_match);

protected:

virtual void defineCategories();

void initializeLeptonCorrections();

void initializeBTagCorrections();

HTTAnalyzer *myAnalyzer;

///Histograms with lepton corrections
TH2F *h2DMuonIdCorrections;
TH3F *h3DMuonIsoCorrections, *h3DMuonTrgCorrections, *h3DMuonXTrgCorrections;
TH1F *h1DMuonTrkCorrections;
TH3F *h3DTauCorrections;
TH2F *h2DTauTrgGenuineCorrections, *h2DTauTrgFakeCorrections;
TH2F *h2DTauXTrgGenuineCorrections, *h2DTauXTrgFakeCorrections;

//For btag calibration
BTagCalibration *calib;
BTagCalibrationReader *reader;
TFile *btagEffFile_;
TH2F *btag_eff_b_, *btag_eff_c_, *btag_eff_oth_;
TRandom3 *rand_;

std::string decayModeName = "None";

std::vector<const HTTAnalysis::eventCategory*> categoryRejester;
HTTAnalysis::eventCategory *jet0, *boosted, *vbf;
HTTAnalysis::eventCategory *antiIso_jet0, *antiIso_boosted, *antiIso_vbf;
HTTAnalysis::eventCategory *mu_pi, *mu_rho, *pi_pi, *pi_rho, *rho_rho;
HTTAnalysis::eventCategory *inclusive;
HTTAnalysis::eventCategory *antiIso_inclusive;
HTTAnalysis::eventCategory *btag, *nobtag;
HTTAnalysis::eventCategory *antiIso_btag, *antiIso_nobtag;

TH1F *tauID_FRSF_mu, *tauID_FRSF_ele;
Float_t bins_mu_[6] = {0,0.4,0.8,1.2,1.7,2.3}, bins_ele_[4] = {0,1.46,1.558,10};

};


#endif
