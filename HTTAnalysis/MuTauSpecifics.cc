#include "MuTauSpecifics.h"
#include "HTTAnalyzer.h"
#include "EventProxyHTT.h"
#include "Tools.h"


MuTauSpecifics::MuTauSpecifics(HTTAnalyzer * aAnalyzer) : ChannelSpecifics(aAnalyzer){

        decayModeName = "MuTau";

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MuTauSpecifics::setAnalysisObjects(const EventProxyHTT & myEventProxy) {

        myAnalyzer->aLeg1 = myAnalyzer->aPair.getMuon();
        myAnalyzer->aLeg2 = myAnalyzer->aPair.getTau();

        myAnalyzer->aGenLeg1 = HTTParticle();
        myAnalyzer->aGenLeg2 = HTTParticle();

        if(myEventProxy.genLeptons &&
           myEventProxy.genLeptons->size()) {
                HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
                if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->aGenLeg1 = aGenTau;
                else if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) myAnalyzer->aGenLeg2 = aGenTau;
                if(myEventProxy.genLeptons->size()>1) {
                        aGenTau =  myEventProxy.genLeptons->at(1);
                        if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->aGenLeg1 = aGenTau;
                        else if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) myAnalyzer->aGenLeg2 = aGenTau;
                }
        }
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<bool, bool> MuTauSpecifics::checkTauDecayMode(const EventProxyHTT & myEventProxy){

        bool goodGenDecayMode = false;
        bool goodRecoDecayMode = false;

        std::vector<std::string> decayNamesGen = HTTAnalysis::getTauDecayName(myAnalyzer->aGenLeg2.getProperty(PropertyEnum::decayMode),
                                                                             myAnalyzer->aGenLeg1.getProperty(PropertyEnum::decayMode));
        std::vector<std::string> decayNamesReco = HTTAnalysis::getTauDecayName(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode),HTTAnalysis::tauDecayMuon);

        for(auto it: decayNamesGen) if(it.find("Lepton1Prong")!=std::string::npos) goodGenDecayMode = true;
        for(auto it: decayNamesReco) if(it.find("Lepton1Prong")!=std::string::npos) goodRecoDecayMode = true;

        return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MuTauSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it: myAnalyzer->categoryDecisions) it = false;

        int tauIDmask = 0;
        for(unsigned int iBit=0; iBit<myAnalyzer->aEvent.ntauIds; iBit++) {
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }
        bool tauID = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool muonKinematics = myAnalyzer->aLeg1.getP4().Pt()>24 && fabs(myAnalyzer->aLeg1.getP4().Eta())<2.1;

        bool trigger = myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
                       myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
                       myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
                       myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1);

        if(myAnalyzer->sampleName!="Data") trigger = true; //MC trigger included in muon SF
        if(!muonKinematics || !tauID || !trigger) return;

        myAnalyzer->nJets30 = 0;
        for(auto itJet: myAnalyzer->aSeparatedJets) {
                if(itJet.getP4(aSystEffect).Pt()>30) ++myAnalyzer->nJets30;
        }

        myAnalyzer->nJetsInGap30 = 0;
        if(myAnalyzer->nJets30>=2) {
                for(unsigned int iJet=2; iJet<myAnalyzer->aSeparatedJets.size(); ++iJet) {
                        if( (myAnalyzer->aSeparatedJets.at(iJet).getP4().Eta()>myAnalyzer->aJet1.getP4().Eta() &&
                             myAnalyzer->aSeparatedJets.at(iJet).getP4().Eta()<myAnalyzer->aJet2.getP4().Eta()) ||
                            (myAnalyzer->aSeparatedJets.at(iJet).getP4().Eta()<myAnalyzer->aJet1.getP4().Eta() &&
                             myAnalyzer->aSeparatedJets.at(iJet).getP4().Eta()>myAnalyzer->aJet2.getP4().Eta()) ) {
                                if(myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Pt()>30) myAnalyzer->nJetsInGap30++;
                        }
                }
        }

        float jetsMass = (myAnalyzer->aJet1.getP4(aSystEffect)+myAnalyzer->aJet2.getP4(aSystEffect)).M();
        float higgsPt =  (myAnalyzer->aLeg1.getP4(aSystEffect) + myAnalyzer->aLeg2.getP4(aSystEffect) + myAnalyzer->aMET.getP4(aSystEffect)).Pt();
        bool mtSelection = myAnalyzer->aPair.getMTMuon(aSystEffect)<50;

        bool jet0_low =  myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>20 && myAnalyzer->aLeg2.getP4(aSystEffect).Pt()<50 && myAnalyzer->nJets30==0;
        bool jet0_high = myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>50 && myAnalyzer->nJets30==0;

        bool jet1_low = (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30==2 && jetsMass<500)) &&
                        (myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>30 && myAnalyzer->aLeg2.getP4(aSystEffect).Pt()<40 ||
                         myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>40 && higgsPt<140);

        bool jet1_high = (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30==2 && jetsMass<500)) &&
                         (myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>40 && higgsPt>140);

        bool vbf_low = myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>20 &&
                       myAnalyzer->nJets30==2 && jetsMass>500 &&
                       (jetsMass<800 || higgsPt<100);

        bool vbf_high = myAnalyzer->aLeg2.getP4(aSystEffect).Pt()>20 &&
                        myAnalyzer->nJets30==2 && jetsMass>800 && higgsPt>100;

        bool cpMuonSelection = myAnalyzer->aLeg1.getPCARefitPV().Perp()>myAnalyzer->nPCAMin_;
        bool cpTauSelection =  myAnalyzer->aLeg2.getPCARefitPV().Mag()>myAnalyzer->nPCAMin_;
        bool cpPi = cpMuonSelection && cpTauSelection && myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero;
        bool cpRho = cpMuonSelection && myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero &&
                     HTTAnalysis::isOneProng(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode));

        //2D categories
        bool jet0 = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>30 && myAnalyzer->nJets30 == 0;
        bool boosted = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>30 && (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30==2 && jetsMass < 300) || myAnalyzer->nJets30 > 2);
        bool vbf = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>30 && myAnalyzer->nJets30==2 && jetsMass>300;

        bool wSelection = myAnalyzer->aPair.getMTMuon(aSystEffect)>80 && myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)<0.15;
        bool ttSelection =  myAnalyzer->aPair.getMTMuon(aSystEffect)>150;
        bool muonAntiIso = myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)>0.15 && myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)<0.30;
        bool muonIso = myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)<0.15;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet0_low] = muonIso && mtSelection && jet0_low;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet0_high] = muonIso && mtSelection && jet0_high;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_low] = muonIso && mtSelection && jet1_low;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_high] = muonIso && mtSelection && jet1_high;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_low] = muonIso && mtSelection && vbf_low;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_high] = muonIso && mtSelection && vbf_high;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet0] = muonIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::CP_Pi] = muonIso && mtSelection && cpPi;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::CP_Rho] = muonIso && mtSelection && cpRho;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::boosted] = muonIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf] = muonIso && mtSelection && vbf;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::wjets_jet0] = muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::wjets_boosted] = muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::wjets_vbf] = muonIso && wSelection && vbf;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::qcd_jet0] = muonAntiIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::qcd_boosted] = muonAntiIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::qcd_vbf] = muonAntiIso && mtSelection && vbf;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::W] = muonIso && wSelection;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::TT] = muonIso && ttSelection;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuTauSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg1.getP4(aSystEffect).Eta(),
                                   myAnalyzer->aLeg1.getP4(aSystEffect).Pt(), HTTAnalysis::hadronicTauDecayModes::tauDecayMuon);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuTauSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg2.getP4(aSystEffect).Eta(), myAnalyzer->aLeg2.getP4(aSystEffect).Pt(),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)));

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
