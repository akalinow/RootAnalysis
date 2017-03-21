#include "TauTauSpecifics.h"
#include "HTTAnalyzer.h"
#include "EventProxyHTT.h"
#include "Tools.h"


TauTauSpecifics::TauTauSpecifics(HTTAnalyzer * aAnalyzer) : ChannelSpecifics(aAnalyzer){

        decayModeName = "TauTau";

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void TauTauSpecifics::setAnalysisObjects(const EventProxyHTT & myEventProxy) {

        myAnalyzer->aLeg1 = myAnalyzer->aPair.getLeg1();
        myAnalyzer->aLeg2 = myAnalyzer->aPair.getLeg2();

        myAnalyzer->aGenLeg1 = HTTParticle();
        myAnalyzer->aGenLeg2 = HTTParticle();

        if(myEventProxy.genLeptons && myEventProxy.genLeptons->size()) {
                HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
                if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) {
                        if(aGenTau.getP4().DeltaR(myAnalyzer->aLeg1.getP4())<0.5)
                                myAnalyzer->aGenLeg1 = aGenTau;
                        else if(aGenTau.getP4().DeltaR(myAnalyzer->aLeg2.getP4())<0.5)
                                myAnalyzer->aGenLeg2 = aGenTau;
                }
                if(myEventProxy.genLeptons->size()>1) {
                        aGenTau =  myEventProxy.genLeptons->at(1);
                        if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) {
                                if(aGenTau.getP4().DeltaR(myAnalyzer->aLeg1.getP4())<0.5)
                                        myAnalyzer->aGenLeg1 = aGenTau;
                                else if(aGenTau.getP4().DeltaR(myAnalyzer->aLeg2.getP4())<0.5)
                                        myAnalyzer->aGenLeg2 = aGenTau;
                        }
                }
        }
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<bool, bool> TauTauSpecifics::checkTauDecayMode(const EventProxyHTT & myEventProxy){

        bool goodGenDecayMode = false;
        bool goodRecoDecayMode = false;

        std::vector<std::string> decayNamesGen = HTTAnalysis::getTauDecayName(myAnalyzer->aGenLeg1.getProperty(PropertyEnum::decayMode),
                                                                              myAnalyzer->aGenLeg2.getProperty(PropertyEnum::decayMode));
        std::vector<std::string> decayNamesReco = HTTAnalysis::getTauDecayName(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode),HTTAnalysis::tauDecayMuon);

        for(auto it : decayNamesGen) if(it.find("1Prong1Prong")!=std::string::npos) goodGenDecayMode = true;
        for(auto it : decayNamesReco) if(it.find("1Prong1Prong")!=std::string::npos) goodRecoDecayMode = true;

        return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void TauTauSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it : myAnalyzer->categoryDecisions) it = false;

        ///This stands for core selection, that is common to all regions.
        bool tau1Kinematics = myAnalyzer->aLeg1.getP4().Pt()>50 && std::abs(myAnalyzer->aLeg1.getP4().Eta())<2.1;
        bool tau2Kinematics = myAnalyzer->aLeg2.getP4().Pt()>40 && std::abs(myAnalyzer->aLeg2.getP4().Eta())<2.1;

        int tauIDmask=0, tauIsoTmask=0, tauIsoMmask=0, tauIsoLmask=0;

        for(unsigned int iBit=0; iBit<myAnalyzer->aEvent.ntauIds; iBit++) {
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIsoTmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="byMediumIsolationMVArun2v1DBoldDMwLT") tauIsoMmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="byLooseIsolationMVArun2v1DBoldDMwLT") tauIsoLmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstMuonLoose3") tauIDmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }

        bool tau1ID = ( (int)myAnalyzer->aLeg1.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool tau1IsoT = ( (int)myAnalyzer->aLeg1.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
        bool tau1IsoM = ( (int)myAnalyzer->aLeg1.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
        bool tau1IsoL = ( (int)myAnalyzer->aLeg1.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;
        bool tau2ID = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool tau2IsoT = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
        bool tau2IsoM = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
        bool tau2IsoL = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;

        bool fullIso = tau1IsoT && tau2IsoT;
        bool relaxedIso = (tau1IsoM && tau2IsoL) || (tau2IsoM && tau1IsoL);
        bool antiIso = (tau1IsoM && tau2IsoL && !tau2IsoT) || (tau2IsoM && tau1IsoL && !tau1IsoT);

        bool mediumIsoTrigger = myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg) &&
                                myAnalyzer->aLeg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);

        bool mediumCombinedIsoTrigger = myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg) &&
                                        myAnalyzer->aLeg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg);

        bool trigger = mediumIsoTrigger || mediumCombinedIsoTrigger;
        if(myAnalyzer->sampleName!="Data") trigger = true; //MC trigger included in SF

        if(!tau1Kinematics || !tau1ID || !tau2Kinematics || !tau2ID || !relaxedIso || !trigger) return;

        myAnalyzer->nJets30 = 0;
        for(auto itJet : myAnalyzer->aSeparatedJets) {
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
        float jetsEta = std::abs(myAnalyzer->aJet1.getP4().Eta() - myAnalyzer->aJet2.getP4().Eta());
        float higgsPt =  (myAnalyzer->aLeg1.getP4(aSystEffect) + myAnalyzer->aLeg2.getP4(aSystEffect) + myAnalyzer->aMET.getP4(aSystEffect)).Pt();

        bool ss = myAnalyzer->aLeg2.getCharge()*myAnalyzer->aLeg1.getCharge() == 1;
        bool os = myAnalyzer->aLeg2.getCharge()*myAnalyzer->aLeg1.getCharge() == -1;

        bool jet0 = myAnalyzer->nJets30==0;

        bool jet1 = (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30>=2 && !(jetsMass>300 && jetsEta>2.5 && myAnalyzer->nJetsInGap30<1)));
        bool jet1_low = jet1 && (higgsPt>100 && higgsPt<170);

        bool jet1_high = jet1 && higgsPt>170;

        bool vbf_1d = myAnalyzer->nJets30>=2 && jetsEta>2.5 && myAnalyzer->nJetsInGap30<1;
        bool vbf_low =  vbf_1d && ((higgsPt<100 && jetsMass>300) || (higgsPt>100 && jetsMass>300 && jetsMass<500));
        bool vbf_high = vbf_1d && (higgsPt>100 && jetsMass>500);
        bool boosted = myAnalyzer->nJets30==1 || (myAnalyzer->nJets30>=2 && !(jetsEta>2.5 && myAnalyzer->nJetsInGap30<1 && higgsPt>100) );
        bool vbf_2d = myAnalyzer->nJets30>=2 && (jetsEta>2.5 && myAnalyzer->nJetsInGap30<1 && higgsPt>100);

        //////////
        // categories by tau decay modes for CP
        bool isPi1 = (myAnalyzer->aLeg1.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero && myAnalyzer->aLeg1.getPCARefitPV().Mag()>myAnalyzer->nPCAMin_);
        bool isPi2 = (myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero && myAnalyzer->aLeg2.getPCARefitPV().Mag()>myAnalyzer->nPCAMin_);
        bool isRho1 = (myAnalyzer->aLeg1.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero && HTTAnalysis::isOneProng(myAnalyzer->aLeg1.getProperty(PropertyEnum::decayMode)) );
        bool isRho2 = (myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero && HTTAnalysis::isOneProng(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)) );

        bool piPi = isPi1 && isPi2;
        bool piRho = (isPi1 && isRho2) || (isPi2 && isRho1);
        bool rhoRho = isRho1 && isRho2;

/*
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_low] = os && jet1_low;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_high] = os && jet1_high;

        myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_low] = os && vbf_low;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_high] = os && vbf_high;
*/

         //Main categories
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->id()] = os && fullIso && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->id()] = os && fullIso && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->id()] = os && fullIso && vbf_2d;

        ///QCD region
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->qcdEstimate()->id()] = os && antiIso && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->qcdEstimate()->id()] = os && antiIso &&  boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->qcdEstimate()->id()] = os && antiIso && vbf_2d;

        ///QCD SF denominator
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->qcdSFDenominator()->id()] = ss && antiIso && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->qcdSFDenominator()->id()] = ss && antiIso &&  boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->qcdSFDenominator()->id()] = ss && antiIso && vbf_2d;

        ///QCD SF numerator
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->qcdSFNumerator()->id()] = ss && fullIso && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->qcdSFNumerator()->id()] = ss && fullIso &&  boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->qcdSFNumerator()->id()] = ss && fullIso && vbf_2d;
/*
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::pi_pi] = os && fullIso && piPi;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::pi_rho] = os && fullIso && piRho;
        myAnalyzer->categoryDecisions[(int)HTTAnalysis::rho_rho] = os && fullIso && rhoRho;
        */
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg1.getP4().Eta(),
                                   myAnalyzer->aLeg1.getP4().Pt(),
                                   myAnalyzer->aLeg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->aLeg1.getProperty(PropertyEnum::decayMode)),true);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg2.getP4().Eta(),
                                   myAnalyzer->aLeg2.getP4().Pt(),
                                   myAnalyzer->aLeg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)),true);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
