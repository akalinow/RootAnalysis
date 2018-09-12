#include "TauTauSpecifics.h"
#include "svfitAnalyzer.h"
#include "EventProxyHTT.h"
#include "Tools.h"


TauTauSpecifics::TauTauSpecifics(svfitAnalyzer * aAnalyzer) : ChannelSpecifics(aAnalyzer){

        decayModeName = "TauTau";
        tauID_FRSF_mu = new TH1F("tauID_FRSF_mu", "", 5, bins_mu_);
        Float_t binContents_mu[5] = {1.010, 1.007, 0.870, 1.154, 2.281};
        for(int i=0; i<tauID_FRSF_mu->GetNbinsX(); i++){
          tauID_FRSF_mu->SetBinContent(i+1, binContents_mu[i]);
          }
        tauID_FRSF_ele = new TH1F("tauID_FRSF_ele", "", 3, bins_ele_);
        Float_t binContents_ele[3] = {1.213, 1, 1.375};
        for(int i=0; i<tauID_FRSF_ele->GetNbinsX(); i++){
          tauID_FRSF_ele->SetBinContent(i+1, binContents_ele[i]);
          }

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void TauTauSpecifics::setAnalysisObjects(const EventProxyHTT & myEventProxy) {

        myAnalyzer->leg1_ = myAnalyzer->pair_.getLeg1();
        myAnalyzer->leg2_ = myAnalyzer->pair_.getLeg2();

        myAnalyzer->genLeg1_ = HTTParticle();
        myAnalyzer->genLeg2_ = HTTParticle();

        if(myEventProxy.genLeptons && myEventProxy.genLeptons->size()) {
                HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
                if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) {
                        if(aGenTau.getP4().DeltaR(myAnalyzer->leg1_.getP4())<0.5)
                                myAnalyzer->genLeg1_ = aGenTau;
                        else if(aGenTau.getP4().DeltaR(myAnalyzer->leg2_.getP4())<0.5)
                                myAnalyzer->genLeg2_ = aGenTau;
                }
                if(myEventProxy.genLeptons->size()>1) {
                        aGenTau =  myEventProxy.genLeptons->at(1);
                        if(aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecayMuon && aGenTau.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecaysElectron) {
                                if(aGenTau.getP4().DeltaR(myAnalyzer->leg1_.getP4())<0.5)
                                        myAnalyzer->genLeg1_ = aGenTau;
                                else if(aGenTau.getP4().DeltaR(myAnalyzer->leg2_.getP4())<0.5)
                                        myAnalyzer->genLeg2_ = aGenTau;
                        }
                }
        }
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<bool, bool> TauTauSpecifics::checkTauDecayMode(const EventProxyHTT & myEventProxy){

        bool goodGenDecayMode = false;
        bool goodRecoDecayMode = false;

        std::vector<std::string> decayNamesGen = HTTAnalysis::getTauDecayName(myAnalyzer->genLeg1_.getProperty(PropertyEnum::decayMode),
                                                                              myAnalyzer->genLeg2_.getProperty(PropertyEnum::decayMode));
        std::vector<std::string> decayNamesReco = HTTAnalysis::getTauDecayName(myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode),HTTAnalysis::tauDecayMuon);

        for(auto it : decayNamesGen) if(it.find("1Prong1Prong")!=std::string::npos) goodGenDecayMode = true;
        for(auto it : decayNamesReco) if(it.find("1Prong1Prong")!=std::string::npos) goodRecoDecayMode = true;

        return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void TauTauSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it : myAnalyzer->categoryDecisions_) it = false;

        ///This stands for core selection, that is common to all regions.
        bool tau1Kinematics = myAnalyzer->leg1_.getP4(aSystEffect).Pt()>50 && std::abs(myAnalyzer->leg1_.getP4(aSystEffect).Eta())<2.1;
        bool tau2Kinematics = myAnalyzer->leg2_.getP4(aSystEffect).Pt()>40 && std::abs(myAnalyzer->leg2_.getP4(aSystEffect).Eta())<2.1;

        int tauIDmask=0, tauIsoTmask=0, tauIsoMmask=0, tauIsoLmask=0;

        for(unsigned int iBit=0; iBit<myAnalyzer->event_.ntauIds; iBit++) {
                if(myAnalyzer->event_.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIsoTmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="byMediumIsolationMVArun2v1DBoldDMwLT") tauIsoMmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="byLooseIsolationMVArun2v1DBoldDMwLT") tauIsoLmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="againstMuonLoose3") tauIDmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }

        bool tau1ID = ( (int)myAnalyzer->leg1_.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool tau1IsoT = ( (int)myAnalyzer->leg1_.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
        bool tau1IsoM = ( (int)myAnalyzer->leg1_.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
        bool tau1IsoL = ( (int)myAnalyzer->leg1_.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;
        bool tau2ID = ( (int)myAnalyzer->leg2_.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;
        bool tau2IsoT = ( (int)myAnalyzer->leg2_.getProperty(PropertyEnum::tauID) & tauIsoTmask) == tauIsoTmask;
        bool tau2IsoM = ( (int)myAnalyzer->leg2_.getProperty(PropertyEnum::tauID) & tauIsoMmask) == tauIsoMmask;
        bool tau2IsoL = ( (int)myAnalyzer->leg2_.getProperty(PropertyEnum::tauID) & tauIsoLmask) == tauIsoLmask;

        bool fullIso = tau1IsoT && tau2IsoT;
        bool relaxedIso = (tau1IsoM && tau2IsoL) || (tau2IsoM && tau1IsoL);
        bool antiIso = (tau1IsoM && tau2IsoL && !tau2IsoT) || (tau2IsoM && tau1IsoL && !tau1IsoT);

        bool trigger1 = myAnalyzer->leg1_.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg) &&
	               myAnalyzer->leg2_.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg);

        bool trigger2 = myAnalyzer->leg1_.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg) &&
                              myAnalyzer->leg2_.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg);

        bool trigger = trigger1 || trigger2;

        unsigned int metFilters = myAnalyzer->event_.getMETFilterDecision();
        unsigned int dataMask = (1<<8) -1;
	dataMask -= (1<<5);//Problem with globalTightHalo2016Filter
        unsigned int mcMask = dataMask - (1<<6) - (1<<7);
        bool metFilterDecision = (metFilters & mcMask) == mcMask;
        if(myAnalyzer->sampleName_=="Data") metFilterDecision = (metFilters & dataMask) == dataMask;
        if(!tau1Kinematics || !tau1ID || !tau2Kinematics || !tau2ID || !relaxedIso || !trigger || !metFilterDecision) return;

        myAnalyzer->nJets30_ = 0;
        for(auto itJet : myAnalyzer->separatedJets_) {
                if(itJet.getP4(aSystEffect).Pt()>30) ++myAnalyzer->nJets30_;
        }

        myAnalyzer->nJetsInGap30_ = 0;
        if(myAnalyzer->nJets30_>=2) {
                for(unsigned int iJet=2; iJet<myAnalyzer->separatedJets_.size(); ++iJet) {
                        if( (myAnalyzer->separatedJets_.at(iJet).getP4(aSystEffect).Eta()>myAnalyzer->jet1_.getP4(aSystEffect).Eta() &&
                             myAnalyzer->separatedJets_.at(iJet).getP4(aSystEffect).Eta()<myAnalyzer->jet2_.getP4(aSystEffect).Eta()) ||
                            (myAnalyzer->separatedJets_.at(iJet).getP4(aSystEffect).Eta()<myAnalyzer->jet1_.getP4(aSystEffect).Eta() &&
                             myAnalyzer->separatedJets_.at(iJet).getP4(aSystEffect).Eta()>myAnalyzer->jet2_.getP4(aSystEffect).Eta()) ) {
                                if(myAnalyzer->separatedJets_.at(iJet).getP4(aSystEffect).Pt()>30) myAnalyzer->nJetsInGap30_++;
                        }
                }
        }

        float jetsEta = std::abs(myAnalyzer->jet1_.getP4().Eta() - myAnalyzer->jet2_.getP4().Eta());
        float higgsPt =  (myAnalyzer->leg1_.getP4(aSystEffect) + myAnalyzer->leg2_.getP4(aSystEffect) + myAnalyzer->MET_.getP4(aSystEffect)).Pt();

        bool ss = myAnalyzer->leg2_.getCharge()*myAnalyzer->leg1_.getCharge() == 1;
        bool os = myAnalyzer->leg2_.getCharge()*myAnalyzer->leg1_.getCharge() == -1;

        bool jet0 = myAnalyzer->nJets30_==0;
        bool boosted = myAnalyzer->nJets30_==1 || (myAnalyzer->nJets30_>=2 && !(jetsEta>2.5 && higgsPt>100));
        bool vbf_2d = myAnalyzer->nJets30_>=2 && (jetsEta>2.5 && higgsPt>100);

        //////////

        //Main categories
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->id()] = os && fullIso && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->id()] = os && fullIso && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->id()] = os && fullIso && vbf_2d;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->id()] = os && fullIso;


        ///QCD region
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->qcdEstimate()->id()] = os && antiIso && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->qcdEstimate()->id()] = os && antiIso &&  boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->qcdEstimate()->id()] = os && antiIso && vbf_2d;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->qcdEstimate()->id()] = os && antiIso;


        ///QCD SF denominator
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->qcdSFDenominator()->id()] = ss && antiIso && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->qcdSFDenominator()->id()] = ss && antiIso &&  boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->qcdSFDenominator()->id()] = ss && antiIso && vbf_2d;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->qcdSFDenominator()->id()] = ss && antiIso;


        ///QCD SF numerator
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->qcdSFNumerator()->id()] = ss && fullIso && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->qcdSFNumerator()->id()] = ss && fullIso &&  boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->qcdSFNumerator()->id()] = ss && fullIso && vbf_2d;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->qcdSFNumerator()->id()] = ss && fullIso;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->leg1_.getP4().Eta(),
                                   myAnalyzer->leg1_.getP4().Pt(),
                                   myAnalyzer->leg1_.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->leg1_.getProperty(PropertyEnum::decayMode)),true,
                                   myAnalyzer->leg1_.getProperty(PropertyEnum::mc_match));
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->leg2_.getP4().Eta(),
                                   myAnalyzer->leg2_.getP4().Pt(),
                                   myAnalyzer->leg2_.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode)),true,
                                   myAnalyzer->leg2_.getProperty(PropertyEnum::mc_match));
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
