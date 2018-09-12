
#include "MuMuSpecifics.h"

#include "AnalysisEnums.h"
#include "svfitAnalyzer.h"
#include "EventProxyHTT.h"
#include "Tools.h"

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
MuMuSpecifics::MuMuSpecifics(svfitAnalyzer * aAnalyzer) : ChannelSpecifics(aAnalyzer){

        decayModeName = "MuMu";
        tauID_FRSF_mu = new TH1F("tauID_FRSF_mu", "", 5, bins_mu_);
        Float_t binContents_mu[5] = {1.263, 1.364, 0.854, 1.712, 2.324};
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
void MuMuSpecifics::setAnalysisObjects(const EventProxyHTT & myEventProxy) {

        myAnalyzer->leg1_ = myAnalyzer->pair_.getMuon();
        myAnalyzer->leg2_ = myAnalyzer->pair_.getTau();

        myAnalyzer->genLeg1_ = HTTParticle();
        myAnalyzer->genLeg2_ = HTTParticle();

        if(myEventProxy.genLeptons &&
           myEventProxy.genLeptons->size()) {
                HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
                if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->genLeg1_ = aGenTau;
                if(myEventProxy.genLeptons->size()>1) {
                        aGenTau =  myEventProxy.genLeptons->at(1);
                        if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->genLeg2_ = aGenTau;
                }
        }

	if(myAnalyzer->leg1_.getP4().Perp() <  myAnalyzer->leg2_.getP4().Perp()) std::swap(myAnalyzer->leg1_, myAnalyzer->leg2_);
	if(myAnalyzer->genLeg1_.getP4().Perp() < myAnalyzer->genLeg2_.getP4().Perp()) std::swap(myAnalyzer->genLeg1_, myAnalyzer->genLeg2_);
	
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<bool, bool> MuMuSpecifics::checkTauDecayMode(const EventProxyHTT & myEventProxy){

        bool goodGenDecayMode = false;
        bool goodRecoDecayMode = false;

        std::vector<std::string> decayNamesGen = HTTAnalysis::getTauDecayName(myAnalyzer->genLeg2_.getProperty(PropertyEnum::decayMode),
                                                                              myAnalyzer->genLeg1_.getProperty(PropertyEnum::decayMode));
        std::vector<std::string> decayNamesReco = HTTAnalysis::getTauDecayName(myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode),HTTAnalysis::tauDecayMuon);

        for(auto it: decayNamesGen) if(it.find("Lepton1Prong")!=std::string::npos) goodGenDecayMode = true;
        for(auto it: decayNamesReco) if(it.find("Lepton1Prong")!=std::string::npos) goodRecoDecayMode = true;

        return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void MuMuSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it: myAnalyzer->categoryDecisions_) it = false;

        int tauIDmask = 0;
        for(unsigned int iBit=0; iBit<myAnalyzer->event_.ntauIds; iBit++) {
                if(myAnalyzer->event_.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
                if(myAnalyzer->event_.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }
        bool tauID = ( (int)myAnalyzer->leg2_.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;

        unsigned int muonIDmask = (1<<7);
        bool muonID = true;
        if(myAnalyzer->event_.getRunId()>278808 || myAnalyzer->event_.getRunId()==1) {
                muonID = ((int)myAnalyzer->leg1_.getProperty(PropertyEnum::muonID) & muonIDmask) == muonIDmask;
        }

        bool muonKinematics = myAnalyzer->leg1_.getP4(aSystEffect).Pt()>20 && fabs(myAnalyzer->leg1_.getP4(aSystEffect).Eta())<2.1;

        bool trigger = myAnalyzer->leg1_.hasTriggerMatch(TriggerEnum::HLT_IsoMu27) ||
                       myAnalyzer->leg1_.hasTriggerMatch(TriggerEnum::HLT_IsoMu24_eta2p1) ||
	               myAnalyzer->leg1_.hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);

        unsigned int metFilters = myAnalyzer->event_.getMETFilterDecision();
        unsigned int dataMask = (1<<8) -1;
	dataMask -= (1<<5);//Problem with globalTightHalo2016Filter
        unsigned int mcMask = dataMask - (1<<6) - (1<<7);
        bool metFilterDecision = (metFilters & mcMask) == mcMask;
        if(myAnalyzer->sampleName_=="Data") metFilterDecision = (metFilters & dataMask) == dataMask;
        metFilterDecision = true;//no met filters should be included - info from Laura Dodd

        if(!muonKinematics || !muonID || !tauID || !trigger || !metFilterDecision) return;

        myAnalyzer->nJets30_ = 0;
        for(auto itJet: myAnalyzer->separatedJets_) {
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

        float jetsMass = (myAnalyzer->jet1_.getP4(aSystEffect)+myAnalyzer->jet2_.getP4(aSystEffect)).M();
        float higgsPt =  (myAnalyzer->leg1_.getP4(aSystEffect) + myAnalyzer->leg2_.getP4(aSystEffect) + myAnalyzer->MET_.getP4(aSystEffect)).Pt();
        bool mtSelection = myAnalyzer->pair_.getMTMuon(aSystEffect)<50;
        bool cpMuonSelection = myAnalyzer->leg1_.getPCARefitPV().Mag()>0.05;
        bool cpTauSelection =  myAnalyzer->leg2_.getPCARefitPV().Mag()>0.05;
        bool cpPi = cpMuonSelection && cpTauSelection && myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero;
        bool cpRho = cpMuonSelection && myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero &&
                     HTTAnalysis::isOneProng(myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode));

        //2D categories
        bool jet0 = myAnalyzer->leg2_.getP4(aSystEffect).Perp()>30 && myAnalyzer->nJets30_ == 0;
        bool boosted = myAnalyzer->leg2_.getP4(aSystEffect).Perp()>30 && (myAnalyzer->nJets30_==1 || (myAnalyzer->nJets30_>=2 && (jetsMass < 300 || higgsPt < 50 || myAnalyzer->leg2_.getP4(aSystEffect).Perp()<40)));
        bool vbf = myAnalyzer->leg2_.getP4(aSystEffect).Perp()>40 && myAnalyzer->nJets30_>=2 && jetsMass>300 && higgsPt > 50;

        bool wSelection = myAnalyzer->pair_.getMTMuon(aSystEffect)>80;
        bool muonAntiIso = myAnalyzer->leg1_.getProperty(PropertyEnum::combreliso)>0.15 && myAnalyzer->leg1_.getProperty(PropertyEnum::combreliso)<0.30;
        bool muonIso = myAnalyzer->leg1_.getProperty(PropertyEnum::combreliso)<0.15;

        bool ss = myAnalyzer->leg2_.getCharge()*myAnalyzer->leg1_.getCharge() == 1;
        bool os = myAnalyzer->leg2_.getCharge()*myAnalyzer->leg1_.getCharge() == -1;
        //b-tag categories
        bool btag = false;
        bool nobtag = false;

        //Main categories
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->id()] = os && muonIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->id()] = os && muonIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->id()] = os && muonIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_jet0->id()] = os && muonAntiIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_boosted->id()] = os && muonAntiIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_vbf->id()] = os && muonAntiIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_pi->id()] = os && muonIso && mtSelection && cpPi;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_rho->id()] = os && muonIso && mtSelection && cpRho;

        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->id()] = os && muonIso && mtSelection;
        //myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_inclusive->id()] = os && muonAntiIso && mtSelection;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::btag->id()] = os && muonIso && mtSelection && btag;
        //myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_btag->id()] = os && muonAntiIso && mtSelection && btag;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::nobtag->id()] = os && muonIso && mtSelection && nobtag;
        //myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_nobtag->id()] = os && muonAntiIso && mtSelection && nobtag;


        //W control region
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->wSF()->id()] = os && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->wSF()->id()] = os && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->wSF()->id()] = os && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_jet0->wSF()->id()] = os && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_boosted->wSF()->id()] = os && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_vbf->wSF()->id()] = os && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_pi->wSF()->id()] = os && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_rho->wSF()->id()] = os && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->wSF()->id()] = os && muonIso && wSelection;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::btag->wSF()->id()] = os && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::nobtag->wSF()->id()] = os && muonIso && wSelection && nobtag;

        ///QCD region in W SF control region
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_jet0->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_boosted->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_vbf->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_pi->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_rho->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::btag->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::nobtag->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && nobtag;

        ///QCD region
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->qcdEstimate()->id()] = ss && muonIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->qcdEstimate()->id()] = ss && muonIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->qcdEstimate()->id()] = ss && muonIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_jet0->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_boosted->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_vbf->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_pi->qcdEstimate()->id()] = ss && muonIso && mtSelection && cpPi;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_rho->qcdEstimate()->id()] = ss && muonIso && mtSelection && cpRho;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->qcdEstimate()->id()] = ss && muonIso && mtSelection;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::btag->qcdEstimate()->id()] = ss && muonIso && mtSelection && btag;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::nobtag->qcdEstimate()->id()] = ss && muonIso && mtSelection && nobtag;

        ///W SF region in QCD region
        myAnalyzer->categoryDecisions_[ChannelSpecifics::jet0->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::boosted->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::vbf->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_jet0->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_boosted->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::antiIso_vbf->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_pi->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::mu_rho->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::inclusive->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::btag->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions_[ChannelSpecifics::nobtag->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && nobtag;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuMuSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->leg1_.getP4(aSystEffect).Eta(),
                                   myAnalyzer->leg1_.getP4(aSystEffect).Pt(),
                                   myAnalyzer->leg1_.getProperty(PropertyEnum::combreliso),
                                   HTTAnalysis::hadronicTauDecayModes::tauDecayMuon, false, myAnalyzer->leg1_.getProperty(PropertyEnum::mc_match));

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuMuSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        bool useTauTrigger = false, useXTrigger = false;
        if(myAnalyzer->leg1_.getP4(aSystEffect).Pt()<23) {useTauTrigger = true; useXTrigger = true;}

        return getLeptonCorrection(myAnalyzer->leg2_.getP4(aSystEffect).Eta(),
                                   myAnalyzer->leg2_.getP4(aSystEffect).Pt(),
                                   myAnalyzer->leg2_.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->leg2_.getProperty(PropertyEnum::decayMode)),
                                   useTauTrigger, myAnalyzer->leg2_.getProperty(PropertyEnum::mc_match), useXTrigger);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
