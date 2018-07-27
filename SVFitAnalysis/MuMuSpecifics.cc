
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

        myAnalyzer->aLeg1 = myAnalyzer->aPair.getMuon();
        myAnalyzer->aLeg2 = myAnalyzer->aPair.getTau();

        myAnalyzer->aGenLeg1 = HTTParticle();
        myAnalyzer->aGenLeg2 = HTTParticle();

        if(myEventProxy.genLeptons &&
           myEventProxy.genLeptons->size()) {
                HTTParticle aGenTau =  myEventProxy.genLeptons->at(0);
                if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->aGenLeg1 = aGenTau;
                if(myEventProxy.genLeptons->size()>1) {
                        aGenTau =  myEventProxy.genLeptons->at(1);
                        if(aGenTau.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecayMuon) myAnalyzer->aGenLeg2 = aGenTau;
                }
        }

	if(myAnalyzer->aLeg1.getP4().Perp() <  myAnalyzer->aLeg2.getP4().Perp()) std::swap(myAnalyzer->aLeg1, myAnalyzer->aLeg2);
	if(myAnalyzer->aGenLeg1.getP4().Perp() < myAnalyzer->aGenLeg2.getP4().Perp()) std::swap(myAnalyzer->aGenLeg1, myAnalyzer->aGenLeg2);
	
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
std::pair<bool, bool> MuMuSpecifics::checkTauDecayMode(const EventProxyHTT & myEventProxy){

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
void MuMuSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it: myAnalyzer->categoryDecisions) it = false;

        int tauIDmask = 0;
        for(unsigned int iBit=0; iBit<myAnalyzer->aEvent.ntauIds; iBit++) {
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
                if(myAnalyzer->aEvent.tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
        }
        bool tauID = ( (int)myAnalyzer->aLeg2.getProperty(PropertyEnum::tauID) & tauIDmask) == tauIDmask;

        unsigned int muonIDmask = (1<<7);
        bool muonID = true;
        if(myAnalyzer->aEvent.getRunId()>278808 || myAnalyzer->aEvent.getRunId()==1) {
                muonID = ((int)myAnalyzer->aLeg1.getProperty(PropertyEnum::muonID) & muonIDmask) == muonIDmask;
        }

        bool muonKinematics = myAnalyzer->aLeg1.getP4(aSystEffect).Pt()>20 && fabs(myAnalyzer->aLeg1.getP4(aSystEffect).Eta())<2.1;

        bool trigger = myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu27) ||
                       myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu24_eta2p1) ||
	               myAnalyzer->aLeg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1);

        unsigned int metFilters = myAnalyzer->aEvent.getMETFilterDecision();
        unsigned int dataMask = (1<<8) -1;
	dataMask -= (1<<5);//Problem with globalTightHalo2016Filter
        unsigned int mcMask = dataMask - (1<<6) - (1<<7);
        bool metFilterDecision = (metFilters & mcMask) == mcMask;
        if(myAnalyzer->sampleName=="Data") metFilterDecision = (metFilters & dataMask) == dataMask;
        metFilterDecision = true;//no met filters should be included - info from Laura Dodd

        if(!muonKinematics || !muonID || !tauID || !trigger || !metFilterDecision) return;

        myAnalyzer->nJets30 = 0;
        for(auto itJet: myAnalyzer->aSeparatedJets) {
                if(itJet.getP4(aSystEffect).Pt()>30) ++myAnalyzer->nJets30;
        }

        myAnalyzer->nJetsInGap30 = 0;
        if(myAnalyzer->nJets30>=2) {
                for(unsigned int iJet=2; iJet<myAnalyzer->aSeparatedJets.size(); ++iJet) {
                        if( (myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Eta()>myAnalyzer->aJet1.getP4(aSystEffect).Eta() &&
                             myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Eta()<myAnalyzer->aJet2.getP4(aSystEffect).Eta()) ||
                            (myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Eta()<myAnalyzer->aJet1.getP4(aSystEffect).Eta() &&
                             myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Eta()>myAnalyzer->aJet2.getP4(aSystEffect).Eta()) ) {
                                if(myAnalyzer->aSeparatedJets.at(iJet).getP4(aSystEffect).Pt()>30) myAnalyzer->nJetsInGap30++;
                        }
                }
        }

        float jetsMass = (myAnalyzer->aJet1.getP4(aSystEffect)+myAnalyzer->aJet2.getP4(aSystEffect)).M();
        float higgsPt =  (myAnalyzer->aLeg1.getP4(aSystEffect) + myAnalyzer->aLeg2.getP4(aSystEffect) + myAnalyzer->aMET.getP4(aSystEffect)).Pt();
        bool mtSelection = myAnalyzer->aPair.getMTMuon(aSystEffect)<50;
        bool cpMuonSelection = myAnalyzer->aLeg1.getPCARefitPV().Mag()>0.05;
        bool cpTauSelection =  myAnalyzer->aLeg2.getPCARefitPV().Mag()>0.05;
        bool cpPi = cpMuonSelection && cpTauSelection && myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)==HTTAnalysis::tauDecay1ChargedPion0PiZero;
        bool cpRho = cpMuonSelection && myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)!=HTTAnalysis::tauDecay1ChargedPion0PiZero &&
                     HTTAnalysis::isOneProng(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode));

        //2D categories
        bool jet0 = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>30 && myAnalyzer->nJets30 == 0;
        bool boosted = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>30 && (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30>=2 && (jetsMass < 300 || higgsPt < 50 || myAnalyzer->aLeg2.getP4(aSystEffect).Perp()<40)));
        bool vbf = myAnalyzer->aLeg2.getP4(aSystEffect).Perp()>40 && myAnalyzer->nJets30>=2 && jetsMass>300 && higgsPt > 50;

        bool wSelection = myAnalyzer->aPair.getMTMuon(aSystEffect)>80;
        bool muonAntiIso = myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)>0.15 && myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)<0.30;
        bool muonIso = myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso)<0.15;

        bool ss = myAnalyzer->aLeg2.getCharge()*myAnalyzer->aLeg1.getCharge() == 1;
        bool os = myAnalyzer->aLeg2.getCharge()*myAnalyzer->aLeg1.getCharge() == -1;
        //b-tag categories
        bool btag = false;
        bool nobtag = false;

        //Main categories
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->id()] = os && muonIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->id()] = os && muonIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->id()] = os && muonIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_jet0->id()] = os && muonAntiIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_boosted->id()] = os && muonAntiIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_vbf->id()] = os && muonAntiIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_pi->id()] = os && muonIso && mtSelection && cpPi;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_rho->id()] = os && muonIso && mtSelection && cpRho;

        myAnalyzer->categoryDecisions[ChannelSpecifics::inclusive->id()] = os && muonIso && mtSelection;
        //myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_inclusive->id()] = os && muonAntiIso && mtSelection;
        myAnalyzer->categoryDecisions[ChannelSpecifics::btag->id()] = os && muonIso && mtSelection && btag;
        //myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_btag->id()] = os && muonAntiIso && mtSelection && btag;
        myAnalyzer->categoryDecisions[ChannelSpecifics::nobtag->id()] = os && muonIso && mtSelection && nobtag;
        //myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_nobtag->id()] = os && muonAntiIso && mtSelection && nobtag;


        //W control region
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->wSF()->id()] = os && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->wSF()->id()] = os && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->wSF()->id()] = os && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_jet0->wSF()->id()] = os && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_boosted->wSF()->id()] = os && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_vbf->wSF()->id()] = os && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_pi->wSF()->id()] = os && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_rho->wSF()->id()] = os && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions[ChannelSpecifics::inclusive->wSF()->id()] = os && muonIso && wSelection;
        myAnalyzer->categoryDecisions[ChannelSpecifics::btag->wSF()->id()] = os && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions[ChannelSpecifics::nobtag->wSF()->id()] = os && muonIso && wSelection && nobtag;

        ///QCD region in W SF control region
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_jet0->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_boosted->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_vbf->wSF()->qcdEstimate()->id()] = ss && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_pi->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_rho->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions[ChannelSpecifics::inclusive->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection;
        myAnalyzer->categoryDecisions[ChannelSpecifics::btag->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions[ChannelSpecifics::nobtag->wSF()->qcdEstimate()->id()] = ss && muonIso && wSelection && nobtag;

        ///QCD region
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->qcdEstimate()->id()] = ss && muonIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->qcdEstimate()->id()] = ss && muonIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->qcdEstimate()->id()] = ss && muonIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_jet0->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_boosted->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_vbf->qcdEstimate()->id()] = ss && muonAntiIso && mtSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_pi->qcdEstimate()->id()] = ss && muonIso && mtSelection && cpPi;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_rho->qcdEstimate()->id()] = ss && muonIso && mtSelection && cpRho;
        myAnalyzer->categoryDecisions[ChannelSpecifics::inclusive->qcdEstimate()->id()] = ss && muonIso && mtSelection;
        myAnalyzer->categoryDecisions[ChannelSpecifics::btag->qcdEstimate()->id()] = ss && muonIso && mtSelection && btag;
        myAnalyzer->categoryDecisions[ChannelSpecifics::nobtag->qcdEstimate()->id()] = ss && muonIso && mtSelection && nobtag;

        ///W SF region in QCD region
        myAnalyzer->categoryDecisions[ChannelSpecifics::jet0->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::boosted->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::vbf->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_jet0->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && jet0;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_boosted->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && boosted;
        myAnalyzer->categoryDecisions[ChannelSpecifics::antiIso_vbf->qcdEstimate()->wSF()->id()] = ss && muonAntiIso && wSelection && vbf;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_pi->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && cpPi;
        myAnalyzer->categoryDecisions[ChannelSpecifics::mu_rho->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && cpRho;
        myAnalyzer->categoryDecisions[ChannelSpecifics::inclusive->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection;
        myAnalyzer->categoryDecisions[ChannelSpecifics::btag->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && btag;
        myAnalyzer->categoryDecisions[ChannelSpecifics::nobtag->qcdEstimate()->wSF()->id()] = ss && muonIso && wSelection && nobtag;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuMuSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg1.getP4(aSystEffect).Eta(),
                                   myAnalyzer->aLeg1.getP4(aSystEffect).Pt(),
                                   myAnalyzer->aLeg1.getProperty(PropertyEnum::combreliso),
                                   HTTAnalysis::hadronicTauDecayModes::tauDecayMuon, false, myAnalyzer->aLeg1.getProperty(PropertyEnum::mc_match));

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float MuMuSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        bool useTauTrigger = false, useXTrigger = false;
        if(myAnalyzer->aLeg1.getP4(aSystEffect).Pt()<23) {useTauTrigger = true; useXTrigger = true;}

        return getLeptonCorrection(myAnalyzer->aLeg2.getP4(aSystEffect).Eta(),
                                   myAnalyzer->aLeg2.getP4(aSystEffect).Pt(),
                                   myAnalyzer->aLeg2.getProperty(PropertyEnum::byIsolationMVArun2v1DBoldDMwLTraw),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)),
                                   useTauTrigger, myAnalyzer->aLeg2.getProperty(PropertyEnum::mc_match), useXTrigger);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
