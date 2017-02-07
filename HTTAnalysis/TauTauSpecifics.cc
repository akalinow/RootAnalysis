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

        for(auto it: decayNamesGen) if(it.find("1Prong1Prong")!=std::string::npos) goodGenDecayMode = true;
        for(auto it: decayNamesReco) if(it.find("1Prong1Prong")!=std::string::npos) goodRecoDecayMode = true;

        return std::pair<bool, bool>(goodGenDecayMode, goodRecoDecayMode);
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void TauTauSpecifics::testAllCategories(const HTTAnalysis::sysEffects & aSystEffect){

        for(auto && it: myAnalyzer->categoryDecisions) it = false;

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
        float jetsEta = std::abs(myAnalyzer->aJet1.getP4().Eta() - myAnalyzer->aJet2.getP4().Eta());
        float higgsPt =  (myAnalyzer->aLeg1.getP4(aSystEffect) + myAnalyzer->aLeg2.getP4(aSystEffect) + myAnalyzer->aMET.getP4(aSystEffect)).Pt();


  bool jet0 = myAnalyzer->nJets30==0;

  bool jet1 = (myAnalyzer->nJets30==1 || (myAnalyzer->nJets30>=2 && !(jetsMass>300 && jetsEta>2.5 && myAnalyzer->nJetsInGap30<1)));
  bool jet1_low = jet1 && (higgsPt>100 && higgsPt<170);

  bool jet1_high = jet1 && higgsPt>170;

  bool vbf_1d = myAnalyzer->nJets30>=2 && jetsEta>2.5 && myAnalyzer->nJetsInGap30<1;
  bool vbf_low =  vbf_1d &&
    ((higgsPt<100 && jetsMass>300) || (higgsPt>100 && jetsMass>300 && jetsMass<500));

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

  myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_low] = jet1_low;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet1_high] = jet1_high;

  myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_low] = vbf_low;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf_high] = vbf_high;

  myAnalyzer->categoryDecisions[(int)HTTAnalysis::jet0] = jet0;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::boosted] = boosted;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::vbf] = vbf_2d;

  myAnalyzer->categoryDecisions[(int)HTTAnalysis::pipi] = piPi;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::pirho] = piRho;
  myAnalyzer->categoryDecisions[(int)HTTAnalysis::rhorho] = rhoRho;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg1Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg1.getP4(aSystEffect).Eta(),
                                   myAnalyzer->aLeg1.getP4(aSystEffect).Pt(), HTTAnalysis::hadronicTauDecayModes::tauDecayMuon);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float TauTauSpecifics::getLeg2Correction(const HTTAnalysis::sysEffects & aSystEffect){

        return getLeptonCorrection(myAnalyzer->aLeg2.getP4(aSystEffect).Eta(), myAnalyzer->aLeg2.getP4(aSystEffect).Pt(),
                                   static_cast<HTTAnalysis::hadronicTauDecayModes>(myAnalyzer->aLeg2.getProperty(PropertyEnum::decayMode)));

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
