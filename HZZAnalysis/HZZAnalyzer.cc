#include <sstream>
#include <bitset>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"

#include "HZZAnalyzer.h"
#include "HZZHistograms.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HZZAnalyzer::HZZAnalyzer(const std::string & aName, const std::string & aDecayMode) : Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HZZAnalyzer::~HZZAnalyzer(){

        if(myHistos_) delete myHistos_;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HZZAnalyzer::clone() const {

        std::string myDecayMode = "4Mu";
        HZZAnalyzer* clone = new HZZAnalyzer(name(),myDecayMode);
        clone->setHistos(myHistos_);
        return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::initialize(TDirectory* aDir,
                             pat::strbitset *aSelections){

        mySelections_ = aSelections;

        myHistos_ = new HZZHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::finalize(){

        myHistos_->finalizeHistograms();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::addBranch(TTree *tree){

        std::string leafList = "mass4Mu/f:massZ1/f:massZ2/f";
        leafList+=":muon1ID/f:muon2ID/f:muon3ID/f:muon4ID/f";
        leafList+=":muon1Isolation/f:muon2Isolation/f:muon3Isolation/f:muon4Isolation/f";
        leafList+=":muon1Pt/f:muon2Pt/f:muon3Pt/f:muon4Pt/f";
        leafList+=":muon1SIP/f:muon2SIP/f:muon3SIP/f:muon4SIP/f";

        tree->Branch("entry",&aEntry,leafList.c_str());

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HZZAnalyzer::analyze(const EventProxyBase& iEvent){

        const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        bool applyCuts = true;

        bool singleMuTrigger = false;
        bool muonID = true;
        unsigned int muonIDmask = (1<<6);
        std::vector<const HTTParticle *> myMuonsPlus, myMuonsMinus, myMuons;
        for(unsigned int iLepton=0; iLepton<myEventProxy.leptons->size(); ++iLepton) {
                const HTTParticle * aLepton = &myEventProxy.leptons->at(iLepton);

                int pdgId = aLepton->getProperty(PropertyEnum::PDGId);
                if(std::abs(pdgId)!=13) continue;
                muonID = (int)aLepton->getProperty(PropertyEnum::muonID) & muonIDmask;

                singleMuTrigger |= aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1);



                if(aLepton->getP4().Perp()<5 || std::abs(aLepton->getP4().Eta())>2.4) continue;
                if(applyCuts) {
                        if(!muonID) continue;
                        if(aLepton->getProperty(PropertyEnum::combreliso)>0.35) continue;
                        if(std::abs(aLepton->getProperty(PropertyEnum::dxy))>0.5) continue;
                        if(std::abs(aLepton->getProperty(PropertyEnum::dz))>1) continue;
                        if(std::abs(aLepton->getProperty(PropertyEnum::SIP))>4) continue;
                }
                if(pdgId==-13) myMuonsPlus.push_back(aLepton);
                if(pdgId==13) myMuonsMinus.push_back(aLepton);
                myMuons.push_back(aLepton);
        }

        if(myMuonsPlus.size()!=2 || myMuonsMinus.size()!=2) return false;

        if(applyCuts) {
                if(fileName.find("Single")==std::string::npos && singleMuTrigger) return false;
                if(fileName.find("Single")!=std::string::npos && !singleMuTrigger) return false;
        }

        float deltaR;
        deltaR = myMuonsPlus[0]->getP4().DeltaR(myMuonsPlus[1]->getP4());
        if(deltaR<0.02) return false;

        deltaR = myMuonsMinus[0]->getP4().DeltaR(myMuonsMinus[1]->getP4());
        if(deltaR<0.02) return false;

        for(auto aMuonPlus: myMuonsPlus) {
                for(auto aMuonMinus: myMuonsMinus) {
                        deltaR = aMuonPlus->getP4().DeltaR(aMuonMinus->getP4());
                        if(deltaR<0.02) return false;
                }
        }

        unsigned int nPt10 = 0, nPt20 = 0;

        if(applyCuts) {
                for(auto it:myMuonsPlus) {
                        if(it->getP4().Perp()>10) ++nPt10;
                        if(it->getP4().Perp()>20) ++nPt20;
                }
                for(auto it:myMuonsMinus) {
                        if(it->getP4().Perp()>10) ++nPt10;
                        if(it->getP4().Perp()>20) ++nPt20;
                }
                if(nPt10<2 || nPt20<1) return false;
        }

        TLorentzVector p4Z1, p4Z2, p4ZZ;
        for(auto it:myMuonsPlus) p4ZZ+=it->getP4();
        for(auto it:myMuonsMinus) p4ZZ+=it->getP4();

        float deltaMass = 9999;
        float deltaTmp = 9999;

        TLorentzVector p4Tmp;
        for(auto aMuonPlus: myMuonsPlus) {
                for(auto aMuonMinus: myMuonsMinus) {
                        p4Tmp = aMuonPlus->getP4() + aMuonMinus->getP4();
                        if(p4Tmp.M()<4) return false;
                        deltaTmp = std::abs(p4Tmp.M() - 91.1896);
                        if(deltaTmp<deltaMass) {
                                deltaMass = deltaTmp;
                                p4Z1 = p4Tmp;
                                p4Z2 = p4ZZ - p4Z1;
                        }
                }
        }

        if(applyCuts) {
                if(p4Z1.M()>120 || p4Z1.M()<12) return false;
                if(p4Z2.M()>120 || p4Z2.M()<12) return false;
                if(p4Z1.M()<40) return false;
        }

        aEntry.mass4Mu = p4ZZ.M();
        aEntry.massZ1 = p4Z1.M();
        aEntry.massZ2 = p4Z2.M();

        aEntry.muon1Pt = myMuons[0]->getP4().Perp();
        aEntry.muon2Pt = myMuons[1]->getP4().Perp();
        aEntry.muon3Pt = myMuons[2]->getP4().Perp();
        aEntry.muon4Pt = myMuons[3]->getP4().Perp();

        aEntry.muon1Isolation = myMuons[0]->getProperty(PropertyEnum::combreliso);
        aEntry.muon2Isolation = myMuons[1]->getProperty(PropertyEnum::combreliso);
        aEntry.muon3Isolation = myMuons[2]->getProperty(PropertyEnum::combreliso);
        aEntry.muon4Isolation = myMuons[3]->getProperty(PropertyEnum::combreliso);

        aEntry.muon1ID = myMuons[0]->getProperty(PropertyEnum::muonID);
        aEntry.muon2ID = myMuons[1]->getProperty(PropertyEnum::muonID);
        aEntry.muon3ID = myMuons[2]->getProperty(PropertyEnum::muonID);
        aEntry.muon4ID = myMuons[3]->getProperty(PropertyEnum::muonID);

        aEntry.muon1SIP = myMuons[0]->getProperty(PropertyEnum::SIP);
        aEntry.muon2SIP = myMuons[1]->getProperty(PropertyEnum::SIP);
        aEntry.muon3SIP = myMuons[2]->getProperty(PropertyEnum::SIP);
        aEntry.muon4SIP = myMuons[3]->getProperty(PropertyEnum::SIP);

        myHistos_->fill1DHistogram("h1DMassZ1",p4Z1.M());
        myHistos_->fill1DHistogram("h1DMassZ2",p4Z2.M());
        myHistos_->fill1DHistogram("h1DMass4Mu",p4ZZ.M());

        myHistos_->fill1DHistogram("h1DPtMu1",myMuons[0]->getP4().Perp());
        myHistos_->fill1DHistogram("h1DPtMu2",myMuons[1]->getP4().Perp());
        myHistos_->fill1DHistogram("h1DPtMu3",myMuons[2]->getP4().Perp());
        myHistos_->fill1DHistogram("h1DPtMu4",myMuons[3]->getP4().Perp());

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
