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
        ntupleFile_ = 0;
        hStatsFromFile_ = 0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::finalize(){

        myHistos_->finalizeHistograms();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::addBranch(TTree *tree){

        std::string leafList = "mass4Mu/f";
        leafList+=":massZ1/f:massZ2/f";
        leafList+=":nMuons/f:";
        leafList+=":muon1Pt/f:muon2Pt/f:muon3Pt/f:muon4Pt/f";
        leafList+=":muon1Eta/f:muon2Eta/f:muon3Eta/f:muon4Eta/f";
        leafList+=":muon1Phi/f:muon2Phi/f:muon3Phi/f:muon4Phi/f";
        leafList+=":muon1Dxy/f:muon2Dxy/f:muon3Dxy/f:muon4Dxy/f";
        leafList+=":muon1Dz/f:muon2Dz/f:muon3Dz/f:muon4Dz/f";
        leafList+=":muon1SIP/f:muon2SIP/f:muon3SIP/f:muon4SIP/f";
        leafList+=":muon1Isol/f:muon2Isol/f:muon3Isol/f:muon4Isol/f";

        tree->Branch("entry",&aEntry,leafList.c_str());

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HZZAnalyzer::getPreselectionEff(const EventProxyHTT & myEventProxy){

    if(omp_get_thread_num()!=0) return;

          if(ntupleFile_!=myEventProxy.getTTree()->GetCurrentFile()) {
                ntupleFile_ = myEventProxy.getTTree()->GetCurrentFile();
                if(hStatsFromFile_) delete hStatsFromFile_;
                hStatsFromFile_ = (TH1F*)ntupleFile_->Get("hStats");

                std::string hName = "h1DStats";
                TH1F *hStats = myHistos_->get1DHistogram(hName,true);

                myHistos_->fill1DHistogram("h1DStats",2,std::abs(hStatsFromFile_->GetBinContent(hStatsFromFile_->FindBin(0))));
                myHistos_->fill1DHistogram("h1DStats",3,std::abs(hStatsFromFile_->GetBinContent(hStatsFromFile_->FindBin(2))));

                //hStats->Fill(1,std::abs(hStatsFromFile_->GetBinContent(hStatsFromFile_->FindBin(1))));
                //hStats->Fill(2,std::abs(hStatsFromFile_->GetBinContent(hStatsFromFile_->FindBin(3))));
        }
}
bool HZZAnalyzer::analyze(const EventProxyBase& iEvent){

        const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        getPreselectionEff(myEventProxy);
        myHistos_->fill1DHistogram("h1DStats",0);

        bool applyCuts = false;

        bool singleMuTrigger = false;
        bool muonID = true;
        unsigned int muonIDmask = (1<<6);
        std::vector<const HTTParticle *> myMuonsPlus, myMuonsMinus, myMuons;

        for(unsigned int iLepton=0; iLepton<myEventProxy.leptons->size(); ++iLepton) {
                const HTTParticle * aLepton = &myEventProxy.leptons->at(iLepton);

                int pdgId = aLepton->getProperty(PropertyEnum::PDGId);
                muonID = (int)aLepton->getProperty(PropertyEnum::muonID) & muonIDmask;

                singleMuTrigger |= aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoMu22) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoMu22_eta2p1) ||
                                   aLepton->hasTriggerMatch(TriggerEnum::HLT_IsoTkMu22_eta2p1);

                if(aLepton->getP4().Perp()<5 || std::abs(aLepton->getP4().Eta())>2.4) continue;
                if(!muonID) continue;
                if(applyCuts) {
                        if(aLepton->getProperty(PropertyEnum::combreliso)>0.35) continue;
                        //if(std::abs(aLepton->getProperty(PropertyEnum::dxy))>0.5) continue;
                        //if(std::abs(aLepton->getProperty(PropertyEnum::dz))>1) continue;
                        if(std::abs(aLepton->getProperty(PropertyEnum::SIP))>4) continue;
                }
                if(pdgId==-13) myMuonsPlus.push_back(aLepton);
                if(pdgId==13) myMuonsMinus.push_back(aLepton);
                myMuons.push_back(aLepton);
        }
        if(!singleMuTrigger) return false;
        if(applyCuts && (myMuonsPlus.size()!=2 || myMuonsMinus.size()!=2)) return false;

        float deltaR = 1;
          if(myMuonsPlus.size()>1){
          deltaR = myMuonsPlus[0]->getP4().DeltaR(myMuonsPlus[1]->getP4());
          if(deltaR<0.02) return false;
        }

        if(myMuonsMinus.size()>1){
          deltaR = myMuonsMinus[0]->getP4().DeltaR(myMuonsMinus[1]->getP4());
          if(deltaR<0.02) return false;
        }

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

        myHistos_->fill1DHistogram("h1DStats",1);

        aEntry.mass4Mu = p4ZZ.M();
        aEntry.massZ1 = p4Z1.M();
        aEntry.massZ2 = p4Z2.M();

        aEntry.nMuonsPlus = myMuonsPlus.size();
        aEntry.nMuonsMinus = myMuonsMinus.size();

        if(myMuons.size()>0){
          aEntry.muon1Pt = myMuons[0]->getP4().Perp();
          aEntry.muon1Eta = myMuons[0]->getP4().Eta();
          aEntry.muon1Phi = myMuons[0]->getP4().Phi();
          aEntry.muon1Isol = myMuons[0]->getProperty(PropertyEnum::combreliso);
          aEntry.muon1Dxy = myMuons[0]->getProperty(PropertyEnum::dxy);
          aEntry.muon1Dz = myMuons[0]->getProperty(PropertyEnum::dz);
          aEntry.muon1SIP = myMuons[0]->getProperty(PropertyEnum::SIP);
        }
        else{
          aEntry.muon1Pt = -99;
          aEntry.muon1Eta = -99;
          aEntry.muon1Phi = -99;
          aEntry.muon1Isol = -99;
          aEntry.muon1Dxy = -99;
          aEntry.muon1Dz = -99;
          aEntry.muon1SIP = -99;
        }
        if(myMuons.size()>1){
          aEntry.muon2Pt = myMuons[1]->getP4().Perp();
          aEntry.muon2Eta = myMuons[1]->getP4().Eta();
          aEntry.muon2Phi = myMuons[1]->getP4().Phi();
          aEntry.muon2Isol = myMuons[1]->getProperty(PropertyEnum::combreliso);
          aEntry.muon2Dxy = myMuons[1]->getProperty(PropertyEnum::dxy);
          aEntry.muon2Dz = myMuons[1]->getProperty(PropertyEnum::dz);
          aEntry.muon2SIP = myMuons[1]->getProperty(PropertyEnum::SIP);
        }
        else{
          aEntry.muon2Pt = -99;
          aEntry.muon2Eta = -99;
          aEntry.muon2Phi = -99;
          aEntry.muon2Isol = -99;
          aEntry.muon2Dxy = -99;
          aEntry.muon2Dz = -99;
          aEntry.muon2SIP = -99;
        }
        if(myMuons.size()>2){
          aEntry.muon3Pt = myMuons[2]->getP4().Perp();
          aEntry.muon3Eta = myMuons[2]->getP4().Eta();
          aEntry.muon3Phi = myMuons[2]->getP4().Phi();
          aEntry.muon3Isol = myMuons[2]->getProperty(PropertyEnum::combreliso);
          aEntry.muon3Dxy = myMuons[2]->getProperty(PropertyEnum::dxy);
          aEntry.muon3Dz = myMuons[2]->getProperty(PropertyEnum::dz);
          aEntry.muon3SIP = myMuons[2]->getProperty(PropertyEnum::SIP);
        }
        else{
          aEntry.muon3Pt = -99;
          aEntry.muon3Eta = -99;
          aEntry.muon3Phi = -99;
          aEntry.muon3Isol = -99;
          aEntry.muon3Dxy = -99;
          aEntry.muon3Dz = -99;
          aEntry.muon3SIP = -99;
        }
        if(myMuons.size()>3){
          aEntry.muon4Pt = myMuons[3]->getP4().Perp();
          aEntry.muon4Eta = myMuons[3]->getP4().Eta();
          aEntry.muon4Phi = myMuons[3]->getP4().Phi();
          aEntry.muon4Isol = myMuons[3]->getProperty(PropertyEnum::combreliso);
          aEntry.muon4Dxy = myMuons[3]->getProperty(PropertyEnum::dxy);
          aEntry.muon4Dz = myMuons[3]->getProperty(PropertyEnum::dz);
          aEntry.muon4SIP = myMuons[3]->getProperty(PropertyEnum::SIP);
        }
        else{
          aEntry.muon4Pt = -99;
          aEntry.muon4Eta = -99;
          aEntry.muon4Phi = -99;
          aEntry.muon4Isol = -99;
          aEntry.muon4Dxy = -99;
          aEntry.muon4Dz = -99;
          aEntry.muon4SIP = -99;
        }

        myHistos_->fill1DHistogram("h1DMassZ1",p4Z1.M());
        myHistos_->fill1DHistogram("h1DMassZ2",p4Z2.M());
        myHistos_->fill1DHistogram("h1DMass4Mu",p4ZZ.M());

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
