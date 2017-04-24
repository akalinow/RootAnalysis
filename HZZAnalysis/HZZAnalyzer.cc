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
void HZZAnalyzer::addBranch(TTree *tree){ /*tree->Branch("muonPt",&muonPt);*/
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HZZAnalyzer::analyze(const EventProxyBase& iEvent){

        const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

        bool muonID = true;
	unsigned int muonIDmask = (1<<6);
        std::vector<HTTParticle> myMuonsPlus, myMuonsMinus, myMuons;
        for(auto aLepton : *myEventProxy.leptons) {
                int pdgId = aLepton.getProperty(PropertyEnum::PDGId);
                if(std::abs(pdgId)!=13) continue;
		muonID = (int)aLepton.getProperty(PropertyEnum::muonID) & muonIDmask;
		if(!muonID) continue;		
                if(aLepton.getP4().Perp()<5 || std::abs(aLepton.getP4().Eta())>2.4) continue;
                if(aLepton.getProperty(PropertyEnum::combreliso)>0.35) continue;
                if(pdgId==-13) myMuonsPlus.push_back(aLepton);
                if(pdgId==13) myMuonsMinus.push_back(aLepton);
		myMuons.push_back(aLepton);
        }

        if(myMuonsPlus.size()!=2 || myMuonsMinus.size()!=2) return true;

        unsigned int nPt10 = 0, nPt20 = 0;

        for(auto it:myMuonsPlus) {
                if(it.getP4().Perp()>10) ++nPt10;
                if(it.getP4().Perp()>20) ++nPt20;
        }
        for(auto it:myMuonsMinus) {
                if(it.getP4().Perp()>10) ++nPt10;
                if(it.getP4().Perp()>20) ++nPt20;
        }

        if(nPt10<2 || nPt20<1) return true;

        TLorentzVector p4Z1, p4Z2;
        TLorentzVector p4Tmp;

        TLorentzVector p4ZZ;
        for(auto it:myMuonsPlus) p4ZZ+=it.getP4();
        for(auto it:myMuonsMinus) p4ZZ+=it.getP4();

        float deltaMass = 9999;
        float deltaTmp = 9999;
	float deltaR;

        deltaR = myMuonsPlus[0].getP4().DeltaR(myMuonsPlus[1].getP4());
	if(deltaR<0.02) return true;

	deltaR = myMuonsMinus[0].getP4().DeltaR(myMuonsMinus[1].getP4());
	if(deltaR<0.02) return true;
		
        for(auto aMuonPlus: myMuonsPlus){
          for(auto aMuonMinus: myMuonsMinus){
            deltaR = aMuonPlus.getP4().DeltaR(aMuonMinus.getP4());
            if(deltaR<0.02) return true;
            p4Tmp = aMuonPlus.getP4() + aMuonMinus.getP4();
            if(p4Tmp.M()<4) return true;
            deltaTmp = std::abs(p4Tmp.M() - 91.1896);
            if(deltaTmp<deltaMass) {
              deltaMass = deltaTmp;
              p4Z1 = p4Tmp;
              p4Z2 = p4ZZ - p4Z1;
            }
          }
        }

        if(p4Z1.M()>120 || p4Z1.M()<12) return true;
        if(p4Z2.M()>120 || p4Z2.M()<12) return true;
        if(p4Z1.M()<40) return true;

        myHistos_->fill1DHistogram("h1DMassZ1",p4Z1.M());
        myHistos_->fill1DHistogram("h1DMassZ2",p4Z2.M());
        myHistos_->fill1DHistogram("h1DMass4Mu",p4ZZ.M());

	myHistos_->fill1DHistogram("h1DMassPtMu1",myMuons[0].getP4().Perp());
	myHistos_->fill1DHistogram("h1DMassPtMu2",myMuons[1].getP4().Perp());
	myHistos_->fill1DHistogram("h1DMassPtMu3",myMuons[2].getP4().Perp());
	myHistos_->fill1DHistogram("h1DMassPtMu4",myMuons[3].getP4().Perp());

	

        return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
