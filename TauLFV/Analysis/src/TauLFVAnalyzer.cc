#include <sstream>
#include <bitset>

#include "TauLFVAnalyzer.h"
#include "TauLFVHistograms.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
TauLFVAnalyzer::~TauLFVAnalyzer(){ if(myHistos_) delete myHistos_; }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* TauLFVAnalyzer::clone() const {

        TauLFVAnalyzer* clone = new TauLFVAnalyzer(name());
        clone->setHistos(myHistos_);
        return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::initialize(TDirectory* aDir,
                             pat::strbitset *aSelections){

  mySelections_ = aSelections;

  myHistos_ = new TauLFVHistograms(aDir);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::finalize(){ myHistos_->finalizeHistograms(); }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::getPreselectionEff(const EventProxyHTT & myEventProxy){

  TFile *ntupleFile = myEventProxy.getTTree()->GetCurrentFile();
  TH1F* hStatsFromFile = (TH1F*)ntupleFile->Get("hStats");
  
  std::string hName = "h1DStats_"+sampleName;
  //Fill bookkeeping histogram. Bin 0 holds number of events analyzed by the TauLFVAnalyzer.
  myHistos_->fill1DHistogram(hName, 0);
  
  TH1F *hStats = myHistos_->get1DHistogram(hName,true);
  if(!hStats) return;
  ///Bin 1 hold number of events read from dataset.
  hStats->SetBinContent(1, hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1)));
  ///Bin 3 holds number of events saved to the ntuple.
  hStats->SetBinContent(3, hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3)));

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::clear(){

  TLorentzVector emptyP4;

  aMuon1.setP4(emptyP4);
  aMuon2.setP4(emptyP4);
  aMuon3.setP4(emptyP4);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::setAnalysisObjects(const EventProxyHTT & myEventProxy){

  clear();
  
  aEvent = *myEventProxy.event;

  //std::cout<<"-------------"<<std::endl;
  
  for(auto aLepton : *myEventProxy.leptons) {
    bool isMuon = std::abs(aLepton.getPDGid())==13;
    bool hasTriggerMach = aLepton.hasTriggerMatch(TriggerEnum::HLT_DoubleMu3_Trk_Tau3mu);

    if(!isMuon || !hasTriggerMach) continue;

    if(aMuon1.getP4().E()<1E-3) aMuon1 = aLepton;
    else if(aMuon2.getP4().E()<1E-3) aMuon2 = aLepton;
    else if(aMuon3.getP4().E()<1E-3) aMuon3 = aLepton;
  }  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void TauLFVAnalyzer::fillControlHistos(const std::string & hNameSuffix, float eventWeight){
       
  const TLorentzVector & the3Mu = aMuon1.getP4() + aMuon2.getP4() + aMuon3.getP4();
	
  myHistos_->fill1DHistogram("h1DMass3Mu_"+hNameSuffix, the3Mu.M(), eventWeight);

  //std::cout<<"3Mu mass: "<<the3Mu.M()<<std::endl;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool TauLFVAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  setAnalysisObjects(myEventProxy);
  sampleName = TauLFVAnalysis::getSampleName(myEventProxy);
  getPreselectionEff(myEventProxy);

  std::string hNameSuffix = sampleName;
  float eventWeight = 1.0;

  bool has3Mu = false;
  if(aMuon1.getP4().E()>0 && aMuon2.getP4().E()>0 && aMuon3.getP4().E()>0) has3Mu =true;
  if(!has3Mu) return true;

  fillControlHistos(hNameSuffix, eventWeight);

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
