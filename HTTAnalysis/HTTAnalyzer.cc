#include <sstream>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::HTTAnalyzer(const std::string & aName):Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTAnalyzer::~HTTAnalyzer(){

  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "HTTAnalyzer"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(&aDir, selectionFlavours_);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  
  float puWeight = myEventProxy.puWeight;
  float genWeight = myEventProxy.wevent->genevtweight();
   
  std::string sampleName = "MC";
  if(myEventProxy.wevent->sample()==0) sampleName = "Data";
  if(myEventProxy.wevent->sample()==1) {
    sampleName = "DY";
    //single event weights are huge, but all the same except sign.
    ///normalise them to 1
    genWeight/=23443.423;
  }
  if(myEventProxy.wevent->sample()==2){
    sampleName = "WJets";
    //single event weights are huge, but all the same except sign.
    ///normalise them to 1
    genWeight/=225892.45;
  }
  if(myEventProxy.wevent->sample()==3){
    sampleName = "TTbar";
    //single event weights are huge, but all the same except sign.
    ///normalise them to 1
    genWeight/=6383;
  }

  float eventWeight = puWeight*genWeight;



  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  ///This stands now for the baseline selection. 
  if(!myEventProxy.wpair->size()) return true;
  
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.wevent->npv(),eventWeight);

  if(!myEventProxy.wpair->size() ||
     !myEventProxy.wtau->size() ||
     !myEventProxy.wmu->size()) return true;

  Wpair aPair = (*myEventProxy.wpair)[0];
  Wtau aTau = (*myEventProxy.wtau)[0];
  Wmu aMuon = (*myEventProxy.wmu)[0];
  
  if(aPair.diq() == 1){
     myHistos_->fill1DHistogram("h1DMassSV"+sampleName+"qcdselSS", aPair.svfit() ,eventWeight);
     myHistos_->fill1DHistogram("h1DMassVis"+sampleName+"qcdselSS", aPair.m_vis() ,eventWeight);
     myHistos_->fill1DHistogram("h1DMassTrans"+sampleName+"qcdselSS", aMuon.mt() ,eventWeight);
     myHistos_->fill1DHistogram("h1DPtMuon"+sampleName+"qcdselSS", aMuon.pt() ,eventWeight);
     myHistos_->fill1DHistogram("h1DEtaMuon"+sampleName+"qcdselSS", aMuon.eta() ,eventWeight);
     myHistos_->fill1DHistogram("h1DPtTau"+sampleName+"qcdselSS", aTau.pt() ,eventWeight);
     myHistos_->fill1DHistogram("h1DEtaTau"+sampleName+"qcdselSS", aTau.eta() ,eventWeight);
     myHistos_->fill1DHistogram("h1DNPV"+sampleName+"qcdselSS", myEventProxy.wevent->npv() ,eventWeight);

     myHistos_->fill1DHistogram("h1DPhiMuon"+sampleName+"qcdselSS", aMuon.phi() ,eventWeight);
     myHistos_->fill1DHistogram("h1DPhiTau"+sampleName+"qcdselSS", aTau.phi() ,eventWeight);
     myHistos_->fill1DHistogram("h1DMtTau"+sampleName+"qcdselSS", aTau.mt() ,eventWeight);
     myHistos_->fill1DHistogram("h1DIsoMuon"+sampleName+"qcdselSS", aMuon.iso() ,eventWeight);
  }


  ///Fill SVfit and visible masses
  myHistos_->fill1DHistogram("h1DMassSV"+sampleName,aPair.svfit(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+sampleName,aPair.m_vis(),eventWeight);

  ///Fill muon pt and eta
  myHistos_->fill1DHistogram("h1DMassTrans"+sampleName,aMuon.mt(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtMuon"+sampleName,aMuon.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+sampleName,aMuon.eta(),eventWeight);

  ///Fill tau pt and eta
  myHistos_->fill1DHistogram("h1DPtTau"+sampleName,aTau.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+sampleName,aTau.eta(),eventWeight);

  myHistos_->fill1DHistogram("h1DIsoMuon"+sampleName,aMuon.iso(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMuon"+sampleName,  aMuon.phi(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+sampleName, aTau.phi() ,eventWeight);
  myHistos_->fill1DHistogram("h1DMtTau"+sampleName,  aTau.mt() ,eventWeight);

  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

