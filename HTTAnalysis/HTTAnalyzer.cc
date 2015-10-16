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

  float eventWeight = puWeight*genWeight;

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  ///This stands now for the baseline selection. 
  if(!myEventProxy.wpair->size()) return true;
  
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.wevent->npv(),eventWeight);

  Wpair aPair = (*myEventProxy.wpair)[0];
  ///Fill SVfit mass
  myHistos_->fill1DHistogram("h1DSVfit"+sampleName,aPair.svfit(),eventWeight);
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

