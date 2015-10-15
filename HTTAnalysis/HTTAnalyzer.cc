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
  float genWeight = myEventProxy.sampleWeight;
  float eventWeight = puWeight*genWeight;

  std::string sampleName = "MC";
  if(genWeight==1) sampleName = "Data";

  if(fabs(fabs(myEventProxy.genDecay/24.0)-13)<1E-5 ||
     fabs(fabs(myEventProxy.genDecay/24.0)-15)<1E-5){
     sampleName = "WJets";
  }
  if(fabs(fabs(myEventProxy.genDecay/23.0)-13)<1E-5 ||
     fabs(fabs(myEventProxy.genDecay/23.0)-15)<1E-5){
    sampleName = "DY";
  }
  if(fabs(genWeight-0.0326672)<1E-5) sampleName = "TT";
  if(fabs(genWeight-0.0031596)<1E-5) sampleName = "Other";

  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  bool baselineSelection = myEventProxy.ptL1>20 && myEventProxy.isPFMuon && myEventProxy.isTightMuon &&
    myEventProxy.ptL2>20 && myEventProxy.muFlag!=1 && myEventProxy.vetoEvent==0 &&
    myEventProxy.tightestHPSMVAWP>=0 && myEventProxy.combRelIsoLeg1DBetav2<0.1 &&
    myEventProxy.diTauCharge==0 && myEventProxy.MtLeg1MVA<40 && 
    myEventProxy.pairIndex<1 && myEventProxy.HLTx==1 && (myEventProxy.run>=163269 || myEventProxy.run==1) &&
    myEventProxy.HLTmatch==1;

  if(!baselineSelection) return true;
    
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.numPV,eventWeight);

  ///Fill SVfit mass
  myHistos_->fill1DHistogram("h1DSVfit"+sampleName,myEventProxy.diTauNSVfitMass,eventWeight);

  ///Fill transverse mass
  myHistos_->fill1DHistogram("h1DMt"+sampleName,myEventProxy.MtLeg1MVA,eventWeight);
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

