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

  ///Load ROOT file with PU histograms.
  std::string filePath = "./RootAnalysis_Weights.root";//FIX ME
  puFile_ = new TFile(filePath.c_str());
  
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
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return "Data";
  if(myEventProxy.wevent->sample()==1) return "DYJets";
  if(myEventProxy.wevent->sample()==2) return "WJets";
  if(myEventProxy.wevent->sample()==3) return "TTbar";

  return "Unknown";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getPUWeight(const EventProxyHTT & myEventProxy){

  ///Load histogram only once,later fetch it from vector<TH1F*>
  ///At the same time divide the histogram to get the weight.
  ///First load Data PU
  if(!hPUVec_.size())  hPUVec_.resize(64);

  if(!hPUVec_[myEventProxy.wevent->sample()]){
    std::string hName = "Summary/h1DNPV"+getSampleName(myEventProxy);
    std::cout<<"Loading PU histogram hName: "<<hName<<std::endl;
    TH1F *hPUData = (TH1F*)puFile_->Get("Summary/h1DNPVData");
    TH1F *hPUSample = (TH1F*)puFile_->Get(hName.c_str());
    ///Normalise both histograms.
    hPUData->Scale(1.0/hPUData->Integral(0,hPUData->GetNbinsX()+1));
    hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
    ///
    hPUData->SetDirectory(0);
    hPUSample->SetDirectory(0);
    hPUData->Divide(hPUSample);
    hPUData->SetName(("h1DPUWeight"+getSampleName(myEventProxy)).c_str());
    ///To get uniform treatment put weight=1.0 for under/overlow bins of
    ///data PU, as nPU for data has a dummy value.
    if(getSampleName(myEventProxy)=="Data"){
      hPUData->SetBinContent(0,1.0);
      hPUData->SetBinContent(hPUData->GetNbinsX()+1,1.0);
    }
    hPUVec_[myEventProxy.wevent->sample()] =  hPUData;
  }

  return  hPUVec_[myEventProxy.wevent->sample()]->GetBinContent(myEventProxy.wevent->npu());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getGenWeight(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return 1.0;
  if(myEventProxy.wevent->sample()==1) return myEventProxy.wevent->genevtweight()/23443.423;  
  if(myEventProxy.wevent->sample()==2) return myEventProxy.wevent->genevtweight()/225892.45;  
  if(myEventProxy.wevent->sample()==3) return myEventProxy.wevent->genevtweight()/6383;

  return 1;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillControlHistos(Wevent & aEvent, 
				    Wpair & aPair, Wtau & aTau, Wmu & aMuon,
				    Wjet & aJet,
				    float eventWeight,
				    const std::string & hNameSuffix){

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+hNameSuffix,aEvent.npv(),eventWeight);

  ///Fill SVfit and visible masses
  myHistos_->fill1DHistogram("h1DMassSV"+hNameSuffix,aPair.svfit(),eventWeight);
  myHistos_->fill1DHistogram("h1DMassVis"+hNameSuffix,aPair.m_vis(),eventWeight);
  
  ///Fill muon pt and eta
  myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix,aMuon.mt(),eventWeight);
  myHistos_->fill1DHistogram("h1DPtMuon"+hNameSuffix,aMuon.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaMuon"+hNameSuffix,aMuon.eta(),eventWeight);

  ///Fill tau pt and eta
  myHistos_->fill1DHistogram("h1DPtTau"+hNameSuffix,aTau.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaTau"+hNameSuffix,aTau.eta(),eventWeight);

  myHistos_->fill1DHistogram("h1DIsoMuon"+hNameSuffix,aMuon.iso(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiMuon"+hNameSuffix,  aMuon.phi(),eventWeight);
  myHistos_->fill1DHistogram("h1DPhiTau"+hNameSuffix, aTau.phi() ,eventWeight);

  ///Fill jets info
  myHistos_->fill1DHistogram("h1DPtLeadingJet"+hNameSuffix,aJet.pt(),eventWeight);
  myHistos_->fill1DHistogram("h1DEtaLeadingJet"+hNameSuffix,aJet.eta(),eventWeight);

  if(aJet.bjet()){
    myHistos_->fill1DHistogram("h1DPtLeadingBJet"+hNameSuffix,aJet.pt(),eventWeight);
    myHistos_->fill1DHistogram("h1DEtaLeadingBJet"+hNameSuffix,aJet.pt(),eventWeight);
  }

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);


  std::string sampleName = getSampleName(myEventProxy);  
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);
  float eventWeight = puWeight*genWeight;
   
  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  ///To spedup processing we load only event with at least one tau pair.
  ///Have to cast away const fromthe event. 
  ///WARNING: needs check with any new ROOT version, as SetBranchStatus
  ///behaviour may change.
  EventProxyHTT & myEventProxyMod = const_cast<EventProxyHTT&>(myEventProxy);
  if(myEventProxy.wpair->size()){
    //myEventProxyMod.enableBranches();
    //myEventProxyMod.reloadEvent();
  }  
  /////////////////////////////////////////////////////////////////////////////


  if(!myEventProxy.wpair->size() || !myEventProxy.wtau->size() || !myEventProxy.wmu->size()) return true;
  

  Wevent aEvent = *myEventProxy.wevent;
  Wpair aPair = (*myEventProxy.wpair)[0];
  Wtau aTau = (*myEventProxy.wtau)[0];
  Wmu aMuon = (*myEventProxy.wmu)[0];
  Wjet aJet;
  if(!myEventProxy.wjet->size()) aJet = (*myEventProxy.wjet)[0];

  ///This stands for core selection, that is common to all regions.
  //TEST if(!myEventProxy.wpair->size() || aTau.pt()<30 || aMuon.pt()<20) return true;
  if(!myEventProxy.wpair->size() || aTau.pt()<25 || aMuon.pt()<20) return true;///Loosen pt cuts

  ///Note: parts of the signal/control region selection are applied in the following code.
  ///FIXME AK: this should be made in a more clear way.
  bool baselineSelection = aPair.diq()==-1 && aMuon.mt()<40 && aMuon.iso()<0.1;
  bool wSelection = aMuon.mt()>60 && aMuon.iso()<0.1;
  bool qcdSelectionSS = aPair.diq()==1;
  bool qcdSelectionOS = aPair.diq()==-1;

  ///Histograms for the baseline selection  
  std::string hNameSuffix = sampleName;
  if(baselineSelection) fillControlHistos(aEvent, aPair, aTau, aMuon, aJet, eventWeight, hNameSuffix);

  ///Histograms for the QCD control region
  if(qcdSelectionSS){
    hNameSuffix = sampleName+"qcdselSS";
    ///SS ans OS isolation histograms are filled only for mT<40 to remove possible contamnation
    //from TT in high mT region.
    if(aMuon.mt()<40) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
    ///Fill SS histos in signal mu isolation region. Those histograms
    ///provide shapes for QCD estimate in signal region and in various control regions.
    ///If control region has OS we still use SS QCD estimate.
    if(aMuon.mt()<40 && aMuon.iso()<0.1) fillControlHistos(aEvent, aPair, aTau, aMuon, aJet, eventWeight, hNameSuffix);
    if(aMuon.mt()>60){
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselSS",aMuon.mt(),eventWeight);    
      myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"wselOS",aMuon.mt(),eventWeight);    
    }
  }
  ///Make QCD shape histograms for specific selection.
  ///Using the same SS/OS scaling factor for now.    
  if(qcdSelectionOS){
    hNameSuffix = sampleName+"qcdselOS";
    if(aMuon.mt()<40) myHistos_->fill1DHistogram("h1DIso"+hNameSuffix,aMuon.iso(),eventWeight);
  }

  ///Histograms for the WJet control region. 
  ///Selection is split into SS and OS regions.
  if(wSelection){
    hNameSuffix = sampleName+"wsel";
    if(aPair.diq()==-1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"OS",aMuon.mt(),eventWeight);
    if(aPair.diq()== 1) myHistos_->fill1DHistogram("h1DMassTrans"+hNameSuffix+"SS",aMuon.mt(),eventWeight);
  }

  ///Histograms for the tt control region


  ///Disable branches before loading next event.
  //myEventProxyMod.disableBranches();
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

