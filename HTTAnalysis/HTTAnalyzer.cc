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

  //myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

  if(myEventProxy.wevent->sample()==0) return "Data";
  if(myEventProxy.wevent->sample()==1) return "DY";
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
bool HTTAnalyzer::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  std::string sampleName = getSampleName(myEventProxy);  
  float puWeight = getPUWeight(myEventProxy);
  float genWeight = getGenWeight(myEventProxy);

  //puWeight = 1.0;
  //genWeight = 1.0;
  float eventWeight = puWeight*genWeight;

  if(myEventProxy.wevent->npv()<2) return true; //Temporary fix against bad PU weights for npv==1
   
  //Fill bookkeeping histogram. Bin 1 holds sum of weights.
  myHistos_->fill1DHistogram("h1DStats"+sampleName,1,eventWeight);

  if(!myEventProxy.wpair->size() ||
     !myEventProxy.wtau->size() ||
     !myEventProxy.wmu->size()) return true;

  Wpair aPair = (*myEventProxy.wpair)[0];
  Wtau aTau = (*myEventProxy.wtau)[0];
  Wmu aMuon = (*myEventProxy.wmu)[0];

  ///This stands now for the baseline selection. 
  ///This stands now for the baseline selection. 
  if(!myEventProxy.wpair->size() || aTau.pt()<30 || aMuon.pt()<20 || aMuon.iso()>0.1) return true;

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.wevent->npv(),eventWeight);
  
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

