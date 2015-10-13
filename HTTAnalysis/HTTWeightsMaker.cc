
#include <sstream>

#include "HTTWeightsMaker.h"
#include "HTTWeightHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTWeightsMaker::HTTWeightsMaker(const std::string & aName):Analyzer(aName){

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTWeightsMaker::~HTTWeightsMaker(){

  if(hPU) delete hPU;
  if(hDatasetPU) delete hDatasetPU;
  if(hPUWeights) delete hPUWeights;
  if(puFile) delete puFile;
  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  mySelections_ = aSelections;

  std::string filePath = "/home/akalinow/scratch/EclipseProjects/RootAnalysis/HTTAnalysis/RootAnalysis_PU.root";
  puFile = new TFile(filePath.c_str());
  hPU = (TH1F*)puFile->Get("Summary/h1DNPVDATA");
  if(hPU){
    hPU->SetName("hPU");
    hPU->Scale(1.0/hPU->Integral(0,hPU->GetNbinsX()+1));
    hPUWeights = (TH1F*)hPU->Clone("hPUWeights");
    hPUWeights->SetName("hPUWeights");
    hPUWeights->Reset();
  }
  hDatasetPU = 0;
  
  ///The histograms for this analyzer will be saved into "HTTWeightsMaker"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTWeightHistograms(&aDir, selectionFlavours_);  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTWeightsMaker::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  
  puWeight = 1.0;
  float genWeight = myEventProxy.genWeight;
  float eventWeight = puWeight*genWeight;

  std::string sampleName = "MC";
  if(myEventProxy.run>1) sampleName = "DATA";
  if(sampleName=="DATA") eventWeight = 1.0;

  std::string hName = "h1DNPV"+sampleName;
  
  if(hPU && (!hDatasetPU || hDatasetPU->GetName()!=hName)){
    hDatasetPU = (TH1F*)puFile->Get(("Summary/"+hName).c_str());
    hDatasetPU->Scale(1.0/hDatasetPU->Integral(0,hDatasetPU->GetNbinsX()+1));
    hPUWeights->Divide(hPU,hDatasetPU);    
  }
    
  if(hPUWeights) puWeight = hPUWeights->GetBinContent(hPUWeights->FindBin(myEventProxy.npv));

  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.npv,eventWeight);
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::addBranch(TTree *tree){

  tree->Branch("PUWeight",&puWeight);
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
