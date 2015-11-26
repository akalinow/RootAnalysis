
#include <sstream>

#include "HTTAnalyzer.h"
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
  if(puFile) delete puFile;
  if(myHistos_) delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::initialize(TFileDirectory& aDir,
				 pat::strbitset *aSelections){

  mySelections_ = aSelections;

  std::string filePath = "./MyDataPileupHistogram.root";//FIX ME
  puFile = new TFile(filePath.c_str());
  hPU = (TH1F*)puFile->Get("pileup");
  
  ///The histograms for this analyzer will be saved into "HTTWeightsMaker"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTWeightHistograms(&aDir, selectionFlavours_);  
  ///Substitute the default PU histogram ranges, with ranges taken from
  ///reference PU histogram.
  TH1F *h = myHistos_->get1DHistogram("h1DNPVTemplate",true);
  TDirectory *tmpDir= h->GetDirectory();
  std::cout<<"path: "<<tmpDir->GetName()<<std::endl;
  *h = TH1F("h1DNPVTemplate","",
	    hPU->GetNbinsX(),
	    hPU->GetXaxis()->GetXmin(),
	    hPU->GetXaxis()->GetXmax());
  h->SetDirectory(tmpDir);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTWeightsMaker::finalize(){ 

  ///Put the correct contents of the DATA PU histogram.
  myHistos_->fill1DHistogram("h1DNPVData",-999,0.0);
  TH1F *hDataPU = myHistos_->get1DHistogram("h1DNPVData",true);
  hDataPU->Reset();
  hDataPU->Add(hPU);

  myHistos_->finalizeHistograms(0,1.0);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTWeightsMaker::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);
  
  ///For making the PU weights we do not use generator level weights,
  ///as PU and generator level weights should be independent.
  float eventWeight = 1.0;

  std::string sampleName  = HTTAnalyzer::getSampleName(myEventProxy);

  ///Do not fill PU histo for the data.
  ///Data histo content is copied from external  histogram 
  ///with nPU. 
  if(sampleName=="Data") return true;

  std::string hName = "h1DNPV"+sampleName;
  
  ///Fill histograms with number of PV.
  myHistos_->fill1DHistogram("h1DNPV"+sampleName,myEventProxy.wevent->npu(),eventWeight);
  
  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
