#include <cstdlib> 

#include "TreeAnalyzer.h"
#include "OTFAnalyzer.h"
#include "EventProxyOTF.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

 const int OTFAnalyzer::color[6] = {kBlack, kRed, kGreen, kBlue, kMagenta, kTeal};
 const float OTFAnalyzer::ptCutsGmt[4] = {0.1, 16, 20, 30};
 const float OTFAnalyzer::ptCutsOtf[4] = {0.1, 14, 18, 25};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::OTFAnalyzer(const std::string & aName):Analyzer(aName){

  clear();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
OTFAnalyzer::~OTFAnalyzer(){

  delete myHistos_;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::initialize(TFileDirectory& aDir,
				       pat::strbitset *aSelections){

  mySelections_ = aSelections;

  ///The histograms for this analyzer will be saved into "TestHistos"
  ///directory of the ROOT file
  myHistos_ = new OTFHistograms(&aDir, selectionFlavours_);

  registerCuts();

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::finalize(){ }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::registerCuts(){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){
    mySelections_->push_back("HLT"+selectionFlavours_[i]);
  /*
 ///Register tree variables
 for(unsigned int i=0;i<triggerItemNames_.size();++i){
   treeVariables_[triggerItemNames_[i]] = -999.0;
 }
 */
 }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::fillTurnOnCurve(float & ptCut, const std::string & sysType){

	bool qualityCut = true;

std::vector<L1Obj> & myL1Coll = theEvent->l1ObjectsGmt;
std::string hName = "h2DGmt";

//if(sysType=="Rpc") myL1Coll = theEvent->l1ObjectsRpc;
if(sysType=="Otf") {
	myL1Coll = theEvent->l1ObjectsOtf;
	qualityCut = myL1Coll.size() && myL1Coll[0].q!=103 &&
		         myL1Coll[0].q!=104 && myL1Coll[0].q!=105 &&
		         myL1Coll[0].q%100>3;
	hName = "h2DOtf";
	}

   bool pass = (myL1Coll.size() && myL1Coll[0].pt>=ptCut && qualityCut) || !myL1Coll.size();


   std::string tmpName = hName+"Pt";
   myHistos_->fill2DHistogram(tmpName,theEvent->pt, pass);
   tmpName = hName+"EtaHit";
   myHistos_->fill2DHistogram(tmpName,theEvent->etaHit, pass);
   tmpName = hName+"PhiHit";
   myHistos_->fill2DHistogram(tmpName,theEvent->phiHit, pass);
   tmpName = hName+"EtaVx";
   myHistos_->fill2DHistogram(tmpName,theEvent->eta, pass);
   tmpName = hName+"PhiVx";
   myHistos_->fill2DHistogram(tmpName,theEvent->phi, pass);

}


bool OTFAnalyzer::analyze(const EventProxyBase& iEvent){

  clear();

  eventWeight_ = 1.0;
  /////////////////
  const EventProxyOTF & myEvent = static_cast<const EventProxyOTF&>(iEvent);
  theEvent = myEvent.events;

  std::string sysType = "Gmt";
  float ptCut = 20;

  sysType = "Gmt";
  fillTurnOnCurve(ptCut,sysType);
  sysType = "Otf";
  fillTurnOnCurve(ptCut,sysType);



  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool  OTFAnalyzer::checkSelections(const std::string & type){

  bool decision = false;

  mySelections_->set("HLT"+type,decision);

 return decision;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void  OTFAnalyzer::addBranch(TTree *tree){

 for(unsigned int i=0;i<selectionFlavours_.size();++i){ 
   std::map<std::string,float>::const_iterator CI = treeVariables_.begin();
    for(;CI!=treeVariables_.end();++CI) tree->Branch(CI->first.c_str(),&treeVariables_[CI->first]);
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void OTFAnalyzer::clear(){

  ///Clear variables
  std::map<std::string,float>::iterator it=treeVariables_.begin();
  for(;it!=treeVariables_.end();++it) it->second = -999;
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
TH1F* OTFAnalyzer::Integrate(TH1F * histoD) {

  TH1F* histoI = (TH1F*)histoD->Clone("hIntegrated");

  Double_t * cont = new Double_t [histoD->GetNbinsX()+2];  //with under+overflow
  Double_t * errs = new Double_t [histoD->GetNbinsX()+2];  //with under+overflow
  histoI->Reset();

  // bin=0 underf
  // bin 1-GetNbinsX() -conten
  // bin GetNbinsX()+1 overflow

  for (Int_t i = 0; i <= histoD->GetNbinsX()+1; i++) {
    cont[i] = (Double_t)(histoD->GetBinContent(i));
    errs[i] = (Double_t)(histoD->GetBinError(i));
  }
  Double_t sum=0.;
  Double_t sume2=0.;
  for (Int_t i = histoD->GetNbinsX()+1; i >= 0; i--) {
    sum+=cont[i];
    sume2+=errs[i]*errs[i];
    histoI->SetBinContent(i,sum);
    histoI->SetBinError(i,sqrt(sume2));
   }
  delete [] cont;
  delete [] errs;
  return histoI;
}

void OTFAnalyzer::plotEffPanel(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("EffVsPt_%s",sysType.c_str()),
			   TString::Format("EffVsPt_%s",sysType.c_str()),
			   460,500);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(0,1);

  TString hName("");
  const float *ptCuts = ptCutsOtf;
  if(sysType=="Gmt") ptCuts = ptCutsGmt;

  for (int icut=0; icut <=3;++icut){

	hName = "";
    TH2F* h2D = myHistos_->get2DHistogram(hName.Data());
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    //TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    TH1D *hEff = 0;


    hEff->SetStats(kFALSE);
    hEff->SetMinimum(0.0001);
    hEff->SetMaximum(1.04);
    hEff->GetXaxis()->SetRange(4,100);
    hEff->SetMarkerStyle(21+icut);
    hEff->SetMarkerColor(color[icut]);
    hEff->SetXTitle("muon p_{T} [GeV/c]");
    hEff->SetYTitle("Efficiency");
    if (icut==0)hEff->DrawCopy();
    else hEff->DrawCopy("same");
    TString nameCut = TString::Format("%d",ptCuts[icut])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/PanelVsPt_%s.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
}
