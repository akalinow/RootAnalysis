#include <iostream>
#include <cmath>

#include "OTFHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"

int nPtBins = 32;
float ptBinsTmp[33]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
  		 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
  		 160. };

const int OTFHistograms::color[6] = {kBlack, kRed, kGreen, kBlue, kMagenta, kTeal};
const float OTFHistograms::ptCutsGmt[4] = {0.1, 16, 20, 30};
const float OTFHistograms::ptCutsOtf[4] = {0.1, 14, 18, 25};


OTFHistograms::OTFHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::OTFHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::OTFHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OTFHistograms::~OTFHistograms(){ 

  std::cout<<__func__<<std::endl;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
	 for(int iCut=0;iCut<4;++iCut){
    //Gmt
    int ptCut = OTFHistograms::ptCutsGmt[iCut];
    add2DHistogram("h2DGmtPt"+std::to_string(ptCut),"",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaHit"+std::to_string(ptCut),"",23,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiHit"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaVx"+std::to_string(ptCut),"",23,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiVx"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    //Otf
    add2DHistogram("h2DOtfPt"+std::to_string(ptCut),"",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaHit"+std::to_string(ptCut),"",23,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiHit"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaVx"+std::to_string(ptCut),"",23,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiVx"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);

    if(iCut==0) continue;
    ptCut = OTFHistograms::ptCutsOtf[iCut];
    add2DHistogram("h2DGmtPt"+std::to_string(ptCut),"",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaHit"+std::to_string(ptCut),"",23,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiHit"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtEtaVx"+std::to_string(ptCut),"",23,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DGmtPhiVx"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    //Otf
    add2DHistogram("h2DOtfPt"+std::to_string(ptCut),"",150,0,150,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaHit"+std::to_string(ptCut),"",23,0.0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiHit"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfEtaVx"+std::to_string(ptCut),"",23,0,2.3,2,-0.5,1.5,file_);
    add2DHistogram("h2DOtfPhiVx"+std::to_string(ptCut),"",64,-3.2,3.2,2,-0.5,1.5,file_);

   histosInitialized_ = true;
	 }
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(int nRuns, float weight){

	plotEffPanel("Gmt");
	plotEffVsEta("Gmt");

	plotEffPanel("Otf");
	plotEffVsEta("Otf");


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* OTFHistograms::Integrate(TH1F * histoD) {

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


TH1D * OTFHistograms::DivideErr(TH1D * h1,
                 TH1D * h2,
                 const char * name,
                 const char * optErr)
{
//
// return h1/h2 with recalculated errors for "B"
//
  if (!h1) std::cout <<"DivideErr called, but histogram (h1) pointer is:"<<h1<<std::endl;
  if (!h2) std::cout <<"DivideErr called, but histogram (h2) pointer is:"<<h1<<std::endl;
  if (!h1 || !h2) return 0;
  TH1D * hout = new TH1D( *h1);
  hout->Reset();
  hout->SetName(name);
//  hout->SetTitleOffset(gStyle->GetTitleXOffset(),"x");
//  hout->SetTitleOffset(gStyle->GetTitleYOffset(),"y");
  hout->Divide(h1,h2,1.,1.,optErr);

  if (strcmp(optErr,"B")==0  || strcmp(optErr,"b")==0) {
    for (int i = 0; i<=hout->GetNbinsX()+1;i++) {
      Float_t tot   = h2->GetBinContent(i) ;
      Float_t tot_e = h2->GetBinError(i);
      Float_t eff = hout->GetBinContent(i) ;
      Float_t Err = 0.;
      if (tot > 0) {
        if (eff == 1.) eff = (tot-1)/tot; //modify efficiency as one even in numerator not fired
        Err = tot_e / tot * sqrt( eff* (1-eff) );
      }
      hout->SetBinError(i, Err);
    }
  } else {
    std::cout << "** Fig--DivideErr ** unknown option ---"<<optErr<<"---"<<std::endl;
  }

  return hout;
}





void OTFHistograms::plotEffPanel(const std::string & sysType){

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
	hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCuts[icut]);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    std::cout<<hName<<" "<<h2D<<std::endl;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
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
    TString nameCut = TString::Format("%d",(int)ptCuts[icut])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/PanelVsPt_%s.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
}


void OTFHistograms::plotEffVsEta(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("EffVsEta_%s",sysType.c_str()),
			   TString::Format("EffVsEta_%s",sysType.c_str()),
			   600,500);

  TLegend l(0.2355932,0.2819957,0.4745763,0.5856833,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TH2F *h2D = new TH2F("h2D","",8*26,0.8,1.25,2,-0.5,1.5);

  const float *ptCuts = ptCutsOtf;
  if(sysType=="Gmt" || sysType=="Rpc") ptCuts = ptCutsGmt;

  int iCut = 2;
 std::string hName = "";
  for (int iType=0; iType<3;++iType){
    hName = "h2D"+sysType+"Type" + std::to_string((int)ptCuts[iType]) + "Pt"+std::to_string((int)ptCuts[iCut]);
    TH2F* h2D = this->get2DHistogram(hName);
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->Print();
    hEff->SetStats(kFALSE);
    hEff->SetMinimum(0.0001);
    hEff->SetMaximum(1.04);
    hEff->SetMarkerStyle(21+iType);
    hEff->SetMarkerColor(color[iType]);
    hEff->SetXTitle("muon #eta");
    hEff->SetYTitle("Efficiency");
    if (iType==0)hEff->DrawCopy();
    else hEff->DrawCopy("same");
    std::string nameCut = std::to_string((int)ptCuts[iCut])+" GeV/c";
    if (iType==0) nameCut = "p_{T}^{#mu}>p_{T}^{cut} + 20 GeV";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 5 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c";
    l.AddEntry(hEff,nameCut.c_str());
  }
  l.SetHeader(TString::Format("p_{T}^{cut} = %d  GeV/c",ptCuts[iCut]).Data());
  l.DrawClone();

   ///OTF eta range used for generating patterns.
  TLine *aLine = new TLine(0,0,0,0);
  aLine->SetLineWidth(2);
  aLine->SetLineColor(2);
  aLine->DrawLine(0.83,0,0.83,1.0);
  aLine->DrawLine(1.24,0,1.24,1.0);

  c->Print(TString::Format("fig_eps/EffVsEta_%s.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVsEta_%s.png",sysType.c_str()).Data());
}

