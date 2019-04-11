#include <iostream>
#include <cmath>

#include "OMTFHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TGaxis.h"

#include "utilsL1RpcStyle.h"

int nPtBins = 32;
const float OMTFHistograms::ptBins[33]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
  		 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
  		 160. };

const int OMTFHistograms::color[6] = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
//Single mu
const int OMTFHistograms::ptCutsGmt[4] =     {0, 14, 19, 20};
const int OMTFHistograms::ptCutsOMTF[4] =     {0, 14, 19, 20};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::OMTFHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::OMTFHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::OMTFHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::~OMTFHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string OMTFHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";

  if(name.find("DeltaEta")!=std::string::npos) templateName = "h1DDeltaEtaTemplate";
  if(name.find("DeltaPhi")!=std::string::npos) templateName = "h1DDeltaPhiTemplate";

  if(name.find("Pt")!=std::string::npos) templateName = "h2DPtTemplate";
  if(name.find("HighPt")!=std::string::npos) templateName = "h2DHighPtTemplate";
  if(name.find("PtRecVsPtGen")!=std::string::npos) templateName = "h2DPtVsPtTemplate";
  
  if(name.find("EtaHit")!=std::string::npos) templateName = "h2DEtaHitTemplate";
  if(name.find("PhiHit")!=std::string::npos) templateName = "h2DPhiHitTemplate";
  if(name.find("EtaVx")!=std::string::npos) templateName = "h2DEtaVxTemplate";
  if(name.find("PhiVx")!=std::string::npos) templateName = "h2DPhiVxTemplate";
  if(name.find("Quality")!=std::string::npos) templateName = "h2DQualityTemplate";
  if(name.find("RateTot")!=std::string::npos) templateName = "h2DRateTotTemplate";
  if(name.find("RateVsEta")!=std::string::npos) templateName = "h2DRateVsEtaTemplate";
  if(name.find("RateVsPt")!=std::string::npos) templateName = "h2DRateVsPtTemplate";
  if(name.find("RateVsQuality")!=std::string::npos) templateName = "h2DRateVsQualityTemplate";
  if(name.find("DeltaPhi")!=std::string::npos) templateName = "h2DDeltaPhiTemplate";
  if(name.find("DeltaPt")!=std::string::npos) templateName = "h2DDeltaPtTemplate";
  if(name.find("GhostsVsProcessor")!=std::string::npos) templateName = "h2DGhostsVsProcessorTemplate";

  if(name.find("LLH")!=std::string::npos) templateName = "h1DLLHTemplate";
  if(name.find("HitsPattern")!=std::string::npos) templateName = "h1DHitsPatternTemplate";
  

  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

 //Make template histos
 add1DHistogram("h1DDeltaEtaTemplate","",11,-0.83,0.83,file_);

 add1DHistogram("h1DDeltaPhiTemplate","",5*32,-M_PI,M_PI,file_);

 ///Efficiency histos
 add2DHistogram("h2DPtTemplate","",150,0,150,2,-0.5,1.5,file_);
 add2DHistogram("h2DHighPtTemplate","",50,50,550,2,-0.5,1.5,file_);

 add2DHistogram("h2DPtVsPtTemplate","",70,0,140,70,0,140,file_);

 add2DHistogram("h2DEtaHitTemplate","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHitTemplate","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaVxTemplate","",40,-1.6,1.6,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiVxTemplate","",4*32,-3.2,3.2,2,-0.5,1.5,file_);

 add2DHistogram("h2DQualityTemplate","",201,-0.5,200.5,2,-0.5,1.5,file_);

 //Rate histos
 add2DHistogram("h2DRateTotTemplate","",400,1,201,142,0,142,file_);
 add2DHistogram("h2DRateVsEtaTemplate","",400,1,201,32*2,-1.6,1.6,file_);

 add2DHistogram("h2DDeltaPhiTemplate","",30,-1,1,2,-0.5,1.5,file_);
 add2DHistogram("h2DDeltaPtTemplate","",21,-0.5,20.5,2,-0.5,1.5,file_);

 add2DHistogram("h2DRateVsPtTemplate","",400,1,201,100,0,50,file_);
 
 add2DHistogram("h2DRateVsQualityTemplate","",400,1,201,201,-0.5,200.5,file_);
 add2DHistogram("h2DGhostsVsProcessorTemplate","",6,-0.5,5.5,5,-0.5,4.5,file_);

 //Likelihood histos
 add1DHistogram("h1DLLHTemplate","",60,60,300,file_);
 add1DHistogram("h1DHitsPatternTemplate","",21,-0.5,20.5,file_);

 
 histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::finalizeHistograms(){

  AnalysisHistograms::finalizeHistograms();

  utilsL1RpcStyle()->cd();

  plotRate("Tot");
  //plotRate("VsEta");
  //plotRate("VsPt");
  plotRate("VsQuality");

  plotEffPanel("OMTF");
  plotEffPanel("kBMTF");
  plotEffPanel("BMTF");

  bool doHigh = true;
  plotEffPanel("OMTF", doHigh);
  plotEffPanel("BMTF", doHigh);
  plotEffPanel("kBMTF", doHigh);
  plotEffVsEta("OMTF");
  plotEffVsEta("EMTF");
  plotEffVsEta("BMTF");
  plotEffVsEta("kBMTF");
  plotEffVsVar("OMTF","EtaVx");
  plotEffVsVar("OMTF","PhiVx");
  plotSingleHistogram("h2DOMTFPtRecVsPtGen");
  plotSingleHistogram("h2DOMTFPtRecVsPtGenLow");
  plotSingleHistogram("h2DkBMTFPtRecVsPtGen");

  plotOMTFVsOther(16,"kBMTF");
  plotOMTFVsOther(18,"kBMTF");
  plotOMTFVsOther(19,"kBMTF");
  plotOMTFVsOther(20,"kBMTF");
  plotOMTFVsOther(20,"BMTF");

  plotOMTFVsOther(16,"EMTF");
  plotOMTFVsOther(18,"EMTF");
  plotOMTFVsOther(19,"EMTF");
  plotOMTFVsOther(20,"EMTF");

  plotSingleHistogram("h1DLLH_Low");
  plotSingleHistogram("h1DLLH_High");
  plotLLH();

  plotSingleHistogram("h1DHitsPattern_Low");
  plotSingleHistogram("h1DHitsPattern_High");

  plotSingleHistogram("h1DDeltaEta_Low");
  plotSingleHistogram("h1DDeltaEta_High");

  plotSingleHistogram("h1DHitsPattern_Low_RefLayer");
  plotSingleHistogram("h1DHitsPattern_High_RefLayer");
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* OMTFHistograms::Integrate(TH1 * histoD) {

  TH1* histoI = (TH1*)histoD->Clone("hIntegrated");

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
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1D * OMTFHistograms::DivideErr(TH1D * h1,
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
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::plotEffPanel(const std::string & sysType, bool doHigh){

  TCanvas* c = new TCanvas(TString::Format("EffVsPt_%s",sysType.c_str()),
			   TString::Format("EffVsPt_%s",sysType.c_str()),
			   460,500);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  //c->SetLogx(1);
  c->SetGrid(0,1);

  TString hName("");
  const int *ptCuts = ptCutsOMTF;
  if(sysType.find("Gmt")!=std::string::npos ||
     sysType=="Rpc" ||
     sysType=="Other") ptCuts = ptCutsGmt;

  for (int icut=0; icut <=3;++icut){
    float ptCut = OMTFHistograms::ptBins[ptCuts[icut]];
    hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
    if(doHigh) hName = "h2D"+sysType+"HighPt"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);    
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->SetStats(kFALSE);
    hEff->SetMinimum(0.0001);
    hEff->SetMaximum(1.04);
    hEff->GetXaxis()->SetRange(1,50);
    hEff->SetMarkerStyle(21+icut);
    hEff->SetMarkerColor(color[icut]);
    hEff->SetXTitle("gen. muon p_{T} [GeV/c]");
    hEff->SetYTitle("Efficiency");
    if (icut==0)hEff->DrawCopy();
    else hEff->DrawCopy("same");
    TString nameCut = TString::Format("%d", (int)OMTFHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
  }
  l.DrawClone();
  if(!doHigh) c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
  else  c->Print(TString::Format("fig_png/PanelVsHighPt_%s.png",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotVar(const std::string & sysType,
			    const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("Var%s_%s",varName.c_str(),sysType.c_str()),
			   TString::Format("Var%s_%s",varName.c_str(),sysType.c_str()),
			   460,500);


  TString hName = "h1D"+varName+sysType;
  TH1F* h1D = this->get1DHistogram(hName.Data());
  h1D->SetLineWidth(3);
  h1D->SetStats(kFALSE);

  c->SetLogy();

  h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
  h1D->SetXTitle("#Delta R(RECO - RECO)");
  h1D->Draw();

  std::cout<<sysType<<" integral total: "<<h1D->Integral(0,h1D->GetNbinsX()+1)<<std::endl;
  std::cout<<sysType<<" integral 0 - 1.0: "<<h1D->Integral(0,h1D->FindBin(1.0))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.5: "<<h1D->Integral(0,h1D->FindBin(0.5))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.05: "<<h1D->Integral(0,h1D->FindBin(0.05))<<std::endl;
  std::cout<<sysType<<" integral 0 - 0.02: "<<h1D->Integral(0,h1D->FindBin(0.02))<<std::endl;

 c->Print(TString::Format("fig_eps/Var%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
 c->Print(TString::Format("fig_png/Var%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotEffVsVar(const std::string & sysType,
		const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   460,500);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TString hName("");
  const int *ptCuts = ptCutsOMTF;
  if(sysType=="Gmt") ptCuts = ptCutsGmt;

  for (int icut=0; icut <=2;++icut){
    if(icut==1) continue;
    float ptCut = OMTFHistograms::ptBins[ptCuts[icut]];
    hName = "h2D"+sysType+varName+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->SetStats(kFALSE);
    hEff->SetMinimum(0.0);
    hEff->SetMaximum(1.04);
    hEff->SetMarkerStyle(21+icut);
    hEff->SetMarkerColor(color[icut]);
    hEff->SetXTitle(varName.c_str());
    hEff->SetYTitle("Efficiency");

    if(sysType=="OMTFDi" && hName.Contains("DeltaPt")){
      hEff->SetMinimum(1E-2);
      c->SetLogy();
    }

    if (icut==0)hEff->DrawCopy("E0");
    else hEff->DrawCopy("same E0");
    std::string nameCut = std::to_string((int)OMTFHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    if(sysType=="OMTFDi" && icut>0) nameCut = "opposite sign";
    if(sysType=="OMTFDi" && icut==0) nameCut = "same sign";

    l.AddEntry(hEff,nameCut.c_str());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotEffVsEta(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("EffVsEta_%s",sysType.c_str()),
			   TString::Format("EffVsEta_%s",sysType.c_str()),
			   800,500);
  
  TLegend l(0.6915995,0.5930233,0.7422325,0.8972868,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.35);

  int iCut = 19;
  std::string hName = "";
  for (int iType=0; iType<3;++iType){
    float ptCut = OMTFHistograms::ptBins[iCut];
    hName = "h2D"+sysType+"Type" + std::to_string(iType) + "EtaVx"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    if(iType==2) hNum->Scale(50.0);
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    hEff->SetStats(kFALSE);
    hEff->SetMinimum(0.0001);
    hEff->SetMaximum(1.04);
    hEff->SetMarkerStyle(21+iType);
    hEff->SetMarkerColor(color[iType]);
    hEff->GetYaxis()->SetTitleOffset(1.0);
    hEff->SetXTitle("generated muon #eta");
    hEff->SetYTitle("Efficiency");
    if (iType==0)hEff->DrawCopy();
    else hEff->DrawCopy("same");
    std::string nameCut = std::to_string((int)OMTFHistograms::ptBins[iCut])+" GeV/c";
    if (iType==0) nameCut = "p_{T}^{#mu}>p_{T}^{cut} + 20 GeV/c";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 5 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c (#epsilon #times 50)";
    l.AddEntry(hEff,nameCut.c_str());
  }
  ///OMTF eta range used for generating patterns.
  TLine *aLine = new TLine(0,0,0,0);
  aLine->SetLineWidth(2);
  aLine->SetLineColor(2);
  aLine->DrawLine(0.83,0,0.83,1.0);
  aLine->DrawLine(-0.83,0,-0.83,1.0);
  aLine->DrawLine(1.24,0,1.24,1.0);
  aLine->DrawLine(-1.24,0,-1.24,1.0);

  l.SetHeader(TString::Format("p_{T} = %d  GeV/c",(int)OMTFHistograms::ptBins[iCut]).Data());
  l.DrawClone();
  c->Print(TString::Format("fig_png/EffVsEta_%s.png",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVsEta_%s.C",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotOMTFVsOther(int iPtCut,
				     const std::string sysType){

  float ptCut = ptBins[iPtCut];

  TCanvas* c = new TCanvas(TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   460,500);

  TLegend l(0.1995614,0.7139831,0.4385965,0.8665254,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(0,1);

  std::string hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
  TH2F* h2D = get2DHistogram(hName);
  if(!h2D) return;
  float  effOther = getEfficiency(h2D,ptCut);
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffOther =DivideErr(hNum,hDenom,"hEffOther","B");

  float effOMTFMatch = 0.0;
  int match = -1;
  float delta = 999.0;
  TH1D *hEffOMTF = 0;
  for (int iCut=0; iCut<=0;++iCut){
    std::string hName = "h2DOMTFPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    TH2F* h2D = get2DHistogram(hName);
    float effOMTF = getEfficiency(h2D,ptCut);
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);    
    hDenom->Add(hNum);

    TH1D* hEffOMTFTmp =DivideErr(hNum,hDenom,"hEffOMTFTmp","B");

    std::cout<<ptCut<<" "<<ptBins[iPtCut+iCut]<<" OMTF: "<<effOMTF
	     <<" "<<sysType<<effOther<<" diff: "<<fabs(effOther-effOMTF)<<std::endl;
    if(fabs(effOther-effOMTF)<delta){
      match = iPtCut+iCut;
      delta = fabs(effOther-effOMTF);
      effOMTFMatch = effOMTF;
      hEffOMTF = hEffOMTFTmp;
    }
  }
  std::cout<<"Other eff: "<<effOther<<std::endl;
  std::cout<<"OMTF eff: "<<effOMTFMatch<<std::endl;

  hEffOther->SetMarkerStyle(23);
  hEffOther->SetMarkerColor(2);

  hEffOMTF->SetXTitle("gen. muon p_{T} [GeV/c]");
  hEffOMTF->SetYTitle("Efficiency");
  hEffOMTF->SetMaximum(1.04);
  hEffOMTF->GetXaxis()->SetRange(2,100);
  hEffOMTF->SetMarkerStyle(8);
  hEffOMTF->SetMarkerColor(1);
  hEffOMTF->SetStats(kFALSE);
  hEffOMTF->DrawCopy();
  hEffOther->DrawCopy("same");

  std::string nameCut(TString::Format("%1.0f GeV/c",ptBins[match]).Data());
  l.AddEntry(hEffOMTF,("OMTF, "+nameCut).c_str());
  std::string tmp = sysType+", %1.0f GeV/c";
  if( int(ptBins[match]*10)%10==5)  tmp = sysType+", %1.1f GeV/c";
  l.AddEntry(hEffOther,TString::Format(tmp.c_str(),ptCut).Data());
  l.DrawClone();

  TLine aLine(0,0,0,0);
  aLine.SetLineColor(2);
  aLine.SetLineWidth(3);
  aLine.DrawLine(ptCut,0,ptCut,1.04);

  c->Print(TString::Format("fig_eps/OMTFVs%s_%d.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OMTFVs%s_%d.png",sysType.c_str(),(int)ptCut).Data());

  c->SetLogy();
  c->Print(TString::Format("fig_eps/OMTFVs%s_%d_log.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OMTFVs%s_%d_log.png",sysType.c_str(),(int)ptCut).Data());

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH2F* OMTFHistograms::makeRateWeights(TH2 *hOrig){

  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
  fIntVxMuRate->SetParameters(-0.235801, -2.82346, 17.162);

  TH2F *hWeights = (TH2F*)hOrig->Clone("hWeights");
  hWeights->Reset();
  TH1D *hPtGen = hOrig->ProjectionX("hPtGen");

  int nEvInBin;
  float ptLow, ptHigh, weight;

  for (int iBin = 1; iBin <= hPtGen->GetNbinsX(); ++iBin){
    ptLow = hPtGen->GetXaxis()->GetBinLowEdge(iBin);
    ptHigh = hPtGen->GetXaxis()->GetBinUpEdge(iBin);
    nEvInBin = hPtGen->GetBinContent(iBin);
    if(nEvInBin<1) nEvInBin = 1;
    weight = (fIntVxMuRate->Eval(ptLow) - fIntVxMuRate->Eval(ptHigh))/nEvInBin;
    for (int iBinY = 1; iBinY<=hOrig->GetNbinsY();++iBinY) hWeights->SetBinContent(iBin,iBinY,weight);
  }
  return hWeights;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1* OMTFHistograms::getRateHisto(std::string sysType,
				 std::string type){

  std::string hName = "h2D"+sysType+"Rate"+type;

  TH2F* h2D_original = (TH2F*)this->get2DHistogram(hName);
  if(! h2D_original) return 0;
  TH2F* h2D = (TH2F*)h2D_original->Clone("h2D");
  if(!h2D) return 0;
  
  //TH2F *hWeights = makeRateWeights(h2D);
  //h2D->Multiply(hWeights);

  TH1D *hRate = h2D->ProjectionY(("hRate"+sysType).c_str());
  if(sysType=="Vx") hRate = h2D->ProjectionX("hRate");

  hRate->SetYTitle("Arbitrary units");
  hRate->SetLineWidth(3);

  if(type=="VsEta") return (TH1*)hRate->Clone("hRateClone");
  if(type=="VsPt") return (TH1*)hRate->Clone("hRateClone");
  if(type=="VsQuality") return (TH1*)hRate->Clone("hRateClone");

  return Integrate(hRate);
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotRate(std::string type){

  TH1 *hRatekBMTF = getRateHisto("kBMTF",type);
  TH1 *hRateVx = getRateHisto("Vx",type);
  TH1 *hRateOMTF = getRateHisto("OMTF",type);

  if(!hRateVx || !hRatekBMTF || !hRateOMTF) return;

  hRateVx->SetLineWidth(3);
  hRatekBMTF->SetLineWidth(3);
  hRateOMTF->SetLineWidth(3);

  hRateVx->SetLineColor(1);
  hRatekBMTF->SetLineColor(2);
  hRateOMTF->SetLineColor(4);

  hRatekBMTF->SetLineStyle(2);

  TCanvas* c = new TCanvas("cRate","Rate",1.5*420,1.5*500);
  c->SetLogy(1);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.60,0.75,0.85,0.85,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type=="Tot"){
    hRateVx->SetAxisRange(2,50);
    hRatekBMTF->SetAxisRange(2,50);
    hRateOMTF->SetAxisRange(2,50);
    hRateVx->SetMinimum(1E1);
    hRateVx->SetMaximum(2E5);
    hRateVx->SetXTitle("p_{T}^{cut} [GeV/c]");

    c->Divide(2);
    TPad *pad1 = (TPad*)c->GetPad(1);
    TPad *pad2 = (TPad*)c->GetPad(2);
    pad1->SetPad(0.01,0.29,0.99,0.99);
    pad2->SetPad(0.01,0.01,0.99,0.29);
    
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    hRateVx->Draw();
    hRatekBMTF->DrawCopy("same");
    hRateOMTF->DrawCopy("same");

    c->cd();
    pad2->Draw();
    pad2->cd();
    
    hRatekBMTF->SetYTitle("kBMTF/OMTF");
    hRatekBMTF->GetXaxis()->SetLabelSize(0.09);
    hRatekBMTF->GetYaxis()->SetLabelSize(0.09);
    hRatekBMTF->GetYaxis()->SetTitleSize(0.09);
    hRatekBMTF->GetYaxis()->SetTitleOffset(0.5);
    hRatekBMTF->Divide(hRateOMTF);
    hRatekBMTF->DrawCopy();
    TLine *aLine = new TLine(0,0,0,0);
    aLine->SetLineWidth(2);
    aLine->SetLineColor(2);
    aLine->DrawLine(2,1.0, 50,1.0);
    c->cd();
  }
  if(type=="VsEta"){
    c->SetLogy(0);
    hRatekBMTF->SetXTitle("muon #eta");
    double max = hRatekBMTF->GetMaximum();
    if(hRateOMTF->GetMaximum()>max) max = hRateOMTF->GetMaximum();
    hRatekBMTF->SetMaximum(1.5*max);
    hRatekBMTF->Draw();
    hRateOMTF->Draw("same");
  }
  if(type=="VsPt"){
    c->SetLogy(1);
    hRatekBMTF->SetXTitle("p_{T}^{gen} [GeV/c]");
    hRatekBMTF->SetAxisRange(2,100);
    double max = hRatekBMTF->GetMaximum();
    if(hRateOMTF->GetMaximum()>max) max = hRateOMTF->GetMaximum();
    hRatekBMTF->SetMaximum(10*max);
    hRatekBMTF->Draw();
    hRateOMTF->Draw("same");
  }
  if(type=="VsQuality"){
    c->SetLogy(0);
    c->SetLeftMargin(0.2);
    c->SetRightMargin(0.3);
    
    TH2F *hEff20 = get2DHistogram("h2DOMTFQuality20");
    if(!hEff20) return;
    TH1F *hRateValues = sortRateHisto((TH1F*)hRateOMTF, hEff20, "rate");
    TH1F *hEffValues = sortRateHisto((TH1F*)hRateOMTF, hEff20, "eff");
    hRateValues->SetAxisRange(0,10);
    hRateValues->GetYaxis()->SetTitleOffset(1.7);
    hRateValues->GetYaxis()->SetTitle("Rate [arbitrary units]");
    hRateValues->Draw();
    c->Update();
    
    Float_t rightmax = 0.2;//1.1*hEffValues->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hEffValues->SetLineColor(kRed);
    hEffValues->Scale(scale);
    hEffValues->Draw("same hist");
    //draw an axis on the right side
    std::cout<<"gPad->GetUxmax(): "<<gPad->GetUxmax()<<std::endl;
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			      gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetTitle("Efficiency for p_{T}^{cut}=20 GeV/c @ p_{T}^{gen}>40 GeV/c");
    axis->SetTitleOffset(1.7);
    axis->SetTitleColor(kRed);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();

    c->Print(("fig_eps/Rate"+type+".eps").c_str());
    c->Print(("fig_png/Rate"+type+".png").c_str());
    return;
  }
 
  leg->AddEntry(hRatekBMTF,"kBMTF");
  leg->AddEntry(hRateOMTF,"OMTF");
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotEffVsRate(int iPtCut){

  TCanvas* c = new TCanvas("cEffVsRate","EffVsRate",460,500);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.2149123,0.6716102,0.4539474,0.8644068,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  TH1 *hRatekBMTF = getRateHisto("kBMTF","Tot");
  TH1 *hRateOMTF = getRateHisto("OMTF","Tot");
  if(!hRatekBMTF || !hRateOMTF) return;

  TGraph *grkBMTF = new TGraph();
  TGraph *grOMTF = new TGraph();

  grkBMTF->SetMarkerColor(2);
  grkBMTF->SetMarkerStyle(20);
  grkBMTF->SetMarkerSize(1.3);

  grOMTF->SetMarkerColor(4);
  grOMTF->SetMarkerStyle(22);
  grOMTF->SetMarkerSize(1.3);

  float effNom=0, rateNom=0;
  for(int iCut=-3;iCut<=2;++iCut){
    ///
    std::string hName = "h2DOMTFPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    float effOMTF = getEfficiency(h2D,ptBins[iPtCut]);
    ///
    hName = "h2DkBMTFPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    float effkBMTF = getEfficiency(h2D,ptBins[iPtCut]);
    ///
    int iBin = hRatekBMTF->FindBin(ptBins[iPtCut+iCut]);
    float ratekBMTF = hRatekBMTF->GetBinContent(iBin);
    grkBMTF->SetPoint(iCut+2,effkBMTF,ratekBMTF);
    if(iCut==0){
      effNom = effkBMTF;
      rateNom = ratekBMTF;
    }
    ///
    float rateOMTF = hRateOMTF->GetBinContent(iBin);
    grOMTF->SetPoint(iCut+2,effOMTF,rateOMTF);
    std::cout<<ptBins[iPtCut+iCut]<<" effOMTF: "<<effOMTF<<" rateOMTF: "<<rateOMTF<<" "
	     <<ptBins[iPtCut+iCut]<<" effkBMTF: "<<effkBMTF<<" ratekBMTF: "<<ratekBMTF<<std::endl;
  }
  
  Double_t maxY, minY, tmp;
  float effMin=0.8, effMax=1.0;
  grkBMTF->GetPoint(0,tmp,maxY);
  effMax = tmp;
  grOMTF->GetPoint(4,tmp,minY);
  maxY*=1.3;
  minY*=0.7;
  effMin = tmp;
  effMin-=0.05;
  effMax+=0.05;

  std::cout<<"effMin: "<<effMin<<" effMax: "<<effMax<<std::endl;

  TH1F *hFrame = new TH1F("hFrame",";Efficiency;Rate",2,effMin,effMax);
  hFrame->SetMinimum(minY);
  hFrame->SetMaximum(maxY);
  hFrame->GetYaxis()->SetTitleOffset(1.7);
  hFrame->GetYaxis()->SetLabelSize(0.04);
  hFrame->GetXaxis()->SetLabelSize(0.04);
  hFrame->SetXTitle(TString::Format("Efficiency for p_{T}^{gen}>%d GeV/c",(int)ptBins[iPtCut]));
  hFrame->Draw();

  grkBMTF->Draw("P");
  grOMTF->Draw("P");

  TMarker *aMarker = new TMarker(effNom,rateNom,4);
  aMarker->SetMarkerSize(2.1);
  aMarker->Draw();

  leg->AddEntry(aMarker,"Nominal GMT cut","p");
  leg->AddEntry(grkBMTF,"GMT","p");
  leg->AddEntry(grOMTF,"OMTF","p");
  leg->Draw();

  c->Print(TString::Format("fig_eps/RateVsEff_%d.eps",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.png",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.C",(int)ptBins[iPtCut]).Data());

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotGhostHistos(const std::string & sysType,
				    const std::string & type){

  TCanvas* c = new TCanvas("cGhostHistos","GhostHistos",460,500);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.2149123,0.6716102,0.4539474,0.8644068,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  std::string hName = "h2DGhostsVsProcessor"+sysType+type;
  TH2F *h2DGhostsVsProcessor = (TH2F*)this->get2DHistogram(hName);

  h2DGhostsVsProcessor->GetYaxis()->SetTitleOffset(1.7);
  h2DGhostsVsProcessor->SetXTitle("Processor number");
  h2DGhostsVsProcessor->SetYTitle("Number of candidates");
  h2DGhostsVsProcessor->SetMarkerSize(2.8);
  h2DGhostsVsProcessor->Draw("text");
  h2DGhostsVsProcessor->GetYaxis()->SetRange(1,5);

  c->Print(("fig_eps/GhostsVsProcessor"+sysType+type+".eps").c_str());
  c->Print(("fig_png/GhostsVsProcessor"+sysType+type+".png").c_str());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
float  OMTFHistograms::getEfficiency(TH2F *h2D, float ptCut){

  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  //Mean eff above pt cut
  int binLow = hEffTmp->FindBin(ptCut);
  int binHigh = hEffTmp->FindBin(ptCut+5);
  float range = hEffTmp->GetBinLowEdge(binHigh+1) - hEffTmp->GetBinLowEdge(binLow);
  float eff = hEffTmp->Integral(binLow,binHigh,"width")/range;

  return eff;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotSingleHistogram(std::string hName){

  TH2F* h2D = get2DHistogram(hName);
  TH1F* h1D = get1DHistogram(hName);
  if(!h2D && !h1D) return;
	
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(h2D) {
    std::cout<<"myDirCopy: "<<myDirCopy<<std::endl;
    h2D->SetDirectory(myDirCopy);
    h2D->SetLineWidth(3);
    h2D->Scale(1.0/h2D->Integral(0,h2D->GetNbinsX()+1));
    h2D->SetXTitle("p_{T}^{GEN}");
    h2D->SetYTitle("p_{T}^{REC}");
    h2D->GetYaxis()->SetTitleOffset(1.4);
    h2D->SetStats(kFALSE);
    h2D->Draw("candle2");
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
  if(h1D) {
    h1D->SetDirectory(myDirCopy);
    h1D->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));    
    h1D->GetXaxis()->SetRange(1,h1D->GetNbinsX()+1);
    h1D->SetXTitle("X");
    h1D->SetYTitle("Y");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw("");
    h1D->Print();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::plotLLH(){

   TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);

   TH1F *h1 = get1DHistogram("h1DLLH_Low");
   TH1F *h2 = get1DHistogram("h1DLLH_High");

   if(!h1 || !h2) return;

   TLegend l(0.18,0.75,0.35,0.87,NULL,"brNDC");
   l.SetTextSize(0.05);
   l.SetFillStyle(4000);
   l.SetBorderSize(0);
   l.SetFillColor(10);

   h1->SetLineWidth(3);
   h2->SetLineWidth(3);

   h1->SetLineColor(2);
   h2->SetLineColor(1);
   
   h1->Scale(1.0/h1->Integral(0,h1->GetNbinsX()+1));
   h2->Scale(1.0/h2->Integral(0,h2->GetNbinsX()+1));

   h1->GetXaxis()->SetRange(1,h1->GetNbinsX()+1);
   h2->GetXaxis()->SetRange(1,h2->GetNbinsX()+1);
   
   h1->SetXTitle("Events");
   h1->SetYTitle("Candidate LLH");
   h1->GetYaxis()->SetTitleOffset(1.4);
   h1->SetStats(kFALSE);
   h2->Draw("");
   h1->Draw("same");

   l.AddEntry(h1,"p_{T}^{GEN}<4 AND p_{T}^{REC}>20 GeV/c");
   l.AddEntry(h2,"p_{T}^{GEN}>20 AND p_{T}^{REC}>20 GeV/c");
   l.Draw();
   
   c->Print("fig_png/LLH_values.png");
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1F * OMTFHistograms::sortRateHisto(TH1F *h1DRate, TH2F *h2DEff, std::string by){

  if(!h1DRate || !h2DEff) return 0;

  TH1F *hRateSorted = (TH1F*)h1DRate->Clone("hRateSorted");    
  hRateSorted->Clear();
  ///Calculate efficiency
  TH1D *hNum = h2DEff->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2DEff->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  hNum->Scale(1.0/hDenom->Integral());
  TH1D* hEff = hNum;
  ////
  std::map<float,std::string> rateMap;
  std::map<float,int> rateMapBin;
  ///Fill map with key = rate, value = bin label. The map will be sorted automatically by the rate value.
  ///Second map has bin number as value. Needed to extract efficiency values.
  TRandom3 aRndm;

  for(int iBin=1;iBin<h1DRate->GetXaxis()->GetNbins();++iBin){
    float rate = h1DRate->GetBinContent(iBin);
    rate+=0.001*aRndm.Uniform();///Randomize rate values, to avoid having two same keys.
    rateMap[rate] = h1DRate->GetXaxis()->GetBinLabel(iBin);
    rateMapBin[rate] = iBin;
  }
  ///Fill histo copy with sorted values
  unsigned int iBin = 1;
  for (auto it = rateMap.rbegin(); it!= rateMap.rend(); ++it){
    if(iBin<20){
      std::cout<<"iBin: "<<iBin
	       <<" Quality: "<<it->second
	       <<" rate: "<<it->first
	       <<" efficiency: "<<hEff->GetBinContent(rateMapBin[it->first])
	       <<std::endl;
    }
    if(by=="rate") hRateSorted->SetBinContent(iBin, it->first);
    if(by=="eff")  hRateSorted->SetBinContent(iBin, hEff->GetBinContent(rateMapBin[it->first]));
    hRateSorted->GetXaxis()->SetBinLabel(iBin, it->second.c_str());
    ++iBin;
  }
  
  return hRateSorted;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
