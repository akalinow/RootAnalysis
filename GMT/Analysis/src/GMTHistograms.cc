#include <iostream>
#include <cmath>

#include "GMTHistograms.h"
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
const float GMTHistograms::ptBins[36]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
		 10., 12., 14., 16., 18., 20., 22,  24, 26, 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
		 160., 200};

const int GMTHistograms::color[6] = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
//Single mu
const int GMTHistograms::ptCutsGmt[4] =     {0, 14, 19, 20};
const int GMTHistograms::ptCutsOMTF[4] =     {0, 14, 19, 20};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::GMTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
GMTHistograms::~GMTHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string GMTHistograms::getTemplateName(const std::string& name){

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
  if(name.find("3DBending")!=std::string::npos) templateName = "h3DBendingTemplate";
  
  if(name.find("LLH")!=std::string::npos) templateName = "h1DLLHTemplate";
  if(name.find("HitsPattern")!=std::string::npos) templateName = "h1DHitsPatternTemplate";
  if(name.find("h2DDeltaVsDelta")!=std::string::npos) templateName = "h2DDeltaVsDeltaTemplate";
  
  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void GMTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

 //Make template histos
 add1DHistogram("h1DDeltaEtaTemplate","",11,-0.83,0.83,file_);

 add1DHistogram("h1DDeltaPhiTemplate","",5*32,-M_PI,M_PI,file_);

 ///Efficiency histos
 add2DHistogram("h2DPtTemplate","",150,0,150,2,-0.5,1.5,file_);
 add2DHistogram("h2DHighPtTemplate","",50,50,550,2,-0.5,1.5,file_);

 add2DHistogram("h2DPtVsPtTemplate","",404,0,202,404,0,202,file_);

 add2DHistogram("h2DEtaHitTemplate","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHitTemplate","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaVxTemplate","",10,0.8,1.3,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiVxTemplate","",4*32,-3.2,3.2,2,-0.5,1.5,file_);

 add2DHistogram("h2DQualityTemplate","",201,-0.5,200.5,2,-0.5,1.5,file_);

 //Rate histos
 add2DHistogram("h2DRateTotTemplate","",404,1,202,404,1,202,file_);
 add2DHistogram("h2DRateVsEtaTemplate","",404,1,202,10,0.8,1.3,file_);

 add2DHistogram("h2DDeltaPhiTemplate","",30,-1,1,2,-0.5,1.5,file_);
 add2DHistogram("h2DDeltaPtTemplate","",21,-0.5,20.5,2,-0.5,1.5,file_);

 add2DHistogram("h2DRateVsPtTemplate","",404,1,202,100,0,50,file_);
 
 add2DHistogram("h2DRateVsQualityTemplate","",404,1,202,201,-0.5,200.5,file_);
 add2DHistogram("h2DGhostsVsProcessorTemplate","",6,-0.5,5.5,5,-0.5,4.5,file_);

 //Likelihood histos
 add1DHistogram("h1DLLHTemplate","",40,0,20,file_);
 add1DHistogram("h1DHitsPatternTemplate","",101,-0.5,100.5,file_);

 //Pattern histos
 add3DHistogram("h3DBendingTemplate","",50, 0, 100, 296,-250.5,45.5, 18, -0.5, 17.5, file_);

 
 add2DHistogram("h2DDeltaVsDeltaTemplate","",176,-150.5,25.5, 176,-150.5,25.5, file_);
 
 histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void GMTHistograms::finalizeHistograms(){

  AnalysisHistograms::finalizeHistograms();
  utilsL1RpcStyle()->cd();

  finaliseGoldenPatterns("h3DBending");
  finaliseGoldenPatterns("h3DBendingRotated");

    for(int iPtCode=1;iPtCode<=30;++iPtCode){
    plotOMTFVsOther(iPtCode,"uGMT_emu");
  }

  std::vector<int> ptCuts = {10, 13, 15, 16, 18, 19, 20, 21, 22, 23};
  for(auto iCut: ptCuts){
    plotEffVsRate(iCut);
  }

  plotRate("VsEta_quality");
  plotEffVsEtaVsQuality();

  plotRate("VsEta");
  plotRate("VsPt");
  plotRate("VsQuality");
  plotRate("VsEta_quality");
    
  bool doHigh = true;
  plotEffPanel("uGMT_emu", doHigh);
  plotEffPanel("uGMT_emu");
  plotEffVsEta("uGMT_emu");

  plotEffVsVar("uGMT_emu","EtaVx");
  plotEffVsVar("uGMT_emu","PhiVx");

  plotSingleHistogram("h2DuGMT_emuPtRecVsPtGen");

  for(int iPtCode=1;iPtCode<=30;++iPtCode){
      plotOMTFVsOther(iPtCode,"uGMT_emu");
  }
   
  plotSingleHistogram("h1DLLH_Low");
  plotSingleHistogram("h1DLLH_High");
  plotLLH();

  ////////////////////////////// For GMT muons, those histograms are not created 
  plotSingleHistogram("h1DHitsPattern_Low_HitCount");  
  plotSingleHistogram("h1DHitsPattern_High_HitCount");

  plotSingleHistogram("h1DHitsPattern_Low_RefLayer");
  plotSingleHistogram("h1DHitsPattern_High_RefLayer");

  plotSingleHistogram("h1DHitsPattern_Low_RefPhi");
  plotSingleHistogram("h1DHitsPattern_High_RefPhi");

  plotSingleHistogram("h1DDeltaEta_Low");
  plotSingleHistogram("h1DDeltaEta_High");
  //////////////////////////////////////////////////

  plotQuantiles("h2DuGMT_emuPtRecVsPtGen");
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* GMTHistograms::Integrate(TH1 * histoD) {

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
  for (Int_t i = histoD->GetNbinsX()+1; i > 0; i--) {
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
TH1D * GMTHistograms::DivideErr(TH1D * h1,
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
void GMTHistograms::plotEffPanel(const std::string & sysType, bool doHigh){

  TCanvas* c = new TCanvas(TString::Format("EffVsPt_%s",sysType.c_str()),
			   TString::Format("EffVsPt_%s",sysType.c_str()),
			   460,500);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TString hName("");
  const int *ptCuts = ptCutsOMTF;
 
  for (int icut=0; icut <=3;++icut){
    float ptCut = GMTHistograms::ptBins[ptCuts[icut]];
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
    TString nameCut = TString::Format("%d", (int)GMTHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
  }
  l.DrawClone();
  if(!doHigh) c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
  else  c->Print(TString::Format("fig_png/PanelVsHighPt_%s.png",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotVar(const std::string & sysType,
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
void GMTHistograms::plotEffVsVar(const std::string & sysType,
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
 
  for (int icut=0; icut<2;++icut){
    float ptCut = GMTHistograms::ptBins[ptCuts[icut]];
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

    if (icut==0)hEff->DrawCopy("E0");
    else hEff->DrawCopy("same E0");
    std::string nameCut = std::to_string((int)ptCut)+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.c_str());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotEffVsEta(const std::string & sysType){

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

  int iCut = 18;
  std::string hName = "";
  for (int iType=0; iType<3;++iType){
    float ptCut = GMTHistograms::ptBins[iCut];
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
    std::string nameCut = std::to_string((int)GMTHistograms::ptBins[iCut])+" GeV/c";
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

  l.SetHeader(TString::Format("p_{T} = %d  GeV/c",(int)GMTHistograms::ptBins[iCut]).Data());
  l.DrawClone();
  c->Print(TString::Format("fig_png/EffVsEta_%s.png",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVsEta_%s.C",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotOMTFVsOther(int iPtCut,
				     const std::string sysType){

  float ptCut = ptBins[iPtCut];

  TCanvas* c = new TCanvas(TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   TString::Format("OMTFVsOther_%d",(int)ptCut).Data(),
			   460,500);

  TLegend l(0.2,0.65,0.44,0.86,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(0,1);

  std::string hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
  TH2F* h2D = get2DHistogram(hName);
  if(!h2D) return;
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffOther =DivideErr(hNum,hDenom,"hEffOther","B");
  hEffOther->SetMarkerStyle(23);
  hEffOther->SetMarkerColor(2);

  hName = "h2DOMTFPt"+std::to_string((int)ptCut);
  h2D = get2DHistogram(hName);
  hNum = h2D->ProjectionX("hNum",2,2);
  hDenom = h2D->ProjectionX("hDenom",1,1);    
  hDenom->Add(hNum);

  TH1D* hEffOMTF =DivideErr(hNum,hDenom,"hEffOMTFTmp","B");
  hEffOMTF->SetXTitle("gen. muon p_{T} [GeV/c]");
  hEffOMTF->SetYTitle("Efficiency");
  hEffOMTF->SetMaximum(1.04);
  hEffOMTF->GetXaxis()->SetRange(2,100);
  hEffOMTF->SetMarkerStyle(8);
  hEffOMTF->SetMarkerColor(1);
  hEffOMTF->SetStats(kFALSE);
  hEffOMTF->DrawCopy();
  hEffOther->DrawCopy("same");

  std::string tmp = "p_{T} #geq ";
  if(int(ptCut*10)%10==5) tmp += "%1.1f GeV/c";
  else   tmp += "%1.0f GeV/c";
  l.AddEntry((TObject*)0, TString::Format(tmp.c_str(),ptCut).Data(), "");
  l.AddEntry((TObject*)0, "", "");
  tmp = sysType;
  if(sysType=="BMTF") tmp = "LUT NN";
  if(sysType=="EMTF") tmp = "TF NN";
  l.AddEntry(hEffOther, tmp.c_str());
  l.AddEntry((TObject*)0, "", "");
  l.AddEntry(hEffOMTF, "Phase1");
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
TH2F* GMTHistograms::makeRateWeights(TH2 *hOrig){

  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","0.1*TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
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
    for (int iBinY = 0; iBinY<=hOrig->GetNbinsY()+1;++iBinY) hWeights->SetBinContent(iBin,iBinY,weight);
  }
  return hWeights;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1* GMTHistograms::getRateHisto(std::string sysType,
				 std::string type){

  std::string hName = "h2D"+sysType+"Rate"+type;

  TH2F* h2D_original = (TH2F*)this->get2DHistogram(hName);
  if(! h2D_original) return 0;
  TH2F* h2D = (TH2F*)h2D_original->Clone("h2D");
  if(!h2D) return 0;

  if(selectionFlavours_.size() &&
     selectionFlavours_[0].find("NU_RATE")==std::string::npos){
    TH2F *hWeights = makeRateWeights(h2D);
    h2D->Multiply(hWeights);
  }
  
  TH1D *hRate = h2D->ProjectionY(("hRate"+sysType).c_str());
  if(sysType=="Vx") hRate = h2D->ProjectionX("hRate");

  hRate->SetYTitle("Arbitrary units");
  hRate->SetLineWidth(3);

  if(type.find("VsEta")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");
  if(type.find("VsPt")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");
  if(type.find("VsQuality")!=std::string::npos) return (TH1*)hRate->Clone("hRateClone");

  return Integrate(hRate);
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotRate(std::string type){
  
  TH1 *hRateuGMT_emu = getRateHisto("uGMT_emu",type);
  TH1 *hRateVx = getRateHisto("Vx",type);
  TH1 *hRateOMTF = getRateHisto("OMTF",type);
  TH1 *hRateEMTF = getRateHisto("EMTF",type);

  if(type.find("VsEta_quality")!=std::string::npos){
   hRateuGMT_emu = getRateHisto("uGMT_emu","VsEta_quality0");
   hRateVx = getRateHisto("uGMT_emu","VsEta_quality1");
   hRateOMTF = getRateHisto("OMTF","VsEta_quality2");
   hRateEMTF = getRateHisto("uGMT_emu","VsEta_quality3");    
  }

  if(!hRateVx || !hRateuGMT_emu || !hRateOMTF || !hRateEMTF) return;

  hRateVx->SetLineWidth(3);
  hRateuGMT_emu->SetLineWidth(3);
  hRateEMTF->SetLineWidth(3);
  hRateOMTF->SetLineWidth(3);

  hRateVx->SetLineColor(1);
  hRateuGMT_emu->SetLineColor(2);
  hRateEMTF->SetLineColor(kGreen-2);
  hRateOMTF->SetLineColor(4);

  hRateuGMT_emu->SetLineStyle(2);
  hRateEMTF->SetLineStyle(1);

  TCanvas* c = new TCanvas("cRate","Rate",1.5*420,1.5*500);
  c->SetLogy(1);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.60,0.75,0.85,0.9,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type.find("Tot")!=std::string::npos){
    hRateVx->SetAxisRange(2,50);
    hRateuGMT_emu->SetAxisRange(2,50);
    hRateEMTF->SetAxisRange(2,50);
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
    pad1->SetGrid(1,0);
    hRateVx->Draw();
    hRateuGMT_emu->DrawCopy("same");
    hRateEMTF->DrawCopy("same");
    hRateOMTF->DrawCopy("same");

    std::cout<<"Rate OMTF @ 20 GeV: "<< hRateOMTF->GetBinContent(hRateOMTF->FindBin(20-0.01))<<std::endl;
    std::cout<<"Rate other @ 20 GeV: "<< hRateuGMT_emu->GetBinContent(hRateuGMT_emu->FindBin(20-0.01))<<std::endl;

    c->cd();
    pad2->Draw();
    pad2->cd();
    
    hRateuGMT_emu->SetYTitle("new model/Phase1");
    hRateuGMT_emu->GetXaxis()->SetLabelSize(0.09);
    hRateuGMT_emu->GetYaxis()->SetLabelSize(0.09);
    hRateuGMT_emu->GetYaxis()->SetTitleSize(0.09);
    hRateuGMT_emu->GetYaxis()->SetTitleOffset(0.5);
    hRateuGMT_emu->Divide(hRateOMTF);
    hRateuGMT_emu->SetMaximum(1.5);
    hRateuGMT_emu->SetMinimum(0.2);
    hRateuGMT_emu->DrawCopy();

    hRateEMTF->Divide(hRateOMTF);
    hRateEMTF->DrawCopy("same");
    TLine *aLine = new TLine(0,0,0,0);
    aLine->SetLineWidth(2);
    aLine->SetLineColor(2);
    aLine->DrawLine(5,1.0, 50,1.0);
    c->cd();
  }  
  if(type.find("VsEta_quality")!=std::string::npos){
    c->SetLogy(0);
    hRateuGMT_emu->SetXTitle("muon #eta");
    double max = hRateuGMT_emu->GetMaximum();
    if(hRateOMTF->GetMaximum()>max) max = hRateOMTF->GetMaximum();
    hRateuGMT_emu->SetMaximum(1.5*max);
    hRateuGMT_emu->Draw();
    hRateVx->Draw("same");
    hRateOMTF->Draw("same");
    hRateEMTF->Draw("same");
  }
  else if(type.find("VsEta")!=std::string::npos){
    c->SetLogy(0);
    hRateuGMT_emu->SetXTitle("muon #eta");
    double max = hRateuGMT_emu->GetMaximum();
    if(hRateOMTF->GetMaximum()>max) max = hRateOMTF->GetMaximum();
    hRateuGMT_emu->SetMaximum(1.5*max);
    hRateuGMT_emu->Draw();
    hRateOMTF->Draw("same");
    hRateEMTF->Draw("same");
  }
  if(type=="VsPt"){
    c->SetLogy(1);
    hRateuGMT_emu->SetXTitle("p_{T}^{gen} [GeV/c]");
    hRateuGMT_emu->SetAxisRange(2,100);
    double max = hRateuGMT_emu->GetMaximum();
    if(hRateOMTF->GetMaximum()>max) max = hRateOMTF->GetMaximum();
    hRateuGMT_emu->SetMaximum(10*max);
    hRateuGMT_emu->Draw();
    hRateOMTF->Draw("same");
    hRateEMTF->Draw("same");
  }
  if(type=="VsQuality"){
    c->SetLogy(0);
    c->SetLeftMargin(0.2);
    c->SetRightMargin(0.3);
    c->SetBottomMargin(0.15);
    
    TH2F *hEff20 = get2DHistogram("h2DEMTFQuality20");    
    if(!hEff20) return;
    TH1F *hRateValues = sortRateHisto((TH1F*)hRateEMTF, hEff20, "rate");
    TH1F *hEffValues = sortRateHisto((TH1F*)hRateEMTF, hEff20, "eff");
    hRateValues->SetAxisRange(0,10);
    hRateValues->GetYaxis()->SetTitleOffset(1.7);
    hRateValues->GetYaxis()->SetTitle("Rate [arbitrary units]");
    hRateValues->SetLineColor(kBlack);
    hRateValues->Draw();
    c->Update();
    
    Float_t rightmax = hEffValues->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    hEffValues->SetLineColor(kRed);
    hEffValues->Scale(scale);
    hEffValues->Draw("same hist");
    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
			      gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetTitle("Efficiency for p_{T}^{cut}=20 GeV/c @ p_{T}^{gen}>20 GeV/c");
    axis->SetTitleOffset(2.0);
    axis->SetTitleColor(kRed);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();

    c->Print(("fig_eps/Rate"+type+".eps").c_str());
    c->Print(("fig_png/Rate"+type+".png").c_str());
    return;
  }

  if(type=="VsEta_quality"){
    leg->AddEntry(hRateuGMT_emu,"Q=0");
    leg->AddEntry(hRateVx,"Q=1");
    leg->AddEntry(hRateOMTF,"Q=2");
    leg->AddEntry(hRateEMTF,"Q=3");
  }
  else{
    leg->AddEntry(hRateOMTF,"Phase 1");
    leg->AddEntry(hRateuGMT_emu,"LUT NN");
    leg->AddEntry(hRateEMTF,"TF NN");
  }
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotEffVsRate(int iPtCut){

  double ptCut = GMTHistograms::ptBins[iPtCut];
  double xMin = 9999.0, xMax = -999.0;
  double yMin = 1E8, yMax = -999.0;

  TGraph *aGraph = new TGraph();
  aGraph->SetMarkerSize(1.3);
  std::map<std::string, TGraph*> sysGraphs;
  std::vector<std::string >sysTypes = {"OMTF", "uGMT_emu", "EMTF"};
  for(auto && sysType: sysTypes){
    for(int iQuality=0;iQuality<4;++iQuality){
      std::string selType = std::string(TString::Format("quality%d",iQuality));      
      TH1 *hRate = getRateHisto(sysType,"Tot_"+selType);
      if(!hRate) return;
      std::string hName = "h2D"+sysType+selType+"Pt"+std::to_string((int)ptCut);
      TH2F* hEff = this->get2DHistogram(hName);
      if(!hEff) return;
      float efficiency = getEfficiency(hEff, ptCut);
      int iBin = hRate->FindBin(ptCut);
      float rate = hRate->GetBinContent(iBin);
      aGraph->SetPoint(iQuality, efficiency, rate);
      if(efficiency<xMin) xMin = efficiency;
      if(efficiency>xMax) xMax = efficiency;
      if(rate<yMin) yMin = rate;
      if(rate>yMax) yMax = rate;
    }
    std::cout<<sysType<<std::endl;
    aGraph->Print();
    if(sysType=="OMTF"){
      aGraph->SetMarkerColor(1);
      aGraph->SetMarkerStyle(21);
    }
    if(sysType=="EMTF"){
      aGraph->SetMarkerColor(4);
      aGraph->SetMarkerStyle(22);
    }
    if(sysType=="uGMT_emu"){
      aGraph->SetMarkerColor(2);
      aGraph->SetMarkerStyle(20);
    }
    sysGraphs[sysType] = (TGraph*)aGraph->Clone(sysType.c_str());
  }
  
  xMin-=0.04;
  xMax+=0.02;
  yMin*=0.95;
  yMax*=1.05;

  TH1F *hFrame = new TH1F("hFrame",";Efficiency;Rate",2,xMin, xMax);
  hFrame->SetMinimum(yMin);
  hFrame->SetMaximum(yMax);
  hFrame->GetYaxis()->SetTitleOffset(1.7);
  hFrame->GetYaxis()->SetLabelSize(0.04);
  hFrame->GetXaxis()->SetLabelSize(0.04);
  hFrame->GetXaxis()->SetNdivisions(505);
  hFrame->SetXTitle(TString::Format("Efficiency for %d < p_{T}^{gen} < 100",(int)ptBins[iPtCut]));

  TCanvas* c = new TCanvas("cEffVsRate","EffVsRate",460,500);
  c->SetGrid(1,1);

  hFrame->Draw();
  sysGraphs["OMTF"]->Draw("P");
  sysGraphs["uGMT_emu"]->Draw("P");
  sysGraphs["EMTF"]->Draw("P");

  TLegend *leg = new TLegend(0.2149123,0.6716102,0.4539474,0.8644068,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetHeader(TString::Format("p_{T}^{REC} #geq %d GeV/c", (int)ptBins[iPtCut]));
  leg->AddEntry(sysGraphs["OMTF"],"Phase 1","p");
  leg->AddEntry(sysGraphs["uGMT_emu"],"LUT NN","p");
  leg->AddEntry(sysGraphs["EMTF"],"TF NN","p");
  leg->Draw();

  c->Print(TString::Format("fig_eps/RateVsEff_%d.eps",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.png",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.C",(int)ptBins[iPtCut]).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotEffVsEtaVsQuality(){

  std::string hName = "h2DuGMT_emu_quality1EtaVx25";
  TH1D * hEffVx = getEfficiencyHisto(hName);

  hName = "h2DuGMT_emu_quality0EtaVx25";
  TH1D * hEffuGMT_emu = getEfficiencyHisto(hName);

  hName = "h2DuGMT_emu_quality2EtaVx25";
  TH1D * hEffOMTF = getEfficiencyHisto(hName);

  hName = "h2DuGMT_emu_quality3EtaVx25";
  TH1D * hEffEMTF = getEfficiencyHisto(hName);

  if(!hEffVx || !hEffuGMT_emu || !hEffOMTF || !hEffEMTF) return;

  hEffVx->SetLineWidth(3);
  hEffuGMT_emu->SetLineWidth(3);
  hEffEMTF->SetLineWidth(3);
  hEffOMTF->SetLineWidth(3);

  hEffVx->SetLineColor(1);
  hEffuGMT_emu->SetLineColor(2);
  hEffEMTF->SetLineColor(9);
  hEffOMTF->SetLineColor(4);

  hEffuGMT_emu->SetLineStyle(2);
  hEffEMTF->SetLineStyle(3);

  TCanvas* c = new TCanvas("EffVsEtaVsQuality","EffVsEtaVsQuality",460,500);
  c->SetGrid(0,1);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);


  hEffuGMT_emu->SetXTitle("muon #eta");
  hEffuGMT_emu->SetYTitle("Efficiency for p_{T}^{cut}<p_{T}^{Gen}<p_{T}^{cut}+50 GeV");  
  hEffuGMT_emu->SetMaximum(1.0);
  hEffuGMT_emu->Draw();
  hEffVx->Draw("same");
  hEffOMTF->Draw("same");
  hEffEMTF->Draw("same");

  l.AddEntry(hEffuGMT_emu,"Q=0");
  l.AddEntry(hEffVx,"Q=1");
  l.AddEntry(hEffOMTF,"Q=2");
  l.AddEntry(hEffEMTF,"Q=3");  
  l.DrawClone();

  c->Print("fig_eps/EffVsEta.eps");
  c->Print("fig_png/EffVsEta.png");
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotGhostHistos(const std::string & sysType,
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
float  GMTHistograms::getEfficiency(TH2F *h2D, float ptCut){

  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  //Mean eff above pt cut
  int binLow = hEffTmp->FindBin(ptCut);
  int binHigh = hEffTmp->FindBin(100);
  float range = hEffTmp->GetBinLowEdge(binHigh+1) - hEffTmp->GetBinLowEdge(binLow);
  float eff = hEffTmp->Integral(binLow,binHigh,"width")/range;
  return eff;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1D * GMTHistograms::getEfficiencyHisto(const std::string & hName){

  TH2F* h2D = this->get2DHistogram(hName);
  if(!h2D) return 0;
  
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);  
  hDenom->Add(hNum);
  TH1D* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  return hEffTmp;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotSingleHistogram(std::string hName){

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
    h2D->Scale(1.0/h2D->Integral());
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
void GMTHistograms::plotLLH(){

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
   
   h1->SetYTitle("Events");
   h1->SetXTitle("Posterior #sigma");
   h1->GetYaxis()->SetTitleOffset(1.4);
   h1->SetStats(kFALSE);
   double histoMax = std::max(h1->GetMaximum(), h2->GetMaximum());
   h1->SetMaximum(1.4*histoMax);
   h1->Draw("");
   h2->Draw("same");

   l.AddEntry(h1,"p_{T}^{GEN}<10 AND p_{T}^{REC} #geq 20 GeV/c");
   l.AddEntry(h2,"p_{T}^{GEN}>20 AND p_{T}^{REC} #geq 20 GeV/c");
   l.Draw();
   
   c->Print("fig_png/LLH_values.png");
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1F * GMTHistograms::sortRateHisto(TH1F *h1DRate, TH2F *h2DEff, std::string by){

  if(!h1DRate || !h2DEff) return 0;

  TH1F *hRateSorted = (TH1F*)h1DRate->Clone("hRateSorted");    
  hRateSorted->Reset();
  ///Calculate efficiency
  TH1D *hNum = h2DEff->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2DEff->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);  
  hNum->Scale(1.0/hDenom->Integral(0,hDenom->GetNbinsX()+1));
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
void GMTHistograms::finaliseGoldenPatterns(std::string hName){

  TH3F *hAllLayers = get3DHistogram(hName, true);
  if(!hAllLayers) return;
  hAllLayers->SetDirectory(myDirCopy);

  TH2F *hLayer = 0;  
  double threshold = 8E-13;
  std::string hNameTmp = "";
  for(int iLayer=0;iLayer<18;++iLayer){
    hAllLayers->GetZaxis()->SetRangeUser(iLayer, iLayer);
    hLayer = (TH2F*)hAllLayers->Project3D("yx");
    hNameTmp = hName+"_iLayer_"+std::to_string(iLayer);
    hLayer->SetName(hNameTmp.c_str());
    hLayer->SetDirectory(myDirCopy);
    double norm = 1.0;
    double normalisedValue = 0.0;
    int nBinsX = hLayer->GetNbinsX();
    int nBinsY = hLayer->GetNbinsY();
    ///Normalise each pT bin.
    for(int iBinX=0;iBinX<=nBinsX;++iBinX){
      norm = hLayer->Integral(iBinX, iBinX, 0, nBinsY);
      for(int iBinY=0;iBinY<=nBinsY;++iBinY){
	normalisedValue = hLayer->GetBinContent(iBinX, iBinY)/norm;
	if(normalisedValue<threshold) normalisedValue = 0.0;
	//else normalisedValue = log(normalisedValue) - log(threshold);
	hLayer->SetBinContent(iBinX, iBinY, normalisedValue);
	hLayer->SetBinError(iBinX, iBinY, 0);
	int iBinZ = hAllLayers->GetZaxis()->FindBin(iLayer);
	hAllLayers->SetBinContent(iBinX, iBinY, iBinZ, normalisedValue);
	hAllLayers->SetBinError(iBinX, iBinY, iBinZ, 0);
      }
    }    
  }

  std::vector<std::string> names = {"h2DDeltaVsDelta_MB2", "h2DDeltaVsDelta_RB1In",
				    "h2DDeltaVsDelta_RB1Out", "h2DDeltaVsDelta_RB2In",
				    "h2DDeltaVsDelta_RB2Out", "h2DDeltaVsDelta_RE23",
				    "h2DDeltaVsDelta_xy", "h2DDeltaVsDelta_xy_badReco"};
  for(auto & hName: names){
    TH2F *hDelta = get2DHistogram(hName, true);
    if(hDelta)hDelta->SetDirectory(myDirCopy);
  }
  
  return;
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void GMTHistograms::plotQuantiles(const std::string & hName){

  TH2F* h2D = get2DHistogram(hName);
  if(!h2D) return;
	
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  double xQ[1], yQ[1];
  xQ[0] = 0.20;

  TH1D *h1D = h2D->ProjectionX();
  h1D->SetXTitle("p_{T}^{GEN}");
  h1D->SetYTitle(TString::Format("p_{T}^{REC}@ CL %d",(int)(xQ[0]*100)));
  h1D->Reset();

  for(int iBinX=1;iBinX<h2D->GetNbinsX();++iBinX){
    TH1D *hPtRec = h2D->ProjectionY("hPtRec",iBinX,iBinX);
    if(hPtRec->Integral()<1) continue;	
    hPtRec->GetQuantiles(1, yQ, xQ);
    h1D->SetBinContent(iBinX, yQ[0]);
  }

  h1D->GetXaxis()->SetRangeUser(0,30);
  h1D->Draw();
  h1D->SetDirectory(myDirCopy);  
  std::cout<<"Fit range: 0-10"<<std::endl;
  h1D->Fit("pol1","W","",5, 10);
  std::cout<<"Fit range: 10-20"<<std::endl;
  h1D->Fit("pol1","W","",10, 10);
  std::cout<<"Fit range: 20-50"<<std::endl;
  h1D->Fit("pol1","W","",20, 50);
  std::cout<<"Fit range: 5-40"<<std::endl;
  h1D->Fit("pol2","W","",5, 40);

  c->Print(TString::Format("fig_png/Quantile%d_%s.png",(int)(xQ[0]*100), hName.c_str()));
}
// ////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////
