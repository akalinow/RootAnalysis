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
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooJohnson.h"
using namespace RooFit ;


#include "utilsL1RpcStyle.h"

int nPtBins = 32;
const float GMTHistograms::ptBins[36]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
		 10., 12., 14., 16., 18., 20., 22,  24, 26, 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
		 160., 200};

const int GMTHistograms::color[6] = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
//Single mu
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
  if(name.find("DiMuonMass")!=std::string::npos) templateName = "h1DDiMuonMassTemplate";
  if(name.find("Pt")!=std::string::npos) templateName = "h2DPtTemplate";
  if(name.find("HighPt")!=std::string::npos) templateName = "h2DHighPtTemplate";
  if(name.find("PtRecVsPtGen")!=std::string::npos) templateName = "h2DPtVsPtTemplate";


  
  if(name.find("EtaHit")!=std::string::npos) templateName = "h2DEtaHitTemplate";
  if(name.find("PhiHit")!=std::string::npos) templateName = "h2DPhiHitTemplate";
  if(name.find("EtauGMT")!=std::string::npos) templateName = "h2DEtauGMTTemplate";
  if(name.find("PhiuGMT")!=std::string::npos) templateName = "h2DPhiuGMTTemplate";
  if(name.find("PtRecVsPtOMTF")!=std::string::npos) templateName = "h2DPtRecVsPtOMTFTemplate";

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
 add1DHistogram("h1DDiMuonMassTemplate", "", 80, 70, 110, file_);
 


 ///Efficiency histos
 add2DHistogram("h2DPtTemplate","",150,0,150,2,-0.5,1.5,file_);
 add2DHistogram("h2DHighPtTemplate","",50,50,550,2,-0.5,1.5,file_);
 add2DHistogram("h2DPtRecVsPtOMTFTemplate", "", 100, 0, 200, 100, 0, 200, file_);
 add2DHistogram("h2DPtVsPtTemplate","",404,0,202,404,0,202,file_);

 add2DHistogram("h2DEtaHitTemplate","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHitTemplate","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtauGMTTemplate","",60,-3,3,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiuGMTTemplate","",4*32,-3.2,3.2,2,-0.5,1.5,file_);

 add2DHistogram("h2DQualityTemplate","",201,-0.5,200.5,2,-0.5,1.5,file_);

 //Rate histos
 add2DHistogram("h2DRateTotTemplate","",404,1,202,404,1,202,file_);
 add2DHistogram("h2DRateVsEtaTemplate","",404,1,202,60,-3,3,file_);

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
  gErrorIgnoreLevel = kError;
  
  //Panel with many turn-on curves
  plotEffPanel("uGMT");
  //Panel with many turn-on curves for high pT range
  plotEffPanel("uGMT", true);

  
  //Efficiency as a function of ete.
  //Lines for selected points on the turn on curve shown
  plotEffVsEta("uGMT");
  plotEffVsVar("uGMT", "Eta");
  plotEffVsVar("uGMT", "Phi");
  plotSingleHistogram("h2DuGMTPtRecVsPtOMTF");
  plotSingleHistogram("h1DDiMuonMass"); 
  //Turn on curves for many pT thresholds.
  ///Lines for reference - Phase2 uGMT, and other algorithm shown
  for(int iPtCode=1;iPtCode<=30;++iPtCode){
      plotGMTVsOther(iPtCode,"uGMT");
  }
   
  
  
  //Event rate plot for selected type
  //Uses some old rate parametrisation for the single muons sample
  ///Works also with neutrino sample
  plotRate("Tot");
   
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
    hEff->SetXTitle("reco. muon p_{T} [GeV/c]");
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
    hName = "h2D"+sysType+"Type" + std::to_string(iType) + "EtauGMT"+std::to_string((int)ptCut);
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
    hEff->SetXTitle("reco muon #eta");
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
void GMTHistograms::plotGMTVsOther(int iPtCut,
				     const std::string sysType){

  float ptCut = ptBins[iPtCut];

  TCanvas* c = new TCanvas(TString::Format("uGMTVsOther_%d",(int)ptCut).Data(),
			   TString::Format("uGMTVsOther_%d",(int)ptCut).Data(),
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

  hName = "h2DuGMTPt"+std::to_string((int)ptCut);
  h2D = get2DHistogram(hName);
  hNum = h2D->ProjectionX("hNum",2,2);
  hDenom = h2D->ProjectionX("hDenom",1,1);    
  hDenom->Add(hNum);

  TH1D* hEffGMT =DivideErr(hNum,hDenom,"hEffuGMTTmp","B");
  hEffGMT->SetXTitle("reco. muon p_{T} [GeV/c]");
  hEffGMT->SetYTitle("Efficiency");
  hEffGMT->SetMaximum(1.04);
  hEffGMT->GetXaxis()->SetRange(2,100);
  hEffGMT->SetMarkerStyle(8);
  hEffGMT->SetMarkerColor(1);
  hEffGMT->SetStats(kFALSE);
  hEffGMT->DrawCopy();
  hEffOther->DrawCopy("same");

  std::string tmp = "p_{T} #geq ";
  if(int(ptCut*10)%10==5) tmp += "%1.1f GeV/c";
  else   tmp += "%1.0f GeV/c";
  l.AddEntry((TObject*)0, TString::Format(tmp.c_str(),ptCut).Data(), "");
  l.AddEntry((TObject*)0, "", "");
  l.AddEntry(hEffOther, sysType.c_str());
  l.AddEntry((TObject*)0, "", "");
  l.AddEntry(hEffGMT, "Phase2 GMT");
  l.DrawClone();

  TLine aLine(0,0,0,0);
  aLine.SetLineColor(2);
  aLine.SetLineWidth(3);
  aLine.DrawLine(ptCut,0,ptCut,1.04);

  c->Print(TString::Format("fig_eps/uGMTVs%s_%d.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/uGMTVs%s_%d.png",sysType.c_str(),(int)ptCut).Data());

  c->SetLogy();
  c->Print(TString::Format("fig_eps/uGMTVs%s_%d_log.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/uGMTVs%s_%d_log.png",sysType.c_str(),(int)ptCut).Data());
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
  
  TH1 *hRateuGMT_emu = getRateHisto("uGMT",type);
  TH1 *hRateVx = getRateHisto("Vx",type);

  if(!hRateVx || !hRateuGMT_emu) return;

  hRateVx->SetLineWidth(3); 
  hRateVx->SetLineColor(1);
  
  hRateuGMT_emu->SetLineWidth(3);
  hRateuGMT_emu->SetLineColor(2);
  hRateuGMT_emu->SetLineStyle(2);
 
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
    hRateVx->SetMinimum(1E1);
    hRateVx->SetMaximum(2E5);
    hRateVx->SetXTitle("p_{T}^{cut} [GeV/c]");
    c->SetLogy();
    c->SetGrid(1,0);
    hRateVx->Draw();
    hRateuGMT_emu->DrawCopy("same");
   
    std::cout<<"Rate uGMT @ 20 GeV: "<< hRateuGMT_emu->GetBinContent(hRateuGMT_emu->FindBin(20-0.01))<<std::endl;
  }  
  else if(type.find("VsEta")!=std::string::npos){
    c->SetLogy(0);
    hRateuGMT_emu->SetXTitle("muon #eta");
    double max = hRateuGMT_emu->GetMaximum();
    hRateuGMT_emu->SetMaximum(1.5*max);
    hRateuGMT_emu->Draw();
  }
  if(type=="VsPt"){
    c->SetLogy(1);
    hRateuGMT_emu->SetXTitle("p_{T}^{gen} [GeV/c]");
    hRateuGMT_emu->SetAxisRange(2,100);
    double max = hRateuGMT_emu->GetMaximum();
    hRateuGMT_emu->SetMaximum(10*max);
    hRateuGMT_emu->Draw();
  }

  leg->AddEntry(hRateuGMT_emu,"uGMT Phase 2");
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());
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
///////////////////////////////////////////////////////////////
void GMTHistograms::plotSingleHistogram(std::string hName){

  TH2F* h2D = get2DHistogram(hName);
  TH1F* h1D = get1DHistogram(hName);
  if(!h2D && !h1D) return;
	
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   800,800);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(h2D) {
    h2D->SetDirectory(myDirCopy);
    h2D->SetLineWidth(3);
    h2D->Scale(1.0/h2D->Integral());
    h2D->SetXTitle("p_{T}^{Reco}");
    h2D->SetYTitle("p_{T}^{OMTF}");
    h2D->GetYaxis()->SetTitleOffset(1.4);
    h2D->SetStats(kFALSE);
    gStyle->SetPalette(kRainBow);
    h2D->Draw("box colz");
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
  if(h1D) {
    h1D->SetDirectory(myDirCopy);
    h1D->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));    
    h1D->GetXaxis()->SetRange(1,h1D->GetNbinsX()+1);
    h1D->SetXTitle("Z(#mu^{+}#mu^{-}) [GeV]");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    gStyle->SetOptStat(0) ;
    gStyle->SetPalette(1) ;
    RooRealVar mass("mass", "Z(#mu^{+}#mu^{-}) (GeV/c^{2})", 70, 110);
    RooDataHist dh("dh", "dh", mass, Import(*h1D));
    RooRealVar mu("mu", "mu", 91.18, 90, 93);
    RooRealVar lambda ("lambda", "lambda",  0.5, 0, 10);
    RooRealVar gamma ("gamma", "gamma",  0., 0.00, 0.1);
    RooRealVar delta ("delta", "delta", 1., 0, 20);
    RooJohnson john ("john", "john", mass, mu, lambda, gamma, delta);    
    RooRealVar conts ("conts","conts", -10.0, 0.0);
    RooExponential expoBg ("expoBg","expoBG",mass,conts);
    RooRealVar nSig("nSig", "nSig",500,0,(int)h1D->GetEntries());
    RooRealVar nBkg("nBkg", "nBkg", 500,0,(int)h1D->GetEntries());
    RooAddPdf zpdf("zpdf","zpdf",RooArgList(john,expoBg), RooArgList(nSig,nBkg));
    RooFitResult* fitRes = zpdf.fitTo(dh,Save(true),NumCPU(8), RooFit::Minimizer("Minuit2", "Migrad"));
    fitRes->Print("v");
    RooPlot* zmassf = mass.frame(Title("Z(#mu^{+}#mu^{-}) (GeV/c^{2})"));
    dh.plotOn(zmassf,DataError(RooAbsData::SumW2));
    zpdf.plotOn(zmassf) ;
    double chisquare_mass = zmassf->chiSquare();
    std::cout<<"Chi square of mass fit is :   "<< chisquare_mass<<"\n";
    zpdf.plotOn(zmassf, RooFit::LineColor(kGreen),RooFit::Components("john"), RooFit::Name("signal"), LineWidth(2), LineStyle(4));
    zpdf.plotOn(zmassf,RooFit::LineColor(kRed),RooFit::Components("expoBg"), RooFit::Name("combinatorial"), LineWidth(2), LineStyle(6));
    TLegend *leg = new TLegend(0.25,0.7,0.45,0.9);
    leg->AddEntry(zmassf->findObject("signal"),"Z^{0}#rightarrow #mu^{+}#mu^{-}","l");
    leg->AddEntry(zmassf->findObject("combinatorial"),"Combi","l");
    zmassf->SetStats(0);
    zmassf->Draw();
    leg->Draw("same");
    //h1D->Draw("");
    //h1D->Print();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
