#include <iostream>
#include <cmath>

#include "OTFHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"

#include "utilsL1RpcStyle.h"

int nPtBins = 32;
const float OTFHistograms::ptBins[33]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
  		 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
  		 160. };

const int OTFHistograms::color[6] = {kBlack, kRed, kGreen, kBlue, kMagenta, kTeal};
const int OTFHistograms::ptCutsGmt[4] = {0, 16, 18, 19};
const int OTFHistograms::ptCutsOtf[4] = {0, 15, 17, 18};


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
bool OTFHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill2DHistogram(name,val1,val2,weight)){
		if(name.find("Pt")!=std::string::npos) hTemplateName = "h2DPt";
		if(name.find("EtaHit")!=std::string::npos) hTemplateName = "h2DEtaHit";
		if(name.find("PhiHit")!=std::string::npos) hTemplateName = "h2DPhiHit";
		if(name.find("EtaVx")!=std::string::npos) hTemplateName = "h2DEtaVx";
		if(name.find("PhiVx")!=std::string::npos) hTemplateName = "h2DPhiVx";
		if(name.find("RateTot")!=std::string::npos) hTemplateName = "h2DRateTot";
		if(name.find("RateVsEta")!=std::string::npos) hTemplateName = "h2DRateVsEta";
		std::cout<<"Adding histogram: "<<name<<std::endl;
		this->add2DHistogram(name,"",
				this->get2DHistogram(hTemplateName)->GetNbinsX(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
				this->get2DHistogram(hTemplateName)->GetNbinsY(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmax(),
				file_);
		return AnalysisHistograms::fill2DHistogram(name,val1,val2,weight);
	}
	return true;
}



void OTFHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
 //Make template histos
 add2DHistogram("h2DPt","",150,0,150,2,-0.5,1.5,file_);
 add2DHistogram("h2DEtaHit","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHit","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);
 add2DHistogram("h2DEtaVx","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiVx","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);
 //Rate histos
 add2DHistogram("h2DRateTot","",400,1,200,142,0,142,file_);
 add2DHistogram("h2DRateVsEta","",400,1,200,2*26,0.8,1.25,file_);

   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(int nRuns, float weight){
  
  utilsL1RpcStyle()->cd();

  plotEffPanel("Gmt");
  plotEffVsEta("Gmt");
  plotEffVsVar("Gmt","EtaVx");
  plotEffVsVar("Gmt","EtaHit");
  plotEffVsVar("Gmt","PhiVx");
  plotEffVsVar("Gmt","PhiHit");
  
  plotEffPanel("Otf");
  plotEffVsEta("Otf");
  plotEffVsVar("Otf","EtaVx");
  plotEffVsVar("Otf","EtaHit");
  plotEffVsVar("Otf","PhiVx");
  plotEffVsVar("Otf","PhiHit");

  plotEffPanel("Rpc");
  plotEffVsEta("Rpc");
  plotEffVsVar("Rpc","EtaVx");
  plotEffVsVar("Rpc","EtaHit");
  plotEffVsVar("Rpc","PhiVx");
  plotEffVsVar("Rpc","PhiHit");

  plotEffPanel("Other");
  plotEffVsEta("Other");
  plotEffVsVar("Other","EtaVx");
  plotEffVsVar("Other","EtaHit");
  plotEffVsVar("Other","PhiVx");
  plotEffVsVar("Other","PhiHit");
  
  plotOtfVsGmt(16);
  plotOtfVsGmt(18);
  plotOtfVsGmt(19);

  //plotRate("Tot");
  plotRate("VsEta");
  //plotEffVsRate(18);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* OTFHistograms::Integrate(TH1 * histoD) {

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
  const int *ptCuts = ptCutsOtf;
  if(sysType=="Gmt") ptCuts = ptCutsGmt;

  for (int icut=0; icut <=3;++icut){
    float ptCut = OTFHistograms::ptBins[ptCuts[icut]];
    hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName.Data());
    if(!h2D) return;
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
    TString nameCut = TString::Format("%d", (int)OTFHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/PanelVsPt_%s.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/PanelVsPt_%s.png",sysType.c_str()).Data());
}


void OTFHistograms::plotEffVsVar(const std::string & sysType,
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
    const int *ptCuts = ptCutsOtf;
    if(sysType=="Gmt") ptCuts = ptCutsGmt;

    for (int icut=0; icut <=1;++icut){
    	float ptCut = OTFHistograms::ptBins[ptCuts[icut]];
    	hName = "h2D"+sysType+varName+std::to_string((int)ptCut);
    	TH2F* h2D = this->get2DHistogram(hName.Data());
	if(!h2D) return;
    	TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    	TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    	hNum->Print();
    	hDenom->Print();
    	hDenom->Add(hNum);
    	hDenom->Print();
    	TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
    	hEff->SetStats(kFALSE);
    	hEff->SetMinimum(0.8);
    	hEff->SetMaximum(1.04);
    	hEff->SetMarkerStyle(21+icut);
    	hEff->SetMarkerColor(color[icut]);
    	hEff->SetXTitle(varName.c_str());
    	hEff->SetYTitle("Efficiency");
    	if (icut==0)hEff->DrawCopy("E0");
    	else hEff->DrawCopy("same E0");
    	hEff->Print();
    	std::string nameCut = std::to_string((int)OTFHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    	if (icut==0) nameCut = "no p_{T} cut";
    	l.AddEntry(hEff,nameCut.c_str());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
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

  const int *ptCuts = ptCutsOtf;
  if(sysType=="Gmt" || sysType=="Rpc" || sysType=="Other") ptCuts = ptCutsGmt;

  int iCut = 2;
 std::string hName = "";
  for (int iType=0; iType<3;++iType){
    float ptCut = OTFHistograms::ptBins[ptCuts[iCut]];
    hName = "h2D"+sysType+"Type" + std::to_string(iType) + "EtaVx"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
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
    std::string nameCut = std::to_string((int)OTFHistograms::ptBins[ptCuts[iCut]])+" GeV/c";
    if (iType==0) nameCut = "p_{T}^{#mu}>p_{T}^{cut} + 20 GeV";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 5 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c";
    l.AddEntry(hEff,nameCut.c_str());
  }
  l.SetHeader(TString::Format("p_{T}^{cut} = %d  GeV/c",(int)OTFHistograms::ptBins[ptCuts[iCut]]).Data());
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

void OTFHistograms::plotOtfVsGmt(int iPtCut){

  float ptCut = ptBins[iPtCut];

  TCanvas* c = new TCanvas(TString::Format("OtfVsGmt_%d",(int)ptCut).Data(),
			   TString::Format("OtfVsGmt_%d",(int)ptCut).Data(),
			   460,500);

  TLegend l(0.1995614,0.7139831,0.4385965,0.8665254,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(0,1);

  std::string hName = "h2DGmtPt"+std::to_string((int)ptCut);
  TH2F* h2D = this->get2DHistogram(hName);
  if(!h2D) return;
  float  effGmt = getEfficiency(h2D,ptCut);
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffGmt =DivideErr(hNum,hDenom,"hEffGmt","B");

  float effOtfMatch = 0.0;
  int match = -1;
  float delta = 999.0;
  TH1D *hEffOtf = 0;
  for (int iCut=-2; iCut<=2;++iCut){
    std::string hName = "h2DOtfPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    TH2F* h2D = this->get2DHistogram(hName);
    float  effOtf = getEfficiency(h2D,ptCut);
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TH1D* hEffOtfTmp =DivideErr(hNum,hDenom,"hEffOtfTmp","B");

    std::cout<<ptCut<<" "<<ptBins[iPtCut+iCut]<<" Otf: "<<effOtf
    		  <<" Gmt: "<<effGmt<<" diff: "<<fabs(effGmt-effOtf)<<std::endl;
    if(fabs(effGmt-effOtf)<delta){
      match = iPtCut+iCut;
      delta = fabs(effGmt-effOtf);
      effOtfMatch = effOtf;
      hEffOtf = hEffOtfTmp;
    }
  }
  std::cout<<"Gmt eff: "<<effGmt<<std::endl;
  std::cout<<"Otf eff: "<<effOtfMatch<<std::endl;

  hEffGmt->SetMarkerStyle(23);
  hEffGmt->SetMarkerColor(13);

  hEffOtf->SetXTitle("muon p_{T} [GeV/c]");
  hEffOtf->SetYTitle("Efficiency");
  hEffOtf->SetMinimum(0.0001);
  hEffOtf->SetMaximum(1.04);
  hEffOtf->GetXaxis()->SetRange(4,100);
  hEffOtf->SetMarkerStyle(8);
  hEffOtf->SetMarkerColor(3);
  hEffOtf->SetStats(kFALSE);
  hEffOtf->DrawCopy();
  hEffGmt->DrawCopy("same");

  std::string nameCut(TString::Format("%1.0f GeV/c",ptBins[match]).Data());
  l.AddEntry(hEffOtf,("Otf, "+nameCut).c_str());
  std::string tmp = "Gmt, %1.0f GeV/c";
  if( int(ptBins[match]*10)%10==5)  tmp = "Gmt, %1.1f GeV/c";
  l.AddEntry(hEffGmt,TString::Format(tmp.c_str(),ptCut).Data());
  l.DrawClone();

  TLine aLine(0,0,0,0);
  aLine.SetLineColor(2);
  aLine.SetLineWidth(3);
  aLine.DrawLine(ptCut,0,ptCut,1.04);

  c->Print(TString::Format("fig_eps/OtfVsGmt_%d.eps",(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OtfVsGmt_%d.png",(int)ptCut).Data());

  c->SetLogy();
  c->Print(TString::Format("fig_eps/OtfVsGmt_%d_log.eps",(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OtfVsGmt_%d_log.png",(int)ptCut).Data());

}

TH2F* OTFHistograms::makeRateWeights(TH2 *hOrig){

  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
  fIntVxMuRate->SetParameters(-0.235801, -2.82346, 17.162);

  TH2F *hWeights = (TH2F*)hOrig->Clone("hWeights");
  hWeights->Reset();
  TH1D *hPtGen = hOrig->ProjectionX("hPtGen");

  int iBin, nEvInBin;
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

TH1* OTFHistograms::getRateHisto(std::string sysType, 
				 std::string type){

  std::string hName = "h2DRate"+type+sysType;
  if(sysType=="Vx") hName = "h2DRate"+type+"Gmt";

  std::cout<<hName<<std::endl;

  TH2F* h2D = (TH2F*)this->get2DHistogram(hName)->Clone("h2D");
  if(!h2D) return 0;
  TH2F *hWeights = makeRateWeights(h2D);
  h2D->Multiply(hWeights);

  TH1D *hRate = h2D->ProjectionY("hRate");
  if(sysType=="Vx") hRate = h2D->ProjectionX("hRate");

  hRate->SetYTitle("Arbitrary units");
  hRate->SetLineWidth(3);

  if(type=="VsEta") return (TH1*)hRate->Clone("hRateClone");
  return Integrate(hRate);    
}

void OTFHistograms::plotRate(std::string type){

  TH1 *hRateVx = getRateHisto("Vx",type);
  TH1 *hRateGmt = getRateHisto("Gmt",type);
  TH1 *hRateOtf = getRateHisto("Otf",type);
  if(!hRateVx || !hRateGmt || !hRateOtf) return;

  hRateGmt->Print();
  hRateOtf->Print();

  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
  fIntVxMuRate->SetParameters(-0.235801, -2.82346, 17.162);
  fIntVxMuRate->SetLineColor(6);
 
  hRateVx->SetLineWidth(3);
  hRateGmt->SetLineWidth(3);
  hRateOtf->SetLineWidth(3);

  hRateGmt->SetLineColor(2);
  hRateOtf->SetLineColor(4);

  hRateGmt->SetLineStyle(2);
  
  TCanvas* c = new TCanvas("cRate","Rate",460,500);
  c->SetLogy(1);  
  c->SetGrid(1,1);
  
  TLegend *leg = new TLegend(0.6074561,0.6800847,0.8464912,0.8728814,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type=="Tot"){
    int iBinMin = hRateVx->FindBin(1); 
    int iBinMax = hRateVx->FindBin(100); 
    hRateVx->GetXaxis()->SetRange(iBinMin,iBinMax);
    hRateVx->SetMinimum(8E6);
    hRateVx->SetXTitle("p_{T}^{cut} [GeV/c]");
    hRateVx->Draw();
    hRateGmt->Draw("same");
    hRateOtf->Draw("same");
    leg->AddEntry(hRateVx,"#mu rate@Vx");
  }
  if(type=="VsEta"){
    c->SetLogy(0);  
    hRateGmt->SetXTitle("muon #eta");
    hRateGmt->Draw();
    hRateOtf->Draw("same");
  }
  
  leg->AddEntry(hRateGmt,"GMT");
  leg->AddEntry(hRateOtf,"OTF");
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());

}

void OTFHistograms::plotEffVsRate(int iPtCut){

  TCanvas* c = new TCanvas("cEffVsRate","EffVsRate",460,500);
  c->SetGrid(1,1);
  
  TLegend *leg = new TLegend(0.2149123,0.6716102,0.4539474,0.8644068,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  TH1 *hRateGmt = getRateHisto("Gmt","Tot");
  TH1 *hRateOtf = getRateHisto("Otf","Tot");
  if(!hRateGmt || !hRateOtf) return;
    
  TGraph *grGmt = new TGraph();
  TGraph *grOtf = new TGraph();

  grGmt->SetMarkerColor(2);
  grGmt->SetMarkerStyle(20);
  grGmt->SetMarkerSize(1.3);

  grOtf->SetMarkerColor(4);
  grOtf->SetMarkerStyle(22);
  grOtf->SetMarkerSize(1.3);

  float effNom, rateNom;
  for(int iCut=-2;iCut<=2;++iCut){
    ///
    std::string hName = "h2DOtfPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    float effOtf = getEfficiency(h2D,ptBins[iPtCut]);    
    ///
    hName = "h2DGmtPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    float effGmt = getEfficiency(h2D,ptBins[iPtCut]);    
    ///
    int iBin = hRateGmt->FindBin(ptBins[iPtCut+iCut]);
    float rateGmt = hRateGmt->GetBinContent(iBin);
    grGmt->SetPoint(iCut+2,effGmt,rateGmt);
    if(iCut==0){
      effNom = effGmt;
      rateNom = rateGmt;
    }
    ///
    float rateOtf = hRateOtf->GetBinContent(iBin);
    grOtf->SetPoint(iCut+2,effOtf,rateOtf);
    std::cout<<ptBins[iPtCut+iCut]<<" effOtf: "<<effOtf<<" rateOtf: "<<rateOtf
	     <<ptBins[iPtCut+iCut]<<" effGmt: "<<effGmt<<" rateGmt: "<<rateGmt<<std::endl;
  }

  Double_t maxY, minY, tmp;
  grGmt->GetPoint(0,tmp,maxY);
  grOtf->GetPoint(4,tmp,minY);
  maxY*=1.3;
  minY*=0.7;
  
  TH1F *hFrame = new TH1F("hFrame",";Efficiency;Rate",2,0.80,1.0);
  hFrame->SetMinimum(minY);
  hFrame->SetMaximum(maxY);
  hFrame->GetYaxis()->SetTitleOffset(1.7);
  hFrame->GetYaxis()->SetLabelSize(0.04);
  hFrame->GetXaxis()->SetLabelSize(0.04);
  hFrame->SetXTitle(TString::Format("Efficiency for p_{T}^{gen}>%d GeV/c",(int)ptBins[iPtCut]));
  hFrame->Draw();

  grGmt->Draw("P");
  grOtf->Draw("P");

  TMarker *aMarker = new TMarker(effNom,rateNom,4);
  aMarker->SetMarkerSize(2.1);
  aMarker->Draw();

  leg->AddEntry(aMarker,"Nominal GMT cut","p");
  leg->AddEntry(grGmt,"GMT","p");
  leg->AddEntry(grOtf,"OTF","p");
  leg->Draw();

  c->Print(TString::Format("fig_eps/RateVsEff_%d.eps",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.png",(int)ptBins[iPtCut]).Data());
  c->Print(TString::Format("fig_png/RateVsEff_%d.C",(int)ptBins[iPtCut]).Data());

}


float  OTFHistograms::getEfficiency(TH2F *h2D, float ptCut){

  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
  hDenom->Add(hNum);
  TH1D* hEffTmp =DivideErr(hNum,hDenom,"hEffTmp","B");
  //Mean eff above pt cut
  int binLow = hEffTmp->FindBin(ptCut+5);
  int binHigh = hEffTmp->FindBin(ptCut+20);
  float range = hEffTmp->GetBinLowEdge(binHigh+1) - hEffTmp->GetBinLowEdge(binLow);
  float eff = hEffTmp->Integral(binLow,binHigh,"width")/range;
  
  return eff;
}
