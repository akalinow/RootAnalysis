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
#include "TRandom3.h"

#include "utilsL1RpcStyle.h"

int nPtBins = 32;
const float OTFHistograms::ptBins[33]={0., 0.1,
  		 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8.,
  		 10., 12., 14., 16., 18., 20., 25., 30., 35., 40., 45.,
  		 50., 60., 70., 80., 90., 100., 120., 140.,
  		 160. };

const int OTFHistograms::color[6] = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
//Single mu
const int OTFHistograms::ptCutsGmt[4] =     {0, 14, 16, 18};
const int OTFHistograms::ptCutsOtf[4] =     {0, 14, 16, 18};
//const int OTFHistograms::ptCutsOtf[4] =     {0, 14, 16, 12};
//Di muon
//const int OTFHistograms::ptCutsGmt[4] =     {0, 5, 7, 14};
//const int OTFHistograms::ptCutsOtf[4] =     {0, 5, 7, 14};
///No scale shift
//const int OTFHistograms::ptCutsOtf[4] =     {0, 15, 16, 18};


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
OTFHistograms::~OTFHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool OTFHistograms::fill1DHistogram(const std::string& name, float val1, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill1DHistogram(name,val1,weight)){
	  if(name.find("DeltaEta")!=std::string::npos) hTemplateName = "h1DDeltaEtaTemplate";
	  if(name.find("DeltaPhi")!=std::string::npos) hTemplateName = "h1DDeltaPhiTemplate";
	  std::cout<<"Adding histogram: "<<name<<" with template: "<<hTemplateName<<std::endl;
	  this->add1DHistogram(name,"",
			       this->get1DHistogram(hTemplateName)->GetNbinsX(),
			       this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
			       this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
			       file_);
	  return AnalysisHistograms::fill1DHistogram(name,val1,weight);
	}
	return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool OTFHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill2DHistogram(name,val1,val2,weight)){
	  
		if(name.find("Pt")!=std::string::npos) hTemplateName = "h2DPtTemplate";
		if(name.find("EtaHit")!=std::string::npos) hTemplateName = "h2DEtaHitTemplate";
		if(name.find("PhiHit")!=std::string::npos) hTemplateName = "h2DPhiHitTemplate";
		if(name.find("EtaVx")!=std::string::npos) hTemplateName = "h2DEtaVxTemplate";
		if(name.find("PhiVx")!=std::string::npos) hTemplateName = "h2DPhiVxTemplate";
		if(name.find("Quality")!=std::string::npos) hTemplateName = "h2DQualityTemplate";
		if(name.find("RateTot")!=std::string::npos) hTemplateName = "h2DRateTotTemplate";
		if(name.find("RateVsEta")!=std::string::npos) hTemplateName = "h2DRateVsEtaTemplate";
		if(name.find("RateVsPt")!=std::string::npos) hTemplateName = "h2DRateVsPtTemplate";
		if(name.find("RateVsQuality")!=std::string::npos) hTemplateName = "h2DRateVsQualityTemplate";                                              
		if(name.find("DeltaPhi")!=std::string::npos) hTemplateName = "h2DDeltaPhiTemplate";
		if(name.find("DeltaPt")!=std::string::npos) hTemplateName = "h2DDeltaPtTemplate";
		if(name.find("GhostsVsProcessor")!=std::string::npos) hTemplateName = "h2DGhostsVsProcessorTemplate";
		std::cout<<"Adding histogram: "<<name<<" with template: "<<hTemplateName<<std::endl;
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
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

 //Make template histos
 add1DHistogram("h1DDeltaEtaTemplate","",100,0,5.0,file_);
   
 add1DHistogram("h1DDeltaPhiTemplate","",5*32,-M_PI,M_PI,file_);

 ///Efficiency histos
 add2DHistogram("h2DPtTemplate","",150,0,150,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaHitTemplate","",8*26,0.8,1.25,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiHitTemplate","",5*32,-M_PI,M_PI,2,-0.5,1.5,file_);

 add2DHistogram("h2DEtaVxTemplate","",40,-1.6,1.6,2,-0.5,1.5,file_);//Full detector
 //add2DHistogram("h2DEtaVxTemplate","",20,0.83,1.24,2,-0.5,1.5,file_);//Overlap 
 //add2DHistogram("h2DEtaVxTemplate","",8*25,-0.1,0.85,2,-0.5,1.5,file_);//Barrel
 //add2DHistogram("h2DEtaVxTemplate","",8*25,1.25,2.65,2,-0.5,1.5,file_);//Endcap
 
 add2DHistogram("h2DPhiVxTemplate","",4*32,-0.2,3.2,2,-0.5,1.5,file_);
 
 add2DHistogram("h2DQualityTemplate","",1000,1.5,1001.5,2,-0.5,1.5,file_);
 
 //Rate histos
 add2DHistogram("h2DRateTotTemplate","",400,1,201,142,0,142,file_);

 add2DHistogram("h2DRateVsEtaTemplate","",400,1,201,32*2,-1.6,1.6,file_);//Full detector
 //add2DHistogram("h2DRateVsEtaTemplate","",400,1,201,25,0.8,1.25,file_);//Overlap
 //add2DHistogram("h2DRateVsEtaTemplate","",400,1,201,25,-0.1,0.85,file_);//Barrel
 //add2DHistogram("h2DRateVsEtaTemplate","",400,1,201,25,1.25,2.7,file_);//Encap
 
 //add2DHistogram("h2DDeltaPhiTemplate","",40,-M_PI,M_PI,2,-0.5,1.5,file_);
 add2DHistogram("h2DDeltaPhiTemplate","",30,-1,1,2,-0.5,1.5,file_);
 add2DHistogram("h2DDeltaPtTemplate","",21,-0.5,20.5,2,-0.5,1.5,file_);

 add2DHistogram("h2DRateVsPtTemplate","",400,1,201,60,0,30,file_);
 add2DHistogram("h2DRateVsQualityTemplate","",400,1,201,17,1.5,18.5,file_);
 add2DHistogram("h2DGhostsVsProcessorTemplate","",6,-0.5,5.5,5,-0.5,4.5,file_);

 histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();
  
  utilsL1RpcStyle()->cd();

  //plotRate("Tot");
  //return;

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

  plotOtfVsGmt(18,"Rpc");
  plotOtfVsGmt(18,"Other");

  plotRate("Tot");
  plotRate("VsEta");
  plotRate("VsPt");
  plotRate("VsQuality");

  plotEffVsRate(18);
  /*
  plotGhostHistos("Gmt","");
  plotGhostHistos("Otf","");
  plotGhostHistos("Otf","SS");
  plotGhostHistos("Otf","OS");
  */
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeDiMuonHistograms(int nRuns, float weight){

  
  plotEffVsVar("OtfDi","DeltaPhi");
  plotEffVsVar("OtfDi","DeltaPt");
  plotEffVsVar("OtfDi","PhiHit");
  plotVar("OtfDiMuon1","DeltaEta");
  plotVar("OtfDiMuon2","DeltaEta");

  /*
  plotEffPanel("GmtiMuon0");
  plotEffPanel("GmtiMuon1");
  
  plotEffPanel("OtfiMuon0");
  plotEffPanel("OtfiMuon1");
  
  plotVar("GmtiMuon00","DeltaEta");
  plotVar("OtfiMuon00","DeltaEta");

  plotVar("GmtiMuon01","DeltaEta");
  plotVar("OtfiMuon01","DeltaEta");
  */
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::finalizeHistograms(){
  finalizeHistograms(0,1.0);
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
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
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
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OTFHistograms::plotEffPanel(const std::string & sysType){

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
  const int *ptCuts = ptCutsOtf;
  if(sysType.find("Gmt")!=std::string::npos || 
     sysType=="Rpc" || 
     sysType=="Other") ptCuts = ptCutsGmt;

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
    hEff->GetXaxis()->SetRange(0,50);
    //hEff->GetXaxis()->SetRange(4,20);
    
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
  c->SetLogy();
  c->Print(TString::Format("fig_eps/PanelVsPt_%s_log.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/PanelVsPt_%s_log.png",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OTFHistograms::plotVar(const std::string & sysType,
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

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
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

  for (int icut=0; icut <=2;++icut){
    if(icut==1) continue;
    float ptCut = OTFHistograms::ptBins[ptCuts[icut]];
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

    if(sysType=="OtfDi" && hName.Contains("DeltaPt")){
      hEff->SetMinimum(1E-2);
      c->SetLogy();
    }

    if (icut==0)hEff->DrawCopy("E0");
    else hEff->DrawCopy("same E0");
    std::string nameCut = std::to_string((int)OTFHistograms::ptBins[ptCuts[icut]])+" GeV/c";
    if (icut==0) nameCut = "no p_{T} cut";
    if(sysType=="OtfDi" && icut>0) nameCut = "opposite sign";
    if(sysType=="OtfDi" && icut==0) nameCut = "same sign";	
    
    l.AddEntry(hEff,nameCut.c_str());
  }
  l.DrawClone();
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
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
    if(iType==2) hNum->Scale(10.0);
    hDenom->Add(hNum);
    TH1D* hEff =DivideErr(hNum,hDenom,"Pt_Int","B");
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
    if (iType==0) nameCut = "p_{T}^{#mu}>24 GeV";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 10 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c (#epsilon #times 10)";
    l.AddEntry(hEff,nameCut.c_str());
  }
 ///OTF eta range used for generating patterns.
  TLine *aLine = new TLine(0,0,0,0);
  aLine->SetLineWidth(2);
  aLine->SetLineColor(2);
  aLine->DrawLine(0.83,0,0.83,1.0);
  aLine->DrawLine(-0.83,0,-0.83,1.0);
  aLine->DrawLine(1.24,0,1.24,1.0);
  aLine->DrawLine(-1.24,0,-1.24,1.0);

  l.SetHeader(TString::Format("p_{T}^{cut %s} = %d  GeV/c",sysType.c_str(),(int)OTFHistograms::ptBins[ptCuts[iCut]]).Data());
  l.DrawClone();

  c->Print(TString::Format("fig_eps/EffVsEta_%s.eps",sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVsEta_%s.png",sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OTFHistograms::plotOtfVsGmt(int iPtCut,
				 const std::string sysType){

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

  std::string hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
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
  //for (int iCut=-2; iCut<=2;++iCut){
  for (int iCut=0; iCut<1;++iCut){//TEST
    std::string hName = "h2DOtfPt"+std::to_string((int)ptBins[iPtCut+iCut]);
    TH2F* h2D = this->get2DHistogram(hName);
    float  effOtf = getEfficiency(h2D,ptCut);
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    TH1D* hEffOtfTmp =DivideErr(hNum,hDenom,"hEffOtfTmp","B");

    std::cout<<ptCut<<" "<<ptBins[iPtCut+iCut]<<" Otf: "<<effOtf
    		  <<" Gmt: "<<effGmt<<" diff: "<<fabs(effGmt-effOtf)<<std::endl;
    //if(fabs(effGmt-effOtf)<delta){
    if(effGmt-effOtf<delta){
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
  hEffOtf->GetXaxis()->SetRange(2,100);
  hEffOtf->SetMarkerStyle(8);
  hEffOtf->SetMarkerColor(3);
  hEffOtf->SetStats(kFALSE);
  hEffOtf->DrawCopy();
  hEffGmt->DrawCopy("same");

  std::string nameCut(TString::Format("%1.0f GeV/c",ptBins[match]).Data());
  l.AddEntry(hEffOtf,("Otf, "+nameCut).c_str());
  std::string tmp = sysType+", %1.0f GeV/c";
  if( int(ptBins[match]*10)%10==5)  tmp = sysType+", %1.1f GeV/c";
  l.AddEntry(hEffGmt,TString::Format(tmp.c_str(),ptCut).Data());
  l.DrawClone();

  TLine aLine(0,0,0,0);
  aLine.SetLineColor(2);
  aLine.SetLineWidth(3);
  aLine.DrawLine(ptCut,0,ptCut,1.04);

  c->Print(TString::Format("fig_eps/OtfVs%s_%d.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OtfVs%s_%d.png",sysType.c_str(),(int)ptCut).Data());

  c->SetLogy();
  c->Print(TString::Format("fig_eps/OtfVs%s_%d_log.eps",sysType.c_str(),(int)ptCut).Data());
  c->Print(TString::Format("fig_png/OtfVs%s_%d_log.png",sysType.c_str(),(int)ptCut).Data());

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TH1* OTFHistograms::getRateHisto(std::string sysType, 
				 std::string type){

  std::string hName = "h2DRate"+type+sysType;
  if(sysType=="Vx") hName = "h2DRate"+type+"Gmt";

  if(!this->get2DHistogram(hName)) return 0;
  TH2F* h2D = (TH2F*)this->get2DHistogram(hName)->Clone("h2D");

  if(!h2D) return 0;
  TH2F *hWeights = makeRateWeights(h2D);
  //h2D->Multiply(hWeights);//comment for TEST

  TH1D *hRate = h2D->ProjectionY("hRate");
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
void OTFHistograms::plotRate(std::string type){

  TH1 *hRateGmt = getRateHisto("Gmt",type);
  TH1 *hRateVx = getRateHisto("Vx",type);
  TH1 *hRateOtf = getRateHisto("Otf",type);

  if(!hRateVx || !hRateGmt || !hRateOtf) return;

 
  TF1 *fIntVxMuRate = new TF1("fIntVxMuRate","TMath::Power(x,[0]*TMath::Log(x))*TMath::Power(x,[1])*TMath::Exp([2])",1,1000);
  fIntVxMuRate->SetParameters(-0.235801, -2.82346, 17.162);
  fIntVxMuRate->SetLineColor(6);
 
  hRateVx->SetLineWidth(3);
  hRateGmt->SetLineWidth(3);
  hRateOtf->SetLineWidth(3);

  hRateVx->SetLineColor(1);
  hRateGmt->SetLineColor(2);
  hRateOtf->SetLineColor(4);

  hRateGmt->SetLineStyle(2);

  TCanvas* c = new TCanvas("cRate","Rate",1.5*420,1.5*500);
  c->SetLogy(1);  
  c->SetGrid(1,1);
  
  TLegend *leg = new TLegend(0.40,0.68,0.65,0.87,NULL,"brNDC");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type=="Tot"){
    int iBinMin = hRateGmt->FindBin(2); 
    int iBinMax = hRateGmt->FindBin(60);    
    hRateVx->GetXaxis()->SetRange(iBinMin,iBinMax);
    hRateGmt->GetXaxis()->SetRange(iBinMin,iBinMax);
    hRateOtf->GetXaxis()->SetRange(iBinMin,iBinMax);
    hRateVx->SetMinimum(1E3);
    hRateVx->SetMaximum(1E6);
    hRateGmt->SetXTitle("p_{T}^{cut} [GeV/c]");

    c->Divide(2);
    TPad *pad1 = (TPad*)c->GetPad(1);
    TPad *pad2 = (TPad*)c->GetPad(2);
    pad1->SetPad(0.01,0.29,0.99,0.99);
    pad2->SetPad(0.01,0.01,0.99,0.29);

     pad1->Draw();
     pad1->cd();
     pad1->SetLogy();
    //TEST hRateVx->Draw();
    hRateGmt->DrawCopy();//TEST
    //hRateGmt->DrawCopy("same");
    hRateOtf->DrawCopy("same");
    //fIntVxMuRate->Draw("same");
    //leg->AddEntry(hRateVx,"#mu rate@Vx");
    c->cd();
    pad2->Draw();
    pad2->cd();
    hRateOtf->SetYTitle("GMT/OMTF");
    hRateOtf->GetXaxis()->SetLabelSize(0.09);
    hRateOtf->GetYaxis()->SetLabelSize(0.09);
    hRateOtf->GetYaxis()->SetTitleSize(0.09);
    hRateOtf->GetYaxis()->SetTitleOffset(0.5);
    hRateOtf->Divide(hRateGmt);
    hRateOtf->DrawCopy();
    c->cd();
  }
  if(type=="VsEta"){
    c->SetLogy(0);  
    hRateGmt->SetXTitle("muon #eta");
    hRateGmt->Draw();
    hRateOtf->Draw("same");
  }
  if(type=="VsPt"){
    c->SetLogy(0);  
    hRateGmt->SetXTitle("p_{T}^{gen} [GeV/c]");
    hRateGmt->Draw();
    hRateOtf->Draw("same");    
  }
  if(type=="VsQuality"){
    c->SetLogy(0);  
    hRateOtf->SetXTitle("");
    //hRateGmt->Draw();
    hRateOtf->Draw();
    hRateOtf->Print();

    TH1D *hRateOtfSorted = (TH1D*)hRateOtf->Clone("hRateOtfSorted");
    TH1D *hEffOtfSorted = (TH1D*)hRateOtf->Clone("hEffOtfSorted");
    hRateOtfSorted->Clear();
    hEffOtfSorted->Clear();
    ///Calculate efficiency    
    TH2F* h2D = this->get2DHistogram("h2DOtfQuality16");
    if(!h2D) return;
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    hDenom->Add(hNum);
    hNum->Scale(1.0/hDenom->Integral());
    TH1D* hEff = hNum;
    hNum->Print();
    ////    
    std::map<float,std::string> rateMap;
    std::map<float,int> rateMapBin;
    ///Fill map with key = rate, value - bin label. The map will be sorted automatically by the rate value.
    ///Second map has bin number as value. Needed to extract efficiency values.
    TRandom3 aRndm;
    
    for(int iBin=1;iBin<hRateOtf->GetXaxis()->GetNbins();++iBin){
      float rate = hRateOtf->GetBinContent(iBin);
      rate+=0.001*aRndm.Uniform();///Randomize rate values, to avoid having two same keys.
      rateMap[rate] = hRateOtf->GetXaxis()->GetBinLabel(iBin);
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
      hRateOtfSorted->SetBinContent(iBin, it->first);
      hEffOtfSorted->SetBinContent(iBin, hEff->GetBinContent(rateMapBin[it->first]));      
      hRateOtfSorted->GetXaxis()->SetBinLabel(iBin, it->second.c_str());
      ++iBin;
    }
    hRateOtfSorted->GetXaxis()->SetRange(1,10);

    hRateOtfSorted->DrawCopy();

    hEffOtfSorted->SetLineColor(2);
    //hEffOtfSorted->Scale(1000);    
    //hEffOtfSorted->DrawCopy("same E0");
    hEffOtfSorted->GetXaxis()->SetRange(1,10);
    delete hRateOtfSorted;
    delete hEffOtfSorted;
  }

  
  leg->AddEntry(hRateGmt,"GMT");
  leg->AddEntry(hRateOtf,"OTF");
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());

}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
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
  for(int iCut=-3;iCut<=2;++iCut){
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
    std::cout<<ptBins[iPtCut+iCut]<<" effOtf: "<<effOtf<<" rateOtf: "<<rateOtf<<" "
	     <<ptBins[iPtCut+iCut]<<" effGmt: "<<effGmt<<" rateGmt: "<<rateGmt<<std::endl;
  }

  Double_t maxY, minY, tmp;
  float effMin=0.8, effMax=1.0;
  grGmt->GetPoint(0,tmp,maxY);
  effMax = tmp;
  grOtf->GetPoint(4,tmp,minY);
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
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OTFHistograms::plotGhostHistos(const std::string & sysType,
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
float  OTFHistograms::getEfficiency(TH2F *h2D, float ptCut){

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
