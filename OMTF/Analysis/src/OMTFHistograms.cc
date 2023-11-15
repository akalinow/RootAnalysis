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
#include "TEfficiency.h"

#include "utilsL1RpcStyle.h"

const std::vector<std::string> OMTFHistograms::algos = {"OMTF", "OMTFDisp", "OMTFDispU", "OMTFDispC", "GMT"};  
const std::vector<double> OMTFHistograms::ptBins = {1., 4, 4.5, 5, 5.5, 6, 6.5, 7, 8.5, 10, 
                                                    12, 14, 16, 18.5, 21, 23, 26, 28, 30, 32, 
                                                    36, 40, 48, 54, 60, 70, 82, 96, 114, 200, 99999};
const std::vector<double> OMTFHistograms::color = {kBlack, kBlue, kRed, kMagenta, kTeal, kGreen};
const std::vector<double> OMTFHistograms::iPtCuts = {0, 3, 10, 16};
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::OMTFHistograms(std::string fileName, int opt){ AnalysisHistograms::init(fileName); }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
OMTFHistograms::OMTFHistograms(TDirectory *myDir){ AnalysisHistograms::init(myDir); }
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
  if(name.find("Pt")!=std::string::npos) templateName = "h2DPtTemplate";
  if(name.find("HighPt")!=std::string::npos) templateName = "h2DHighPtTemplate";
  if(name.find("PtRecVsPtGen")!=std::string::npos) templateName = "h2DPtVsPtTemplate";  
  if(name.find("EtaVx")!=std::string::npos) templateName = "h2DEtaVxTemplate";
  if(name.find("PhiVx")!=std::string::npos) templateName = "h2DPhiVxTemplate";
  if(name.find("dxyVsPhiB")!=std::string::npos) templateName = "h2DdxyVsPhiBTemplate";
  else if(name.find("dxy")!=std::string::npos) templateName = "h2DdxyTemplate";
  if(name.find("dz")!=std::string::npos) templateName = "h2DdzTemplate";
  if(name.find("RateTot")!=std::string::npos) templateName = "h2DRateTotTemplate";
  if(name.find("RateVsEta")!=std::string::npos) templateName = "h2DRateVsEtaTemplate";
  if(name.find("RateVsPt")!=std::string::npos) templateName = "h2DRateVsPtTemplate";

  if(name.find("GenPt")!=std::string::npos) templateName = "h1DGenPtTemplate";
  if(name.find("GenEta")!=std::string::npos) templateName = "h1DGenEtaTemplate";
  if(name.find("GenDxy")!=std::string::npos) templateName = "h1DGenDxyTemplate";
  if(name.find("GenDz")!=std::string::npos) templateName = "h1DGenDzTemplate";
  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::defineHistograms(){

 if(!histosInitialized_){
 ///Efficiency histos
 add2DHistogram("h2DPtTemplate",";gen. muon p_{T} [GeV/c];Efficiency",60,0,60,2,-0.5,1.5,file_);
 add2DHistogram("h2DHighPtTemplate",";gen. muon p_{T} [GeV/c];Efficiency",50,50,550,2,-0.5,1.5,file_);
 add2DHistogram("h2DPtVsPtTemplate","",41,0,205,41,0,205,file_);
 add2DHistogram("h2DEtaVxTemplate",";gen. muon #eta;Efficiency",10,0.83,1.24,2,-0.5,1.5,file_);
 add2DHistogram("h2DPhiVxTemplate",";gen. muon #varphi;Efficiency",128,-3.2,3.2,2,-0.5,1.5,file_);
 add2DHistogram("h2DdxyTemplate",";d_{xy} [cm];Efficiency",50,0,500,2,-0.5,1.5,file_);
 add2DHistogram("h2DdzTemplate",";d_{z} [cm];Efficiency",70,0,700,2,-0.5,1.5,file_);
 add2DHistogram("h2DdxyVsPhiBTemplate",";d_{xy} [cm];#varphi bending;",100,-500,500,100,-500,500,file_);
 //Rate histos
 add2DHistogram("h2DRateTotTemplate","",402,1,202,402,1,202,file_);
 add2DHistogram("h2DRateVsEtaTemplate","",402,1,202,10,0.8,1.3,file_);
 add2DHistogram("h2DRateVsPtTemplate","",402,1,202,100,0,50,file_);
 //Misc histos
 add1DHistogram("h1DGenPtTemplate",";p_{T} [MeV/c];",60,0,120,file_);
 add1DHistogram("h1DGenEtaTemplate",";#eta;",300,-3.0,3.0,file_);
 add1DHistogram("h1DGenDxyTemplate",";d_{xy} [cm];",50,0,500,file_);
 add1DHistogram("h1DGenDzTemplate",";d_{z} [cm];",140,-700,700,file_);

 ///////////////////
 histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::finalizeHistograms(){

  AnalysisHistograms::finalizeHistograms();
  utilsL1RpcStyle()->cd();

  plotOMTFVsOther(OMTFHistograms::iPtCuts.at(3),"OMTFDispU");
  plotOMTFVsOther(OMTFHistograms::iPtCuts.at(3),"OMTFDispC");
  plotOMTFVsOther(OMTFHistograms::iPtCuts.at(3),"GMT");
  
  for(auto & anAlgo : algos){    
    plotEffPanel(anAlgo,"Pt");
    plotEffPanel(anAlgo,"HighPt");
    plotEffPanel(anAlgo,"dxy");
    plotEffPanel(anAlgo,"dz");
    plotEffVsVar(anAlgo,"EtaVx");
    plotEffVsVar(anAlgo,"PhiVx");
    plotEffVsVar(anAlgo,"dxy");
    plotEffVsVar(anAlgo,"dz");
    plotEffVsEta(anAlgo);
  }
   
  plotSingleHistogram("h2DOMTFPtRecVsPtGen");
  plotSingleHistogram("h2DOMTFDispUPtRecVsPtGen");
  plotSingleHistogram("h2DOMTFDispCPtRecVsPtGen");

  plotSingleHistogram("h2DOMTFdxyVsPhiB");
  plotSingleHistogram("h2DOMTFDispUdxyVsPhiB");
  plotSingleHistogram("h2DOMTFDispUdxyVsPhiBRefLayer0");
  plotSingleHistogram("h2DOMTFDispUdxyVsPhiBRefLayer2");
  
  plotSingleHistogram("h1DGenPt");
  plotSingleHistogram("h1DGenEta");
  plotSingleHistogram("h1DGenDxy");
  plotSingleHistogram("h1DGenDz");

  plotSingleHistogram("h1DGenEtaAll");
  plotSingleHistogram("h1DGenDxyAll");
  plotSingleHistogram("h1DGenDzAll");
  
  plotRate("Tot");
  plotRate("VsEta");
  plotRate("VsPt"); 
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void OMTFHistograms::plotEffPanel(const std::string & sysType, const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("EffVsPt_%s",sysType.c_str()),
			   TString::Format("EffVsPt_%s",sysType.c_str()),
			   460,500);
	c->SetGrid(0,1);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");  
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  
  TH1F hFrame("hFrame","",1,0,5000);    
  hFrame.SetStats(kFALSE);
  hFrame.SetMinimum(0.0001);
  hFrame.SetMaximum(1.04);
  hFrame.SetXTitle("gen. muon p_{T} [GeV/c]");
  hFrame.SetYTitle("Efficiency");
  hFrame.Draw();

  std::string hName(""); 
  for (int iCut=0; iCut<=3;++iCut){
    float ptCut = OMTFHistograms::ptBins[iPtCuts[iCut]];
    hName = "h2D"+sysType+varName+std::to_string((int)ptCut);
    TEfficiency *hEff = getEfficiency(hName);
    if(!hEff) return;

    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    hFrame.GetXaxis()->Set(1, h2D->GetXaxis()->GetXmin(), h2D->GetXaxis()->GetXmax());
    if(varName=="dxy") hFrame.SetXTitle("d_{xy} [cm]");
    else if(varName=="dz") hFrame.SetXTitle("d_{z} [cm]");

    hEff->SetMarkerStyle(21+iCut);
    hEff->SetMarkerColor(color[iCut]);
    hEff->Draw("same P");  
    
    TString nameCut = TString::Format("%d", (int)OMTFHistograms::ptBins[iPtCuts[iCut]])+" GeV/c";
    if (iCut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.Data(),"lp");
  }
  l.Draw();
  c->Print(TString::Format("fig_png/PanelVs%s_%s.png",varName.c_str(),sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotEffVsVar(const std::string & sysType,
		const std::string & varName){

  TCanvas* c = new TCanvas(TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   TString::Format("EffVs%s_%s",varName.c_str(),sysType.c_str()),
			   460,500);
	c->SetGrid(0,1);

  TLegend l(0.6513158,0.1673729,0.8903509,0.470339,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);
  
  TH1F hFrame("hFrame","",1,0,1);    
  hFrame.SetStats(kFALSE);
  hFrame.SetMinimum(0.0);
  hFrame.SetMaximum(1.04);
  hFrame.SetXTitle(varName.c_str());
  hFrame.SetYTitle("Efficiency");
  hFrame.Draw();

  std::string hName(""); 
  for (int iCut=2; iCut<=2;++iCut){
    float ptCut = OMTFHistograms::ptBins[iPtCuts[iCut]];
    hName = "h2D"+sysType+"Type0"+varName+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    hFrame.GetXaxis()->Set(1, h2D->GetXaxis()->GetXmin(), h2D->GetXaxis()->GetXmax());
    TEfficiency *hEff = getEfficiency(hName);
    hEff->SetMarkerStyle(21+iCut);
    hEff->SetMarkerColor(color[iCut]);
    hEff->Draw("same P");  
    
    std::string nameCut = std::to_string((int)ptCut)+" GeV/c";
    if (iCut==0) nameCut = "no p_{T} cut";
    l.AddEntry(hEff,nameCut.c_str(),"lp");
  }
  l.Draw();
  c->Print(TString::Format("fig_eps/EffVs%s_%s.eps",varName.c_str(), sysType.c_str()).Data());
  c->Print(TString::Format("fig_png/EffVs%s_%s.png",varName.c_str(), sysType.c_str()).Data());
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotEffVsEta(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("EffVsEta_%s",sysType.c_str()),
			   TString::Format("EffVsEta_%s",sysType.c_str()),
			   800,500);
	c->SetGrid(0,1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.35);
  
  TLegend l(0.6915995,0.5930233,0.7422325,0.8972868,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
 
  TH1F hFrame("hFrame","",1,0,1);    
  hFrame.SetStats(kFALSE);
  hFrame.SetMinimum(0.0);
  hFrame.SetMaximum(1.04);
  hFrame.GetYaxis()->SetTitleOffset(1.0);
  hFrame.SetXTitle("generated muon #eta");
  hFrame.SetYTitle("Efficiency");
  hFrame.Draw();

  int iCut = 14;
  std::string hName = "";
  for (int iType=0; iType<3;++iType){
    float ptCut = OMTFHistograms::ptBins[iCut];
    hName = "h2D"+sysType+"Type" + std::to_string(iType) + "EtaVx"+std::to_string((int)ptCut);
    TH2F* h2D = this->get2DHistogram(hName);
    if(!h2D) return;
    hFrame.GetXaxis()->Set(1, h2D->GetXaxis()->GetXmin(), h2D->GetXaxis()->GetXmax());
    
    TH1D *hNum = h2D->ProjectionX("hNum",2,2);
    TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);
    if(iType==2) hNum->Scale(50.0);
    hDenom->Add(hNum);
    TEfficiency *hEff = new TEfficiency(*hNum,*hDenom);
    hEff->SetMarkerStyle(21+iType);
    hEff->SetMarkerColor(color[iType]);
    hEff->Draw("same P");
    
    std::string nameCut = std::to_string((int)OMTFHistograms::ptBins[iCut])+" GeV/c";
    if (iType==0) nameCut = "p_{T}^{#mu}>p_{T}^{cut} + 20 GeV/c";
    if (iType==1) nameCut = "p_{T}^{cut}<p_{T}^{#mu}<#dot p_{T}^{cut} + 5 GeV/c";
    if (iType==2) nameCut = "p_{T}^{#mu}<10 GeV/c (#epsilon #times 50)";
    l.AddEntry(hEff,nameCut.c_str(),"lp");
  }
  TLine *aLine = new TLine(0,0,0,0);
  aLine->SetLineWidth(2);
  aLine->SetLineColor(2);
  aLine->DrawLine(0.83,0,0.83,1.0);
  aLine->DrawLine(-0.83,0,-0.83,1.0);
  aLine->DrawLine(1.24,0,1.24,1.0);
  aLine->DrawLine(-1.24,0,-1.24,1.0);

  l.SetHeader(TString::Format("p_{T}^{cut} = %d  GeV/c",(int)OMTFHistograms::ptBins[iCut]).Data());
  l.Draw();
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

  TLegend l(0.2,0.65,0.44,0.86,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetLogx(1);
  c->SetGrid(0,1);

  std::string hName = "h2D"+sysType+"Pt"+std::to_string((int)ptCut);
  TEfficiency *hEffOther = getEfficiency(hName);
  if(!hEffOther) return;
  hEffOther->SetMarkerStyle(23);
  hEffOther->SetMarkerColor(2);

  hName = "h2DOMTFPt"+std::to_string((int)ptCut);
  TEfficiency *hEffOMTF = getEfficiency(hName);
  hEffOMTF->SetMarkerStyle(21);
  hEffOMTF->SetMarkerColor(1);
  if(!hEffOMTF) return;
  
  TH1F hFrame("hFrame","",1,1,50);    
  hFrame.SetStats(kFALSE);
  hFrame.SetMinimum(0.0001);
  hFrame.SetMaximum(1.04);
  hFrame.SetXTitle("gen. muon p_{T} [GeV/c]");
  hFrame.SetYTitle("Efficiency");
  hFrame.Draw();
  
  hEffOMTF->Draw("same P");
  hEffOther->Draw("same P");
    
  std::string tmp = "p_{T} #geq ";
  if(int(ptCut*10)%10==5) tmp += "%1.1f GeV/c";
  else   tmp += "%1.0f GeV/c";
  l.AddEntry((TObject*)0, TString::Format(tmp.c_str(),ptCut).Data(), "");
  l.AddEntry((TObject*)0, "", "");
  tmp = sysType;
  if(sysType=="OMTFDispU") tmp = "unconstrained p_{T}";
  if(sysType=="OMTFDispC") tmp = "constrained p_{T}";
  l.AddEntry(hEffOther, tmp.c_str(), "lp");
  l.AddEntry(hEffOMTF, "CMSSW default", "lp");
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
TH1* OMTFHistograms::getRateHisto(std::string sysType,
				 std::string type){

  std::string hName = "h2D"+sysType+"Rate"+type;

  TH2F* h2D_original = (TH2F*)this->get2DHistogram(hName);
  if(!h2D_original) return 0;
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
  return Integrate(hRate);
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
void OMTFHistograms::plotRate(std::string type){
  
  TH1 *hRateVx = getRateHisto("Vx",type);
  TH1 *hRateOMTF = getRateHisto("OMTF",type);
  TH1 *hRateOMTFdisp = getRateHisto("OMTFDisp",type);
  TH1 *hRateOMTFdispC = getRateHisto("OMTFDispC",type);
  
  if(!hRateVx || !hRateOMTF || !hRateOMTFdispC || !hRateOMTFdisp) return;

  int iPtCut = OMTFHistograms::iPtCuts.at(3);
  float ptCut = OMTFHistograms::ptBins.at(iPtCut);

  hRateVx->SetLineWidth(3);
  hRateOMTF->SetLineWidth(3);
  hRateOMTFdispC->SetLineWidth(3);
  hRateOMTFdisp->SetLineWidth(3);

  hRateVx->SetLineColor(1);
  hRateOMTF->SetLineColor(4);
  hRateOMTFdispC->SetLineColor(2);
  hRateOMTFdisp->SetLineColor(kYellow+4);

  hRateOMTFdispC->SetLineStyle(2);
  hRateOMTFdisp->SetLineStyle(1);

  TCanvas* c = new TCanvas("cRate","Rate", 600, 650);
  c->SetLogy(1);
  c->SetGrid(1,1);

  TLegend *leg = new TLegend(0.43,0.8,0.8,0.95,NULL,"brNDC");
  leg->SetTextSize(0.04);
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);

  if(type.find("Tot")!=std::string::npos){
    hRateVx->GetXaxis()->SetRangeUser(2,50);
    hRateVx->SetMinimum(1E1);
    hRateVx->SetMaximum(2E5);
    
    c->Divide(2);
    TPad *pad1 = (TPad*)c->GetPad(1);
    TPad *pad2 = (TPad*)c->GetPad(2);
    pad1->SetPad(0.01,0.29,0.99,0.99);
    pad2->SetPad(0.01,0.01,0.99,0.29);

    pad1->SetTopMargin(0.0);
    pad1->SetBottomMargin(0.008);
    
    pad2->SetTopMargin(0.0); 
    pad2->SetBottomMargin(0.28);
    
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy();
    pad1->SetGrid(1,0);
    hRateVx->Draw();
    hRateOMTF->DrawCopy("same");
    hRateOMTFdispC->DrawCopy("same");
    hRateOMTFdisp->DrawCopy("same");
    
    std::cout<<"Rate OMTF CMSSW default @ "<<ptCut<<" GeV "<< hRateOMTF->GetBinContent(hRateOMTF->FindBin(ptCut+0.01))<<std::endl;
    std::cout<<"Rate OMTFdispC @ "<<ptCut<<" GeV "<< hRateOMTFdispC->GetBinContent(hRateOMTFdispC->FindBin(ptCut+0.01))<<std::endl;

    c->cd();
    pad2->Draw();
    pad2->cd();

    hRateOMTFdispC->GetXaxis()->SetRangeUser(2,50);
    hRateOMTFdispC->SetXTitle("p_{T}^{cut} [GeV/c]");
    hRateOMTFdispC->SetYTitle("new/default");
    hRateOMTFdispC->GetXaxis()->SetLabelSize(0.1);
    hRateOMTFdispC->GetXaxis()->SetTitleSize(0.1);
    hRateOMTFdispC->GetYaxis()->SetLabelSize(0.1);
    hRateOMTFdispC->GetYaxis()->SetTitleSize(0.1);    
    hRateOMTFdispC->GetYaxis()->SetTitleOffset(0.5);
    hRateOMTFdispC->Divide(hRateOMTF);
    hRateOMTFdisp->Divide(hRateOMTF);
    hRateOMTFdispC->SetMaximum(2.6);
    hRateOMTFdispC->SetMinimum(0.7);
    hRateOMTFdispC->DrawCopy();
    hRateOMTFdisp->DrawCopy("same");

    TLine *aLine = new TLine(0,0,0,0);
    aLine->SetLineWidth(2);
    aLine->SetLineColor(2);
    aLine->DrawLine(2,1.0, 50,1.0);
    c->cd();
  }  
  else if(type.find("VsEta")!=std::string::npos){
    c->SetLogy(0);
    hRateOMTFdispC->SetXTitle("muon #eta");
    double max = std::max(hRateOMTFdispC->GetMaximum(), hRateOMTF->GetMaximum());
    max = std::max(max, hRateOMTFdisp->GetMaximum());
    hRateOMTFdispC->SetMaximum(1.5*max);
    hRateOMTFdispC->Draw();
    hRateOMTFdisp->Draw("same");
    hRateOMTF->Draw("same");
    leg->SetHeader(TString::Format("p_{T}^{cut} = %d  GeV/c", int(ptCut)).Data());
  }
  else if(type=="VsPt"){
    c->SetLogy(1);
    hRateOMTFdispC->SetXTitle("p_{T}^{gen} [GeV/c]");
    hRateOMTFdispC->SetAxisRange(2,100);
    double max = std::max(hRateOMTFdispC->GetMaximum(), hRateOMTF->GetMaximum());
    max = std::max(max, hRateOMTFdisp->GetMaximum());
    hRateOMTFdispC->SetMaximum(10*max);
    hRateOMTFdispC->Draw();
    hRateOMTFdisp->Draw("same");
    hRateOMTF->Draw("same");
    leg->SetHeader(TString::Format("p_{T}^{cut} = %d  GeV/c",int(ptCut)).Data());
  }
 
  leg->AddEntry(hRateOMTF,"CMSSW default");
  leg->AddEntry(hRateOMTFdispC,"p_{T} constr.");
  leg->AddEntry(hRateOMTFdisp,"max(unconstr., constr.)");
  leg->Draw();

  c->Print(("fig_eps/Rate"+type+".eps").c_str());
  c->Print(("fig_png/Rate"+type+".png").c_str());  
}
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
TEfficiency * OMTFHistograms::getEfficiency(const std::string & hName){

  TH2F* h2D = this->get2DHistogram(hName);
  if(!h2D) return 0;
  
  TH1D *hNum = h2D->ProjectionX("hNum",2,2);
  TH1D *hDenom = h2D->ProjectionX("hDenom",1,1);  
  hDenom->Add(hNum);
  return new TEfficiency(*hNum, *hDenom);
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

  if(h2D && hName.find("PtRecVsPtGen")!=std::string::npos) {
    h2D->SetDirectory(myDirCopy);
    h2D->SetLineWidth(3);
    h2D->Scale(1.0/h2D->Integral());
    h2D->SetXTitle("p_{T}^{GEN}");
    h2D->SetYTitle("AA p_{T}^{REC}");
    h2D->GetYaxis()->SetTitleOffset(1.4);
    h2D->SetStats(kFALSE);
    h2D->GetYaxis()->SetRangeUser(0,150);
    h2D->GetXaxis()->SetRangeUser(0,100);
    h2D->Draw("candlex(000311)");
    TF1 fDiag("fDiag","x",0,100);
    fDiag.DrawCopy("same");
  }
  else if(h2D){
    h2D->SetXTitle("d_{xy} [cm]");
    h2D->SetYTitle("#varphi bending");
    h2D->SetZTitle("entries");
    h2D->Draw("col");
  }
  else if(h1D) {
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
  }
  if(hName.find("dxyVsPhiBRefLayer0")!=std::string::npos){
    TF1 fDiag("fDiag","-3.8/3.0*x",-300,300);
    fDiag.DrawCopy("same");
  }
  if(hName.find("dxyVsPhiBRefLayer2")!=std::string::npos){
    TF1 fDiag("fDiag","-3.0/3.0*x",-300,300);
    fDiag.DrawCopy("same");
  }
   c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
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
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////