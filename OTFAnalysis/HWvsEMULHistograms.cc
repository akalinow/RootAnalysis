#include <iostream>
#include <cmath>

#include "HWvsEMULHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TRandom3.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HWvsEMULHistograms::HWvsEMULHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HWvsEMULHistograms::HWvsEMULHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HWvsEMULHistograms::HWvsEMULHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HWvsEMULHistograms::~HWvsEMULHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HWvsEMULHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";

  if(name.find("iProcessor")!=std::string::npos) templateName = "h1DiProcessorTemplate";
  if(name.find("Phi")!=std::string::npos) templateName = "h1DPhiTemplate";
  if(name.find("Eta")!=std::string::npos) templateName = "h1DEtaTemplate";
  if(name.find("1DPt")!=std::string::npos) templateName = "h1DPtTemplate";
  if(name.find("Quality")!=std::string::npos) templateName = "h1DQualityTemplate";

  if(name.find("DeltaPt")!=std::string::npos) templateName = "h1DDeltaPtTemplate";
  else if(name.find("Delta")!=std::string::npos) templateName = "h1DDeltaTemplate";
  
  if(name.find("2DPt")!=std::string::npos) templateName = "h2DPtTemplate";
  
  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::defineHistograms(){

  using namespace std;
  
  if(!histosInitialized_){
    
    //Make template histos
    add1DHistogram("h1DiProcessorTemplate","",12,-6,6,file_);
    add1DHistogram("h1DQualityTemplate","",13,-0.5,12.5,file_);
    add1DHistogram("h1DPtTemplate","",251,-0.5,250.5,file_);
    add1DHistogram("h1DPhiTemplate","",40,-3.141,3.141,file_);
    add1DHistogram("h1DEtaTemplate","",481,-240.5,240.5,file_);
    add1DHistogram("h1DDeltaTemplate","",81,-20.5,20.5,file_);
    add1DHistogram("h1DDeltaPtTemplate","",500,-250,250,file_);
    
    histosInitialized_ = true;
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  plotHWvsEMUL("iProcessor");
  plotHWvsEMUL("Quality");
  plotHWvsEMUL("Eta");
  plotHWvsEMUL("PtER");
  
  plotPhi();
  plotVariable("PtHW");
  plotVariable("EtaHW");
  plotVariable("QualityHW");
  plotDelta("Pt");
  plotDelta("Phi");
  plotDelta("Eta");
  plotDelta("Charge");
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::finalizeHistograms(){
  finalizeHistograms(0,1.0);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::plotHWvsEMUL(std::string type){

  TCanvas* c = new TCanvas("Counts","",460,500);

  TLegend l(0.6,0.5,0.8,0.7,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);
  c->SetLeftMargin(0.15);

  TH1F* hHW = this->get1DHistogram("h1D"+type+"HW");
  TH1F* hEMUL = this->get1DHistogram("h1D"+type+"EMUL");

  if(type=="iProcessor") hHW->SetXTitle("iProcessor * sgn(iEta)");
  if(type=="Eta"){
    hHW->GetXaxis()->SetRangeUser(70,115);
    //hHW->GetXaxis()->SetRangeUser(100,115);
    hHW->SetXTitle("iEta");
  }
  if(type=="PtER"){
    hHW->SetXTitle("p_{T} [GeV]");
    hHW->GetXaxis()->SetRangeUser(0,20);
  }
  if(type=="Quality") hHW->SetXTitle("Quality");
  
  hHW->SetYTitle("Candidate count");
  hHW->GetYaxis()->SetTitleOffset(1.7);
  hHW->SetStats(kFALSE);

  hHW->SetLineWidth(3);
  hHW->SetLineColor(1);
  hHW->SetLineStyle(1);

  hEMUL->SetLineWidth(3);
  hEMUL->SetLineColor(2);
  hEMUL->SetLineStyle(2);

  if(hEMUL->GetMaximum()>hHW->GetMaximum()) hHW->SetMaximum(1.05*hEMUL->GetMaximum());
  if(hEMUL->GetMinimum()<hHW->GetMinimum()) hHW->SetMinimum(0.95*hEMUL->GetMinimum());

  hHW->Draw();
  hEMUL->Draw("same");

  l.AddEntry(hHW,"Hardware");
  l.AddEntry(hEMUL,"Emulator");
  l.Draw();

  c->Print(TString::Format("fig_png/HWvsEMUL_%s.png",type.c_str()).Data());

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::plotPhi(){

  TCanvas* c = new TCanvas("Counts","",460,500);

  TLegend l(0.5,0.5,0.7,0.6,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TH1F* hHWPositive = this->get1DHistogram("h1DPhiHWPositive");
  TH1F* hHWNegative = this->get1DHistogram("h1DPhiHWNegative");

  hHWPositive->SetXTitle("#varphi [radians]");
  hHWPositive->SetYTitle("Candidate count");
  hHWPositive->GetYaxis()->SetTitleOffset(1.55);
  hHWPositive->SetStats(kFALSE);

  hHWPositive->SetLineWidth(3);
  hHWPositive->SetLineColor(1);
  hHWPositive->SetLineStyle(1);

  hHWNegative->SetLineWidth(3);
  hHWNegative->SetLineColor(2);
  hHWNegative->SetLineStyle(2);

  hHWPositive->Rebin(4);
  hHWNegative->Rebin(4);

  if(hHWNegative->GetMaximum()>hHWPositive->GetMaximum()) hHWPositive->SetMaximum(1.05*hHWNegative->GetMaximum());
  hHWPositive->SetMinimum(0);

  hHWPositive->Draw();
  hHWNegative->Draw("same");

  l.AddEntry(hHWPositive,"Positive eta");
  l.AddEntry(hHWNegative,"Negative eta");
  l.Draw();

  c->Print("fig_png/PhiHW.png");

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::plotVariable(std::string type){

  TCanvas* c = new TCanvas("Counts","",460,500);

  TLegend l(0.6,0.55,0.8,0.7,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);
  c->SetLeftMargin(0.15);
  

  TH1F* hHW = this->get1DHistogram(("h1D"+type).c_str());

  hHW->SetXTitle(type.c_str());
  hHW->SetYTitle("Candidate count");
  hHW->GetYaxis()->SetTitleOffset(1.7);
  hHW->SetStats(kFALSE);

  hHW->SetLineWidth(3);
  hHW->SetLineColor(1);
  hHW->SetLineStyle(1);

  if(type=="EtaHW") hHW->GetXaxis()->SetRangeUser(-116,116);
  if(type=="PtHW") hHW->GetXaxis()->SetRangeUser(0,20);
  hHW->Draw();
  c->Print(TString::Format("fig_png/%s.png",type.c_str()).Data());

  if(type.find("PtHW")!=std::string::npos){
    hHW->GetXaxis()->SetRangeUser(0,250);
    hHW->Draw();
    c->SetLogy();
    c->Print(TString::Format("fig_png/%s_Range.png",type.c_str()).Data());
  }
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::plotDelta(std::string type){

  TCanvas* c = new TCanvas("Delta","",460,500);

  TLegend l(0.6,0.55,0.8,0.7,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TH1F* hDelta = this->get1DHistogram("h1DDelta"+type);

  hDelta->SetXTitle(("#Delta(HW-EMUL) "+type).c_str());
  hDelta->SetYTitle("Candidate count");
  hDelta->GetYaxis()->SetTitleOffset(1.6);
  if(type=="Pt"){
    c->SetLogy();
  }
  if(type=="Phi") hDelta->GetXaxis()->SetRangeUser(-5,5);
  if(type=="Eta") hDelta->GetXaxis()->SetRangeUser(-20,20);
  if(type=="Charge") hDelta->GetXaxis()->SetRangeUser(-3,3);
  
  hDelta->SetStats(kFALSE);

  hDelta->SetLineWidth(3);
  hDelta->SetLineColor(1);
  hDelta->SetLineStyle(1);

  if(type=="Pt") hDelta->Draw("LH");
  else hDelta->Draw();

  c->Print(TString::Format("fig_png/Delta_%s.png",type.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
