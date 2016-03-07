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
bool HWvsEMULHistograms::fill1DHistogram(const std::string& name, float val1, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill1DHistogram(name,val1,weight)){
	  if(name.find("iProcessor")!=std::string::npos) hTemplateName = "h1DiProcessorTemplate";
	  if(name.find("Phi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
	  if(name.find("Eta")!=std::string::npos) hTemplateName = "h1DEtaTemplate";
	  if(name.find("Pt")!=std::string::npos) hTemplateName = "h1DPtTemplate";
	  if(name.find("Quality")!=std::string::npos) hTemplateName = "h1DQualityTemplate";
	  if(name.find("Delta")!=std::string::npos) hTemplateName = "h1DDeltaTemplate";
	  this->add1DHistogram(name,"",
			       this->get1DHistogram(hTemplateName)->GetNbinsX(),
			       this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
			       this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
			       file_);
	  this->get1DHistogram(name,true)->SetDirectory(this->get1DHistogram(hTemplateName,true)->GetDirectory());
	  return AnalysisHistograms::fill1DHistogram(name,val1,weight);
	}
	return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool HWvsEMULHistograms::fill2DHistogram(const std::string& name, float val1, float val2, float weight){

	std::string hTemplateName = "";
	if(!AnalysisHistograms::fill2DHistogram(name,val1,val2,weight)){
	  
		if(name.find("Pt")!=std::string::npos) hTemplateName = "h2DPtTemplate";
	
		this->add2DHistogram(name,"",
				this->get2DHistogram(hTemplateName)->GetNbinsX(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
				this->get2DHistogram(hTemplateName)->GetNbinsY(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmin(),
				this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmax(),
				file_);
		this->get2DHistogram(name,true)->SetDirectory(this->get2DHistogram(hTemplateName,true)->GetDirectory());
		return AnalysisHistograms::fill2DHistogram(name,val1,val2,weight);
	}
	return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::defineHistograms(){

  using namespace std;
  
  if(!histosInitialized_){
    
    //Make template histos
    add1DHistogram("h1DiProcessorTemplate","",11,-5.5,5.5,file_);
    add1DHistogram("h1DQualityTemplate","",6,-0.5,5.5,file_);
    add1DHistogram("h1DPtTemplate","",21,-0.5,20.5,file_);
    add1DHistogram("h1DPhiTemplate","",40,-3,0,file_);
    add1DHistogram("h1DEtaTemplate","",30,-2.2,2.2,file_);
    add1DHistogram("h1DDeltaTemplate","",201,-10,10,file_);
    
    histosInitialized_ = true;
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HWvsEMULHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  plotHWvsEMUL("iProcessor");
  plotHWvsEMUL("Quality");
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

  TLegend l(0.6,0.65,0.8,0.8,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TH1F* hHW = this->get1DHistogram("h1D"+type+"HW");
  TH1F* hEMUL = this->get1DHistogram("h1D"+type+"EMUL");

  if(type=="iProcessor") hHW->SetXTitle("iProcessor * sgn(iEta)");
  hHW->SetYTitle("Candidate count");
  hHW->GetYaxis()->SetTitleOffset(1.5);
  hHW->SetStats(kFALSE);

  hHW->SetLineWidth(3);
  hHW->SetLineColor(1);
  hHW->SetLineStyle(1);

  hEMUL->SetLineWidth(3);
  hEMUL->SetLineColor(2);
  hEMUL->SetLineStyle(2);

  if(hEMUL->GetMaximum()>hHW->GetMaximum()) hHW->SetMaximum(1.05*hEMUL->GetMaximum());

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

  TLegend l(0.6,0.55,0.8,0.7,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  c->SetGrid(0,1);

  TH1F* hHWPositive = this->get1DHistogram("h1DPhiHWPositive");
  TH1F* hHWNegative = this->get1DHistogram("h1DPhiHWNegative");

  hHWPositive->SetXTitle("#varphi [radians]");
  hHWPositive->SetYTitle("Candidate count");
  hHWPositive->GetYaxis()->SetTitleOffset(1.5);
  hHWPositive->SetStats(kFALSE);

  hHWPositive->SetLineWidth(3);
  hHWPositive->SetLineColor(1);
  hHWPositive->SetLineStyle(1);

  hHWNegative->SetLineWidth(3);
  hHWNegative->SetLineColor(2);
  hHWNegative->SetLineStyle(2);

  if(hHWNegative->GetMaximum()>hHWPositive->GetMaximum()) hHWPositive->SetMaximum(1.05*hHWNegative->GetMaximum());

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

  TH1F* hHW = this->get1DHistogram(("h1D"+type).c_str());

  hHW->SetXTitle(type.c_str());
  hHW->SetYTitle("Candidate count");
  hHW->GetYaxis()->SetTitleOffset(1.5);
  hHW->SetStats(kFALSE);

  hHW->SetLineWidth(3);
  hHW->SetLineColor(1);
  hHW->SetLineStyle(1);

  hHW->Draw();

  c->Print(TString::Format("fig_png/%s.png",type.c_str()).Data());
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
  if(type=="Pt") hDelta->GetXaxis()->SetRangeUser(-6,2);
  if(type=="Phi") hDelta->GetXaxis()->SetRangeUser(-1,1);
  if(type=="Eta") hDelta->GetXaxis()->SetRangeUser(-2,2);
  
  hDelta->SetStats(kFALSE);

  hDelta->SetLineWidth(3);
  hDelta->SetLineColor(1);
  hDelta->SetLineStyle(1);

  hDelta->Draw();

  c->Print(TString::Format("fig_png/Delta_%s.png",type.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
