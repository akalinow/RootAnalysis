#include <iostream>
#include <cmath>

#include "commonUtils.h"
#include "TestHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TMath.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TestHistograms::TestHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TestHistograms::TestHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TestHistograms::TestHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TestHistograms::~TestHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool TestHistograms::fill1DHistogram(const std::string& name, float val, float weight){

  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DX")!=std::string::npos) hTemplateName = "h1DXTemplate";
    if(name.find("h1DY")!=std::string::npos) hTemplateName = "h1DYTemplate";
    std::cout<<"Adding histogram: "<<name<<" "<<file_<<" "<<file_->fullPath()<<std::endl;
    this->add1DHistogram(name,"",
			 this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmin(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmax(),
			 file_);
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TestHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   std::cout<<"Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;
   add1DHistogram("h1DXTemplate","",20,-200,200,file_);
   add1DHistogram("h1DYTemplate","",20,-200,200,file_);
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TestHistograms::finalizeHistograms(int nRuns, float weight){

  plotAnyHistogram("h1DXSignal");
  plotAnyHistogram("h1DYSignal");
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TestHistograms::plotAnyHistogram(const std::string & hName){
  
   TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TH1F* h1D = this->get1DHistogram(hName.c_str());

  if(h1D){
    h1D->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    if(hName.find("DeltaR")!=std::string::npos) c->SetLogy();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
