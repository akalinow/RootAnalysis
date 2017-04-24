#include <iostream>
#include <cmath>

#include "commonUtils.h"
#include "HZZHistograms.h"
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
HZZHistograms::HZZHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HZZHistograms::HZZHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HZZHistograms::HZZHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HZZHistograms::~HZZHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HZZHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";

  if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";

  return templateName;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HZZHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   add1DHistogram("h1DMassTemplate","",50,0,200,file_);
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HZZHistograms::finalizeHistograms(){

  plotAnyHistogram("h1DMass4Mu");
  plotAnyHistogram("h1DMassZ1");
  plotAnyHistogram("h1DMassZ2");

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HZZHistograms::plotAnyHistogram(const std::string & hName){

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
    h1D->SetYTitle("Events");
    h1D->SetXTitle("m [GeV/c^{2}]");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->SetLineColor(1);
    h1D->SetFillColor(1);
    h1D->SetMarkerStyle(20);
    h1D->Draw("HIST p");
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
