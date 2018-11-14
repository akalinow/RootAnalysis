#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "HSCPHistograms.h"

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
#include "TProfile2D.h"


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HSCPHistograms::HSCPHistograms(TDirectory *myDir){
        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HSCPHistograms::~HSCPHistograms(){
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HSCPHistograms::getTemplateName(const std::string& name){

        std::string templateName = "TemplateNotFound: "+name;
        if(name.find("h1DPt")!=std::string::npos) templateName = "h1DPtTemplate";

        return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HSCPHistograms::defineHistograms(){

  using namespace std;

  if(!histosInitialized_) {
    add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
    add1DHistogram("h1DPtTemplate",";p_{T} [GeV/c]; Events",100,0,500,file_);
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HSCPHistograms::finalizeHistograms(){

  std::cout<<"HSCPHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  plotSingleHistogram("h1DPt1");
  plotSingleHistogram("h1DPt2");

  std::cout<<"HSCPHistograms::finalizeHistograms() END"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HSCPHistograms::plotSingleHistogram(std::string hName){

        TH1F* h1D = get1DHistogram(hName);
        if(!h1D) return;

        h1D->SetDirectory(myDirCopy);

        TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);

        TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
        l.SetTextSize(0.05);
        l.SetFillStyle(4000);
        l.SetBorderSize(0);
        l.SetFillColor(10);

        if(h1D) {
                h1D->Print();
                h1D->SetLineWidth(3);
                h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
                h1D->SetYTitle("Events");
                h1D->GetYaxis()->SetTitleOffset(1.4);
                //h1D->SetStats(kFALSE);
                h1D->Draw();
                c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HSCPHistograms::plotSingleHistogram2D(std::string hName){

        TH2F* h2D = get2DHistogram(hName);
        if(!h2D) return;

        h2D->SetDirectory(myDirCopy);

        TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);

        TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
        l.SetTextSize(0.05);
        l.SetFillStyle(4000);
        l.SetBorderSize(0);
        l.SetFillColor(10);

        if(h2D) {
                h2D->SetLineWidth(3);
                h2D->Scale(1.0/h2D->Integral(0,h2D->GetNbinsX()+1, 0,h2D->GetNbinsY()+1));
                h2D->SetZTitle("Events");
                h2D->GetYaxis()->SetTitleOffset(1.4);
                h2D->SetStats(kFALSE);
                h2D->Draw("colz");
                c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HSCPHistograms::plot3DProfile(std::string hName, std::string option){

        TH3F* h3D = get3DHistogram(hName);
        if(!h3D) return;

        h3D->SetDirectory(myDirCopy);

        TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);
	
	TProfile2D *h2DProf = h3D->Project3DProfile(option.c_str());
	h2DProf->SetStats(kFALSE);
	h2DProf->Draw("colz");
	c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
