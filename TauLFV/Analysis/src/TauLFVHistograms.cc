#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "TauLFVHistograms.h"

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
TauLFVHistograms::TauLFVHistograms(TDirectory *myDir){

        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TauLFVHistograms::~TauLFVHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string TauLFVHistograms::getTemplateName(const std::string& name){

        std::string templateName = "TemplateNotFound: "+name;
        if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";

        return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TauLFVHistograms::defineHistograms(){

        using namespace std;

        if(!histosInitialized_) {
                add1DHistogram("h1DMassTemplate",";mass [GeV/c^{2}]; Events",20,0,5,file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TauLFVHistograms::finalizeHistograms(){

        std::cout<<"TauLFVHistograms::finalizeHistograms() START"<<std::endl;

        AnalysisHistograms::finalizeHistograms();

	plotSingleHistogram("h1DMass3Mu_DsToTau");

        std::cout<<"TauLFVHistograms::finalizeHistograms() END"<<std::endl;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void TauLFVHistograms::plotSingleHistogram(std::string hName){

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
    h2D->Scale(1.0/h2D->Integral(0,h2D->GetNbinsX()+1));
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
