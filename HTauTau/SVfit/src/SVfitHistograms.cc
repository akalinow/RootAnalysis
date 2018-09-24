#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "SVfitHistograms.h"
#include "AnalysisEnums.h"
#include "Tools.h"

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
SVfitHistograms::SVfitHistograms(TDirectory *myDir){
        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
SVfitHistograms::SVfitHistograms(TDirectory *myDir, const std::vector<std::string> & flavours, std::string channel){

        selectionFlavours_ = flavours;
        myChannel_ = channel;

        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
SVfitHistograms::~SVfitHistograms(){
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string SVfitHistograms::getTemplateName(const std::string& name){

        std::string templateName = "TemplateNotFound: "+name;
        if(name.find("hProf")!=std::string::npos && name.find("VsMag")!=std::string::npos) templateName = "hProfVsMagTemplate";
        else if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";
        else if(name.find("h1DDelta")!=std::string::npos) templateName = "h1DDeltaTemplate";
        else if(name.find("h1DCpuTime")!=std::string::npos) templateName = "h1DCpuTimeTemplate";
        return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::defineHistograms(){

        using namespace std;

        if(!histosInitialized_) {
                add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
                add1DHistogram("h1DMassTemplate",";mass [GeV/c^{2}]; Events",100,0,500,file_);
                add1DHistogram("h1DDeltaTemplate","",201,-5,5,file_);
                add1DHistogram("h1DCpuTimeTemplate","",200,0.0,0.01,file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

  std::cout<<"SVfitHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  std::vector<std::string> names = {"Pythia8"};

  for(unsigned int i=0;i<names.size();++i){
    std::string hNameSuffix = names[i];

    plotSingleHistogram("h1DMassVis"+hNameSuffix);
    plotSingleHistogram("h1DMassGen"+hNameSuffix);
    plotSingleHistogram("h1DMassFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DDeltaFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeFast"+hNameSuffix);
  }
  std::cout<<"SVfitHistograms::finalizeHistograms() END"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::plotSingleHistogram(std::string hName){

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
void SVfitHistograms::plotSingleHistogram2D(std::string hName){

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
