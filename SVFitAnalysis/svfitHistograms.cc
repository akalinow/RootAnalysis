#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "svfitHistograms.h"
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
svfitHistograms::svfitHistograms(TDirectory *myDir){
        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
svfitHistograms::svfitHistograms(TDirectory *myDir, const std::vector<std::string> & flavours, std::string channel){

        selectionFlavours_ = flavours;
        myChannel_ = channel;

        AnalysisHistograms::init(myDir);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
svfitHistograms::~svfitHistograms(){
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string svfitHistograms::getTemplateName(const std::string& name){

        std::string templateName = "TemplateNotFound: "+name;
        if(name.find("hProf")!=std::string::npos && name.find("VsMag")!=std::string::npos) templateName = "hProfVsMagTemplate";
        else if(name.find("h1DFlightPath")!=std::string::npos) templateName = "h1DFlightPathTemplate";
        else if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";
        else if(name.find("h1DDeltaR")!=std::string::npos) templateName = "h1DDeltaRTemplate";
        else if(name.find("h1DDelta")!=std::string::npos) templateName = "h1DDeltaTemplate";
        else if(name.find("h1DLLH")!=std::string::npos) templateName = "h1DLLHTemplate";
        else if(name.find("h1DCpuTime")!=std::string::npos) templateName = "h1DCpuTimeTemplate";

        else if(name.find("h2DFlightPathVsDeltaR")!=std::string::npos) templateName = "h2DFlightPathVsDeltaRTemplate";
        else if(name.find("h2DDelta")!=std::string::npos) templateName = "h2DDeltaTemplate";

        return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void svfitHistograms::defineHistograms(){

        using namespace std;

        if(!histosInitialized_) {
                add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
                add1DHistogram("h1DMassTemplate",";mass [GeV/c^{2}]; Events",100,0,500,file_);
                add1DHistogram("h1DFlightPathTemplate",";flight path [cm]; Events",100,0,0.1,file_);
                add1DHistogram("h1DDeltaRTemplate","",21,-6.0,6.0,file_);
                add1DHistogram("h1DDeltaTemplate","",50,-2.0,3.0,file_);

                add1DHistogram("h1DLLHTemplate","",200,-0.1,3,file_);
                add1DHistogram("h1DCpuTimeTemplate","",25,0.2,0.4,file_);

                add2DHistogram("h2DFlightPathVsDeltaRTemplate",";flight path [cm]; #Delta R",20,0,4, 50,0,0.05,file_);
                add2DHistogram("h2DDeltaTemplate","",20,-3.0,3.0, 20, -1.0, 1.0, file_);

        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void svfitHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

  std::cout<<"svfitHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  std::vector<std::string> names = {"ggHTT125", "ggHTT140", "ggHTT200",
                                    "ggHTT250", "ggHTT500", "ggHTT1000",
                                    "DYAllJetsMatchT","WAllJets"};

  plotSingleHistogram("h1DDeltaR_1");
  plotSingleHistogram("h1DDeltaR_2");
  plotSingleHistogram("h1DDeltaR_3");
  plotSingleHistogram("h1DDeltaR_4");
  plotSingleHistogram("h1DLLH_1");
  plotSingleHistogram("h1DLLH_2");
  plotSingleHistogram("h1DLLH_3");
  plotSingleHistogram("h1DLLH_4");

  plotSingleHistogram2D("h2DDelta_1");
  plotSingleHistogram2D("h2DDelta_2");

  for(unsigned int i=0;i<names.size();++i){
    std::string hNameSuffix = names[i];

    plotSingleHistogram("h1DMassSVCA"+hNameSuffix);
    plotSingleHistogram("h1DMassSVClassic"+hNameSuffix);
    plotSingleHistogram("h1DMassSVFast"+hNameSuffix);
    plotSingleHistogram("h1DMassSVStandalone"+hNameSuffix);

    plotSingleHistogram("h1DCpuTimeFast"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeClassic"+hNameSuffix);

    plotSingleHistogram("h1DDeltaPhiFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtFast"+hNameSuffix);

    plotSingleHistogram("h1DDeltaPhiStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtClassic"+hNameSuffix);
    /*
    plotSingleHistogram("h1DFlightPathRec"+hNameSuffix);
    plotSingleHistogram("h1DFlightPathPCARec"+hNameSuffix);
    plotSingleHistogram("h1DFlightPathPCARecLeg1"+hNameSuffix);
    plotSingleHistogram("h1DFlightPathGen"+hNameSuffix);
    plotSingleHistogram("h1DFlightPathPCAGen"+hNameSuffix);

    plotSingleHistogram2D("h2DFlightPathVsDeltaRGen"+hNameSuffix);
    */
  }
  std::cout<<"svfitHistograms::finalizeHistograms() END"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void svfitHistograms::plotSingleHistogram(std::string hName){

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
void svfitHistograms::plotSingleHistogram2D(std::string hName){

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
                h2D->Draw();
                c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
