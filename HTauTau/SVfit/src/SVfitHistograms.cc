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
#include "TMVA/ROCCalc.h"


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
	else if(name.find("h1DTauID")!=std::string::npos) templateName = "h1DTauIDTemplate";
	else if(name.find("h2DDelta")!=std::string::npos) templateName = "h2DDeltaTemplate";	
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
		add1DHistogram("h1DTauIDTemplate","",200,-1, 1,file_);

		add2DHistogram("h2DDeltaTemplate","",10,0,100, 41,-10,10,file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

  std::cout<<"SVfitHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  std::vector<std::string> names = {"Pythia8","ggHTT125", "DY0Jets"};

  for(unsigned int i=0;i<names.size();++i){
    std::string hNameSuffix = names[i];

    plotSingleHistogram("h1DMassVis"+hNameSuffix);
    plotSingleHistogram("h1DMassGen"+hNameSuffix);
    plotSingleHistogram("h1DMassFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DMassCA"+hNameSuffix);
    plotSingleHistogram("h1DDeltaFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DDeltaCA"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMET_X_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMET_Y_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg1_E_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_E_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PX_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PY_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PZ_Res"+hNameSuffix);

    plotSingleHistogram2D("h2DDeltaMET_X_Res_Vs_Mass"+hNameSuffix);
  }

  names = {"DPFTau_2016_v1tauVSall", "deepTau2017v1tauVSjet", "MVArun2v1DBoldDMwLTraw2017v2"};
  for(unsigned int i=0;i<names.size();++i){
    std::string hNamePrefix = "h1DTauID_"+names[i];
    plotSingleHistogram(hNamePrefix+"ggHTT125");
    plotSingleHistogram(hNamePrefix+"QCD_MC");
  }
  plotROC("h1DTauID_DPFTau", "ggHTT125", "QCD_MC");
  
  
  std::cout<<"SVfitHistograms::finalizeHistograms() END"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::plotROC(std::string hName, std::string signal, std::string background){

  std::string tmpName = hName+signal;
  TH1F* h1DSignal = get1DHistogram(tmpName);
  if(!h1DSignal) return;

  hName = hName+background;
  TH1F* h1DBackground = get1DHistogram(tmpName);
  if(!h1DBackground) return;

  TMVA::ROCCalc aROCCalculator(h1DSignal, h1DBackground);
  TH1D * h1DROC = aROCCalculator.GetROC();

  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   460,500);
  h1DROC->Draw();
  tmpName = "fig_png/"+hName+"_"+signal+"_"+background+".png";
  c->Print(tmpName.c_str());
  
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
