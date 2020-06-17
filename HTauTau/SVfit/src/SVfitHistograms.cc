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
#include "TProfile2D.h"

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
	else if(name.find("h2DDelta")!=std::string::npos) templateName = "h2DDeltaTemplate";
	else if(name.find("h3DDelta")!=std::string::npos) templateName = "h3DDeltaTemplate";
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
                add1DHistogram("h1DCpuTimeTemplate","",10,0.0,0.2,file_);
		add1DHistogram("h1DTauIDTemplate","",200,-1, 1,file_);
		///
		add2DHistogram("h2DDeltaTemplate","",10,0,100, 41,-10,10,file_);
		///
		add3DHistogram("h3DDeltaTemplate","",20,-0.5,0.5,30,-3,3,30, 20, 100, file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

  std::cout<<"SVfitHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  std::vector<std::string> names = {"Pythia8","ggHTT125", "DY0Jets", "DY0JetsMatchT"};

  for(unsigned int i=0;i<names.size();++i){
    std::string hNameSuffix = names[i];
    
    plotSingleHistogram("h1DMassVis"+hNameSuffix);
    plotSingleHistogram("h1DMassGen"+hNameSuffix);
    plotSingleHistogram("h1DMassFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DMassSVClassic"+hNameSuffix);
    plotSingleHistogram("h1DMassCA"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMCA"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMSVClassic"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeSVClassic"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMET_X_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMET_Y_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg1_E_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_E_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PX_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PY_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2_PZ_Res"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPhiFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaFastMTT"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtFastMTT"+hNameSuffix);    
    plotSingleHistogram("h1DDeltaPhiClassic"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaClassic"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtClassic"+hNameSuffix);

    plotSingleHistogram2D("h2DDeltaMET_X_Res_Vs_Mass"+hNameSuffix);
    
    plot3DProfile("h3DDeltaLeg2_E_Res_Vs_Eta_Vs_E"+hNameSuffix,"yz");
  }

  names = {"training", "DPFTau_2016_v0tauVSall", "DPFTau_2016_v1tauVSall", "deepTau2017v1tauVSjet", "MVArun2v1DBnewDMwLTraw2017v2"};
  for(unsigned int i=0;i<names.size();++i){
    std::string hNamePrefix = "h1DTauID_"+names[i];
    plotSingleHistogram(hNamePrefix+"ggHTT125");
    plotSingleHistogram(hNamePrefix+"QCD_MC");
  }
  plotROC("ggHTT125", "QCD_MC");
  
  std::cout<<"SVfitHistograms::finalizeHistograms() END"<<std::endl;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void SVfitHistograms::plotROC(std::string signal, std::string background){

  TH1D *h1 = plotROC("h1DTauID_MVArun2v1DBnewDMwLTraw2017v2", "ggHTT125", "QCD_MC");
  TH1D *h2 = plotROC("h1DTauID_DPFTau_2016_v1tauVSall", "ggHTT125", "QCD_MC");
  TH1D *h3 = plotROC("h1DTauID_deepTau2017v1tauVSjet", "ggHTT125", "QCD_MC");
  TH1D *h4 = plotROC("h1DTauID_training", "ggHTT125", "QCD_MC");

  if(!h1 || !h2 || !h3 || !h4) return;
  
  TCanvas* c = new TCanvas("ROC","AnyHistogram", 460,500);
  c->SetLogy();
  c->SetGrid();
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  TLegend l(0.15,0.68,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
  l.AddEntry(h1,"MVArun2v1DBnewDMwLTraw2017v2");
  l.AddEntry(h2,"DPFTau_2016_v1tauVSall");
  l.AddEntry(h3,"deepTau2017v1tauVSjet");
  l.AddEntry(h4,"training");    
  l.Draw();
  
  c->Print("fig_png/tauID_ROC.png");

  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1D* SVfitHistograms::plotROC(std::string hName, std::string signal, std::string background){

  std::string tmpName = hName+signal;
  TH1F* h1DSignal = get1DHistogram(tmpName);
  if(!h1DSignal) return 0;

  tmpName = hName+background;
  TH1F* h1DBackground = get1DHistogram(tmpName);
  if(!h1DBackground) return 0;

  TMVA::ROCCalc aROCCalculator(h1DSignal, h1DBackground);

  //aROCCalculator.Root() method is private in ROOT. Needs a hack in ROOT code.
  /*double effS, effSError;
  std::cout<<"Discriminator: "<<hName<<std::endl;
  effS = aROCCalculator.GetEffSForEffBof(1E-3, effSError);
  std::cout<<"60% S efficiency threshold: "<<aROCCalculator.Root(0.6)
  	   <<" S efficiency: "<<effS<<std::endl;
  effS = aROCCalculator.GetEffSForEffBof(2E-3, effSError);
  std::cout<<"70% S efficiency threshold: "<<aROCCalculator.Root(0.7)
  	   <<" S efficiency: "<<effS<<std::endl;
  */
  
  TH1D * h1DROC = aROCCalculator.GetROC();
  h1DROC->SetName((hName+"_ROC").c_str());
  h1DROC->SetDirectory(myDirCopy);
  double value = 0.0;
  for(int iBinX=0;iBinX<=h1DROC->GetNbinsX()+1;++iBinX){
    value = h1DROC->GetBinContent(iBinX);
    h1DROC->SetBinContent(iBinX, 1.0 - value);
  }
  
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
			   460,500);
  c->SetLogy();
  h1DROC->SetYTitle("#epsilon_{B}");
  h1DROC->SetAxisRange(0.15, 1.0, "X");
  h1DROC->SetMinimum(2E-4);
  h1DROC->SetMaximum(0.2);
  h1DROC->SetStats(kFALSE);
  h1DROC->SetLineWidth(2);
  h1DROC->Draw();
  tmpName = "fig_png/"+hName+"_"+signal+"_"+background+"_ROC.png";
  c->Print(tmpName.c_str());

  return (TH1D*)h1DROC->Clone();
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
void SVfitHistograms::plot3DProfile(std::string hName, std::string option){

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
