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
        else if(name.find("h1DDelta")!=std::string::npos) templateName = "h1DDeltaTemplate";
        else if(name.find("h1DLLH")!=std::string::npos) templateName = "h1DLLHTemplate";
        else if(name.find("h1DCpuTime")!=std::string::npos) templateName = "h1DCpuTimeTemplate";
	else if(name.find("h1DCosGJ")!=std::string::npos) templateName = "h1DCosGJTemplate";

        else if(name.find("h2DFlightPathVsDeltaR")!=std::string::npos) templateName = "h2DFlightPathVsDeltaRTemplate";
        else if(name.find("h2DDelta")!=std::string::npos) templateName = "h2DDeltaTemplate";
	else if(name.find("h2DLikelihood")!=std::string::npos) templateName = "h2DLikelihoodTemplate";
	else if(name.find("h2DMVisRatio")!=std::string::npos) templateName = "h2DMVisRatioTemplate";
	else if(name.find("h2DMHvs")!=std::string::npos) templateName = "h2DMHvsTemplate";

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
                add1DHistogram("h1DDeltaTemplate","",201,-5,5,file_);
                add1DHistogram("h1DLLHTemplate","",200,-0.1,3,file_);
                add1DHistogram("h1DCpuTimeTemplate","",200,0.0,0.5,file_);
		add1DHistogram("h1DCosGJTemplate","",200,-0.999,1.001,file_);

                add2DHistogram("h2DFlightPathVsDeltaRTemplate",";flight path [cm]; #Delta R",20,0,4, 50,0,0.05,file_);
		add2DHistogram("h2DDeltaTemplate","",101,-1,1, 101,0,3, file_);
		add2DHistogram("h2DLikelihoodTemplate","",200,-1,1, 200,-1,1, file_);
		add2DHistogram("h2DMVisRatioTemplate","",100,0, 200, 60, 0,2, file_);
		add2DHistogram("h2DMHvsTemplate","",101,0,300, 20,0,200, file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void svfitHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

  std::cout<<"svfitHistograms::finalizeHistograms() START"<<std::endl;

  AnalysisHistograms::finalizeHistograms();

  std::vector<std::string> names = {"ggHTT125", "ggHTT140", "ggHTT200",
                                    "ggHTT250", "ggHTT500", "ggHTT1000",
                                    "DYAllJetsMatchT","DYAllJetsMatchL","WAllJets"};

  plotSingleHistogram2D("h2DDelta_1");
  plotSingleHistogram2D("h2DDelta_2");
  plotSingleHistogram2D("h2DMVisRatio");

  plotSingleHistogram("h1DDeltaCorrectSolution1");
  plotSingleHistogram("h1DDeltaCorrectSolution2");

  plotSingleHistogram("h1DDeltaMETRatioS1");
  plotSingleHistogram("h1DDeltaMETRatioS2");

  for(unsigned int i=0;i<names.size();++i){
    std::string hNameSuffix = names[i];

    plotSingleHistogram("h1DMassMET"+hNameSuffix);
    plotSingleHistogram("h1DMassSVCA"+hNameSuffix);
    plotSingleHistogram("h1DMassSVClassic"+hNameSuffix);
    plotSingleHistogram("h1DMassSVFast"+hNameSuffix);
    plotSingleHistogram("h1DMassSVSuperFast"+hNameSuffix);
    
    plotSingleHistogram("h1DMassSV3ProngFast"+hNameSuffix);
    plotSingleHistogram("h1DMassSV3ProngFastE"+hNameSuffix);
    
    plotSingleHistogram("h1DMassSVStandalone"+hNameSuffix);
    plotSingleHistogram("h1DMassSV3ProngStandalone"+hNameSuffix);

    plotSingleHistogram("h1DDeltaLeg1"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2"+hNameSuffix);

    plotSingleHistogram("h1DDeltaMassResolutionStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMassResolutionFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMassResolutionFastE"+hNameSuffix);

    plotSingleHistogram("h1DDeltaLeg1E"+hNameSuffix);
    plotSingleHistogram("h1DDeltaLeg2E"+hNameSuffix);

    plotSingleHistogram("h1DMassVis"+hNameSuffix);
    plotSingleHistogram("h1DMassGen"+hNameSuffix);

    plotSingleHistogram("h1DCpuTimeFast"+hNameSuffix);
    plotSingleHistogram("h1DCpuTimeClassic"+hNameSuffix);

    plotSingleHistogram("h1DDeltaPhiFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPxFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPyFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPzFast"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEFast"+hNameSuffix);

    plotSingleHistogram("h1DDeltaM"+hNameSuffix);
    plotSingleHistogram2D("h2DDeltaM"+hNameSuffix);

    plotSingleHistogram("h1DDeltaPhiStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEtaStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPtStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPxStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPyStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPzStandalone"+hNameSuffix);
    plotSingleHistogram("h1DDeltaEStandalone"+hNameSuffix);

    plotSingleHistogram("h1DDeltaPtClassic"+hNameSuffix);

    plotSingleHistogram("h1DDeltaMVisMSV"+hNameSuffix);

    plotSingleHistogram("h1DDeltaSV_X"+hNameSuffix);
    plotSingleHistogram("h1DDeltaSV_Y"+hNameSuffix);
    plotSingleHistogram("h1DDeltaSV_Z"+hNameSuffix);

    plotSingleHistogram("h1DDelta_FlightPathLeg1_Gen"+hNameSuffix);
    plotSingleHistogram("h1DDelta_FlightPathLeg1_Gen_CMS"+hNameSuffix);
    plotSingleHistogram("h1DDelta_FlightPathLeg2_Gen"+hNameSuffix);
    plotSingleHistogram("h1DDelta_FlightPathLeg2_Reco"+hNameSuffix);
    
    plotSingleHistogram("h1DCosGJReco"+hNameSuffix);
    plotSingleHistogram("h1DCosGJGen"+hNameSuffix);
    plotSingleHistogram("h1DDeltaReco"+hNameSuffix);
    plotSingleHistogram("h1DDeltaGen"+hNameSuffix);
    plotSingleHistogram("h1DDeltaDelta"+hNameSuffix);
    
    plotSingleHistogram("h1DDeltaCosGJ"+hNameSuffix);
    plotSingleHistogram("h1DDeltaSinGJ"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMA1"+hNameSuffix);
    plotSingleHistogram("h1DDeltaPA1"+hNameSuffix);
    
    plotSingleHistogram("h1DDeltaSolution1"+hNameSuffix);
    plotSingleHistogram("h1DDeltaSolution2"+hNameSuffix);
    plotSingleHistogram("h1DDeltaDelta"+hNameSuffix);
    plotSingleHistogram("h1DDeltaRatioGJMax"+hNameSuffix);

    plotSingleHistogram("h1DDeltaELeg1"+hNameSuffix);
    plotSingleHistogram("h1DDeltaELeg2"+hNameSuffix);
    plotSingleHistogram("h1DDeltaMET"+hNameSuffix);
    
    plotSingleHistogram2D("h2DDeltaELeg1"+hNameSuffix);
    plotSingleHistogram2D("h2DDeltaELeg2"+hNameSuffix);

    plotSingleHistogram2D("h2DDeltaM"+hNameSuffix);
    plotSingleHistogram2D("h2DLikelihoodMapNorm"+hNameSuffix);
    plotSingleHistogram2D("h2DLikelihoodMap"+hNameSuffix);
    plotSingleHistogram2D("h2DLikelihoodMapE"+hNameSuffix);
    plotSingleHistogram2D("h2DMVisRatio"+hNameSuffix);
    plotSingleHistogram2D("h2DMHvsMVis"+hNameSuffix);
    
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
                h2D->Draw("colz");
                c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
