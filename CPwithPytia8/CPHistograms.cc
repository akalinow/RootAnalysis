#include <iostream>
#include <cmath>

#include "CPHistograms.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TMath.h"

CPHistograms::CPHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::CPHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::CPHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
CPHistograms::~CPHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool CPHistograms::fill1DHistogram(const std::string& name, float val, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DPhi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
    if(name.find("h1DRho")!=std::string::npos) hTemplateName = "h1DRhoTemplate";
    std::cout<<"Adding histogram: "<<name<<std::endl;
    this->add1DHistogram(name,"",
			 this->get1DHistogram(hTemplateName)->GetNbinsX(),
			 this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
			 this->get1DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
			 file_);
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   add1DHistogram("h1DPhiTemplate",";#phi^{*} [rad]; Events",32,0,M_PI,file_);
   add1DHistogram("h1DRhoTemplate",";#rho^{*} [rad]; Events",32,0.5*M_PI,M_PI,file_);
   add1DHistogram("h1DDPhiTemplate",";#Delta#phi^{*} [rad]; Events",32,-0.5*M_PI,0.5*M_PI,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::finalizeHistograms(int nRuns, float weight){

  plotHistograms();
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotHistograms(){


  std::string sysType = "";
  TCanvas* c = new TCanvas(TString::Format("Phi_%s",sysType.c_str()),
			   TString::Format("Phi_%s",sysType.c_str()),
			   460,500);

  TString hName = "h1DPhi";
  TH1F* h1D = this->get1DHistogram(hName.Data());

  h1D->Draw();

  c->Print(TString::Format("fig_png/Phi_%s.png",sysType.c_str()).Data());
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
