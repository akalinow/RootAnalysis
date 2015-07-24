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
#include "TLatex.h"

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
    if(name.find("h1DIP")!=std::string::npos) hTemplateName = "h1IPTemplate";
    if(name.find("h1DVxPull")!=std::string::npos) hTemplateName = "h1DVxPullTemplate";
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
   add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",50,-0.1,0.1,file_);
   add1DHistogram("h1DPhiTemplate",";#phi^{*} [rad]; Events",18,0,M_PI,file_);
   //add1DHistogram("h1DRhoTemplate",";#rho^{*} [rad]; Events",18,0,M_PI,file_);
   add1DHistogram("h1DRhoTemplate",";#rho^{*} [rad]; Events",18,M_PI-0.15,M_PI,file_);
   add1DHistogram("h1DDPhiTemplate",";#Delta#phi^{*} [rad]; Events",32,-0.5*M_PI,0.5*M_PI,file_);
   add1DHistogram("h1IPTemplate",";#phi^{*} [rad]; Events",20,0,1.0,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::finalizeHistograms(int nRuns, float weight){

  std::string hName = "h1DPhi_nVectors";
  std::string sysType;
  for(auto it:my1Dhistograms_){
    if(it.first.find(hName)!=std::string::npos &&
       it.first.find("Template")==std::string::npos){
      sysType = it.first.substr(hName.size());
      std::cout<<sysType<<std::endl;
      plotHistograms(sysType);
    }
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotHistograms(const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("Phi_%s",sysType.c_str()),
			   TString::Format("Phi_%s",sysType.c_str()),
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TLatex aLatex(0,0,"");

  TString hName = "h1DPhi"+sysType;
  TString hName1 = "h1DPhi_nVectors"+sysType;
  TH1F* h1D = this->get1DHistogram(hName.Data());
  TH1F* h1DExp = this->get1DHistogram(hName1.Data());
  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.15);
    h1D->SetXTitle("#phi^{*}");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    h1DExp->SetLineColor(2);
    h1DExp->Draw("same");
    ///
    l.AddEntry(h1D,"tau momentum");
    l.AddEntry(h1DExp,"impact param. vectors");
    l.Draw();
    aLatex.DrawLatex(0.05,0.12,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/Phi_%s.png",sysType.c_str()).Data());
  }

  hName = "h1DRho"+sysType;
  hName1 = "h1DRho_nVectors"+sysType;
  h1D = this->get1DHistogram(hName.Data());
  h1DExp = this->get1DHistogram(hName1.Data());
  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));    
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.2);
    h1D->SetXTitle("#rho^{*}");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    h1DExp->SetLineColor(2);
    h1DExp->Draw("same");
    ///
    l.Draw();
    aLatex.DrawLatex(0.05,0.8,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/Rho_%s.png",sysType.c_str()).Data());
  }

  hName = "h1DPhi_Rho_y1y2Minus"+sysType;
  hName1 = "h1DPhi_Rho_y1y2Plus"+sysType;
  h1D = this->get1DHistogram(hName.Data());
  h1DExp = this->get1DHistogram(hName1.Data());
  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.14);
    h1D->SetXTitle("#rho^{*}");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    h1DExp->SetLineColor(2);
    h1DExp->Draw("same");
    ///
    l.Clear();
    l.AddEntry(h1D,"y_{1} #cdot y_{2}<0");
    l.AddEntry(h1DExp,"y_{1} #cdot y_{2}>0");
    aLatex.DrawLatex(0.05,0.12,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/Rho_y1y2_%s.png",sysType.c_str()).Data());
  }


  hName = "h1DPhi_Rho_y1y2MinusLAB"+sysType;
  hName1 = "h1DPhi_Rho_y1y2PlusLAB"+sysType;
  h1D = this->get1DHistogram(hName.Data());
  h1DExp = this->get1DHistogram(hName1.Data());
  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.14);
    h1D->SetXTitle("#rho^{*}");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    h1DExp->SetLineColor(2);
    h1DExp->Draw("same");
    ///
    l.Clear();
    l.AddEntry(h1D,"y_{1}^{LAB} #cdot y_{2}^{LAB}<0");
    l.AddEntry(h1DExp,"y_{1}^{LAB} #cdot y_{2}^{LAB}>0");
    aLatex.DrawLatex(0.05,0.12,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/Rho_y1y2LAB_%s.png",sysType.c_str()).Data());
  }
  

  hName = "h1DIPPlus"+sysType;
  h1D = this->get1DHistogram(hName.Data());  
  hName = "h1DIPMinus"+sysType;
  TH1F* h1DMinus = this->get1DHistogram(hName.Data());
  if(h1D && h1DMinus){
    h1D->SetLineWidth(3);
    h1DMinus->SetLineWidth(3);
    h1D->SetXTitle("3D IP [mm]");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1DMinus->SetLineColor(2);
    h1D->Draw();
    h1DMinus->Draw("same");
    l.Clear();
    l.AddEntry(h1D,"#tau^{+}");
    l.AddEntry(h1DMinus,"#tau^{-}");
    l.Draw();  
    c->Print(TString::Format("fig_png/3DIP_%s.png",sysType.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
