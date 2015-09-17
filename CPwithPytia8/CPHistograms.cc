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
    if(name.find("h1DDeltaR")!=std::string::npos) hTemplateName = "h1DDeltaRTemplate";
    if(name.find("h1DCosPhi")!=std::string::npos) hTemplateName = "h1DCosPhiTemplate";
    if(name.find("h1DPhi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
    if(name.find("h1DRho")!=std::string::npos) hTemplateName = "h1DRhoTemplate";
    if(name.find("h1DIP")!=std::string::npos) hTemplateName = "h1IPTemplate";
    if(name.find("h1DVxPull")!=std::string::npos) hTemplateName = "h1DVxPullTemplate";
    //std::cout<<"Adding histogram: "<<name<<std::endl;
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
   add1DHistogram("h1DDeltaRTemplate",";#Delta R; Events",100,0,0.5,file_);
   add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",50,-0.1,0.1,file_);
   add1DHistogram("h1DPhiTemplate",";#phi^{*} [rad]; Events",18,0,M_PI,file_);
   //add1DHistogram("h1DPhiTemplate",";#phi^{*} [rad]; Events",2*18,-1,1,file_);
   add1DHistogram("h1DRhoTemplate",";#rho^{*} [rad]; Events",18,M_PI-0.15,M_PI,file_);
   add1DHistogram("h1DDPhiTemplate",";#Delta#phi^{*} [rad]; Events",32,-0.5*M_PI,0.5*M_PI,file_);
   add1DHistogram("h1IPTemplate",";#phi^{*} [rad]; Events",100,0,1.0,file_);
   add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",100,0.0,1.0,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::finalizeHistograms(int nRuns, float weight){

  plotPCAResolution("h1DCosPhi_nGenRecoMinus");

  //plotAnyHistogram("h1DCosPhi_nGenRecoMinus");
  //plotAnyHistogram("h1DCosPhi_nGenRecoPlus");
  
  std::string hName = "h1DPhi_nVectors";
  std::string sysType;
  for(auto it:my1Dhistograms_){
    if(it.first.find(hName)!=std::string::npos &&
       it.first.find("Template")==std::string::npos){
      sysType = it.first.substr(hName.size());
      plotHistograms(sysType);
      std::string aType = sysType.substr(4,sysType.size());
      if(sysType.find("Z0")!=std::string::npos){
	plot_HAZ_Histograms("h1DPhi",aType);
	plot_HAZ_Histograms("h1DPhi_nVectors",aType);
	plot_HAZ_Histograms("h1DRho",aType);
	//plot_HAZ_Histograms("h1DPhi_Rho_y1y2Minus",aType);
	//plot_HAZ_Histograms("h1DPhi_Rho_y1y2MinusLAB",aType);
	//plot_HAZ_Histograms("h1DPhi_Rho_y1y2Plus",aType);
	//plot_HAZ_Histograms("h1DPhi_Rho_y1y2PlusLAB",aType);
      }
    }
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plot_HAZ_Histograms(const std::string & hName,
				       const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("%s_%s",hName.c_str(), sysType.c_str()),
			   TString::Format("%s_%s",hName.c_str(), sysType.c_str()),
			   460,500);

  TLegend l(0.15,0.7,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TLatex aLatex(0,0,"");

  TString name = hName+"_h0_"+sysType;  
  TH1F* h_h = this->get1DHistogram(name.Data());
  name = hName+"_A0_"+sysType;
  TH1F* h_A = this->get1DHistogram(name.Data());
  name = hName+"_Z0_"+sysType;
  TH1F* h_Z = this->get1DHistogram(name.Data());

  if(h_h && h_A && h_Z){
    h_h->SetLineWidth(3);
    h_A->SetLineWidth(3);
    h_Z->SetLineWidth(3);

    h_h->SetLineStyle(1);
    h_A->SetLineStyle(2);
    h_Z->SetLineStyle(3);

    h_h->Scale(1.0/h_h->Integral(0,h_h->GetNbinsX()+1));
    h_A->Scale(1.0/h_A->Integral(0,h_A->GetNbinsX()+1));
    h_Z->Scale(1.0/h_Z->Integral(0,h_Z->GetNbinsX()+1));

    h_h->SetMinimum(0.0);
    h_h->SetMaximum(0.15);
    h_h->SetXTitle("#phi^{*}");
    h_h->SetYTitle("Events");
    h_h->GetYaxis()->SetTitleOffset(1.4);
    h_h->SetStats(kFALSE);
    h_A->SetLineColor(2);
    h_Z->SetLineColor(3);
    h_h->Draw();
    h_A->Draw("same");
    h_Z->Draw("same");
    ///
    l.AddEntry(h_h,"SM h(125)");
    l.AddEntry(h_A,"MSSM A(125)");
    l.AddEntry(h_Z,"SM Z(91)");
    l.Draw();
    aLatex.DrawLatex(0.05,0.10,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/%s_h_A_Z_%s.png",hName.c_str(), sysType.c_str()).Data());
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

  //TString hName = "h1DCosPhi_collinearMinus"+sysType;
  //TString hName1 = "h1DCosPhi_collinearPlus"+sysType;
  
  TH1F* h1D = this->get1DHistogram(hName.Data());
  TH1F* h1DExp = this->get1DHistogram(hName1.Data());
  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.15);
    ///TEST
    //h1D->SetMinimum(1E-4);
    //h1D->SetMaximum(1.0);
    ///////////

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

    //c->SetLogy(); //TEST
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
  /*
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
    l.AddEntry(h1D,"y_{1} #dot y_{2}<0");
    l.AddEntry(h1DExp,"y_{1} #dot y_{2}>0");
    l.Draw();
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
    l.Draw();
    aLatex.DrawLatex(0.05,0.12,sysType.c_str());
    ///
    c->Print(TString::Format("fig_png/Rho_y1y2LAB_%s.png",sysType.c_str()).Data());
  }
  */

  hName = "h1DIP_PCA"+sysType;
  h1D = this->get1DHistogram(hName.Data());  
  hName = "h1DIP_3DIP"+sysType;
  TH1F* h1DMinus = this->get1DHistogram(hName.Data());
  if(h1D && h1DMinus){
    h1D->SetLineWidth(3);
    h1DMinus->SetLineWidth(3);
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DMinus->Scale(1.0/h1DMinus->Integral(0,h1DMinus->GetNbinsX()+1));
    h1D->SetXTitle("distance [cm]");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1DMinus->SetLineColor(2);
    h1D->Draw();
    h1DMinus->Draw("same");
    l.Clear();
    l.AddEntry(h1D,"PCA distance");
    l.AddEntry(h1DMinus,"3D flight path");
    l.Draw();
    c->SetLogy();
    c->Print(TString::Format("fig_png/3DIP_%s.png",sysType.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotPCAResolution(const std::string & hName){

  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.68,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TH1F* h1DAOD = this->get1DHistogram("h1DCosPhi_PCA_AOD");
  TH1F* h1DGen = this->get1DHistogram("h1DCosPhi_PCA_Gen");
  TH1F* h1DRefit = this->get1DHistogram("h1DCosPhi_PCA_Refit");

  if(h1DGen){
    h1DGen->SetLineWidth(3);
    h1DAOD->SetLineWidth(3);
    h1DRefit->SetLineWidth(3);
    //
    h1DGen->Scale(1.0/h1DGen->Integral(0,h1DGen->GetNbinsX()+1));
    h1DAOD->Scale(1.0/h1DAOD->Integral(0,h1DAOD->GetNbinsX()+1));
    h1DRefit->Scale(1.0/h1DRefit->Integral(0,h1DRefit->GetNbinsX()+1));
    ///
    h1DAOD->SetLineColor(2);
    h1DRefit->SetLineColor(4);
    ///
    h1DGen->SetYTitle("Events");
    h1DGen->SetXTitle("#hat{n}_{GEN}^{#pi^{+}} #bullet #hat{n}_{RECO}^{#pi^{+}}");
    h1DGen->GetYaxis()->SetTitleOffset(1.4);
    h1DGen->SetStats(kFALSE);
    h1DGen->Draw();
    h1DAOD->Draw("same");
    h1DRefit->Draw("same");

    l.AddEntry(h1DGen,"Generator PV");
    l.AddEntry(h1DAOD,"AOD PV");
    l.AddEntry(h1DRefit,"Refitted PV");
    l.Draw();
    //c->SetLogy();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotAnyHistogram(const std::string & hName){
  
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
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw();
    //c->SetLogy();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
