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
#include "TStyle.h"

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
bool CPHistograms::fillProfile(const std::string& name, float x, float val, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fillProfile(name,x, val,weight)){
    if(name.find("hProfPhiVsMag")!=std::string::npos) hTemplateName = "hProfPhiVsMagTemplate";
    
    this->addProfile(name,"",
		     this->getProfile(hTemplateName,true)->GetNbinsX(),
		     this->getProfile(hTemplateName,true)->GetXaxis()->GetXmin(),
		     this->getProfile(hTemplateName,true)->GetXaxis()->GetXmax(),
		     file_);
    
    return AnalysisHistograms::fillProfile(name,x,val,weight);
  }
  
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool CPHistograms::fill1DHistogram(const std::string& name, float val, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DDecayMode")!=std::string::npos) hTemplateName = "h1DDecayModeTemplate";
    if(name.find("h1DDeltaR")!=std::string::npos) hTemplateName = "h1DDeltaRTemplate";
    if(name.find("h1DCosPhi")!=std::string::npos) hTemplateName = "h1DCosPhiTemplate";
    if(name.find("h1DPhi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
    if(name.find("h1DRho")!=std::string::npos) hTemplateName = "h1DRhoTemplate";
    if(name.find("h1DIP")!=std::string::npos) hTemplateName = "h1IPTemplate";
    if(name.find("h1DVxPull")!=std::string::npos) hTemplateName = "h1DVxPullTemplate";
    if(get1DHistogram(hTemplateName,true)->GetXaxis()->IsVariableBinSize()){
      Float_t* binsArray = new Float_t[this->get1DHistogram(hTemplateName,true)->GetNbinsX()+1];
      for(unsigned int iBin=0;iBin<=this->get1DHistogram(hTemplateName,true)->GetNbinsX();++iBin){
	binsArray[iBin] = this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXbins()->At(iBin);
      }
      this->add1DHistogram(name,"",this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			   binsArray, file_);
      delete binsArray;      
    }
    else{
      this->add1DHistogram(name,"",
			   this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			   this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmin(),
			   this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmax(),
			   file_);
    }   
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool CPHistograms::fill2DHistogram(const std::string& name,
				   float valX, float valY, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill2DHistogram(name,valX,valY,weight)){
    if(name.find("h2DVxPullVsNTrack")!=std::string::npos) hTemplateName = "h2DVxPullVsNTrackTemplate";
    this->add2DHistogram(name,"",
			 this->get2DHistogram(hTemplateName)->GetNbinsX(),
			 this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmin(),
			 this->get2DHistogram(hTemplateName)->GetXaxis()->GetXmax(),
			 this->get2DHistogram(hTemplateName)->GetNbinsY(),
			 this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmin(),
			 this->get2DHistogram(hTemplateName)->GetYaxis()->GetXmax(),
			 file_);
    return AnalysisHistograms::fill2DHistogram(name,valX,valY,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   add1DHistogram("h1DDecayModeTemplate",";Decay mode; Events",19,-0.5,18.5,file_);
   add1DHistogram("h1DDeltaRTemplate",";#Delta R; Events",20,0,0.05,file_);
   add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",11,-0.01,0.01,file_);
   add1DHistogram("h1DPhiTemplate",";#phi^{*} [rad]; Events",4,0,M_PI,file_);//
   add1DHistogram("h1DRhoTemplate",";#rho^{*} [rad]; Events",18,M_PI-0.15,M_PI,file_);
   add1DHistogram("h1DDPhiTemplate",";#Delta#phi^{*} [rad]; Events",32,-0.5*M_PI,0.5*M_PI,file_);
   add1DHistogram("h1IPTemplate",";#phi^{*} [rad]; Events",100,0,1.0,file_);
   add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",10,-1.0,1.0,file_);

   add2DHistogram("h2DVxPullVsNTrackTemplate","",21,-0.5,20.5,11,-0.01,0.01,file_);

   addProfile("hProfPhiVsMagTemplate","",10,0,0.05,file_);

   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  plotProfiles("_h0");
  plotProfiles("_A0");
  plotProfiles("_Z0");
  
  plotPCAResolution("_h0");
  plotPCAResolution("_A0");
  plotPCAResolution("_Z0");

  plotVerticesPulls("h1DVxPullX_h0");
  plotVerticesPulls("h1DVxPullX_A0");
  plotVerticesPulls("h1DVxPullX_Z0");

  plotVerticesPulls("h1DVxPullY_h0");
  plotVerticesPulls("h1DVxPullY_A0");
  plotVerticesPulls("h1DVxPullY_Z0");

  plotVerticesPulls("h1DVxPullZ_h0");
  plotVerticesPulls("h1DVxPullZ_A0");
  plotVerticesPulls("h1DVxPullZ_Z0");

  plotVerticesPulls("h2DVxPullVsNTrackTrans_Z0");
  plotVerticesPulls("h2DVxPullVsNTrackLong_Z0");

  plotAnyHistogram("h1DDeltaRPlus_h0");
  plotAnyHistogram("h1DDeltaRPlus_A0");

  plotAnyHistogram("h1DDecayModePlus_h0");
  plotAnyHistogram("h1DDecayModeMinus_h0");
  
  std::string hName = "h1DPhi_nVectors";
  std::string sysType;
  for(auto it:my1Dhistograms_[0]){//FIX ME access to merged histograms map.
    if(it.first.find(hName)!=std::string::npos &&
       it.first.find("Template")==std::string::npos){
      sysType = it.first.substr(hName.size());
      plotHistograms(sysType);      
      std::string aType = sysType.substr(4,sysType.size());
      if(sysType.find("Z0")!=std::string::npos){
	plot_HAZ_Histograms("h1DPhi",aType);
	plot_HAZ_Histograms("h1DPhi_nVectors",aType);
	plot_HAZ_Histograms("h1DCosPhiNN",aType);
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

    float max = h_h->GetMaximum();
    if(h_A->GetMaximum()>max) max = h_A->GetMaximum();
    if(h_Z->GetMaximum()>max) max = h_Z->GetMaximum();	
    h_h->SetMinimum(0.0);
    h_h->SetMaximum(1.1*max);
    
    h_h->SetXTitle("#phi^{*}");
    if(name.Contains("CosPhiNN")) h_h->SetXTitle("#hat{n}_{RECO}^{#pi^{+}} #bullet #hat{n}_{RECO}^{#pi^{-}}");
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

  TLegend l(0.15,0.7,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TLatex aLatex(0,0,"");

  TString hName = "h1DPhi"+sysType;
  TString hName1 = "h1DPhi_nVectors"+sysType;
  TString hName2 = "h1DPhi_nVectors"+sysType;
  TString hName3 = "h1DPhi_nVectors"+sysType;
  
  if(hName.Contains("RECO")){
    hName = "h1DPhi_nVectors"+sysType;
    hName.ReplaceAll("RECO","GEN");
    hName2.ReplaceAll("RECO","AOD");
    hName3.ReplaceAll("RECO","RECOGEN");
  }
  
 
  TH1F* h1D = this->get1DHistogram(hName.Data());
  TH1F* h1DExp = this->get1DHistogram(hName1.Data());
  TH1F* h1DExp2 = this->get1DHistogram(hName2.Data());
  TH1F* h1DExp3 = this->get1DHistogram(hName3.Data());

  if(sysType.find("GEN")!=std::string::npos){
    h1DExp2 = 0;
    h1DExp3 = 0;
  }

  if(h1D){
    h1D->SetLineWidth(3);
    h1DExp->SetLineWidth(3);
    if(h1DExp2) h1DExp2->SetLineWidth(3);
    if(h1DExp3) h1DExp3->SetLineWidth(3);	
    h1D->Scale(1.0/h1D->Integral(0,h1D->GetNbinsX()+1));
    h1DExp->Scale(1.0/h1DExp->Integral(0,h1DExp->GetNbinsX()+1));    
    if(h1DExp2) h1DExp2->Scale(1.0/h1DExp2->Integral(0,h1DExp2->GetNbinsX()+1));
    if(h1DExp3) h1DExp3->Scale(1.0/h1DExp3->Integral(0,h1DExp3->GetNbinsX()+1));
    h1D->SetMinimum(0.0);
    h1D->SetMaximum(0.2);
    h1D->SetXTitle("#phi^{*}");
    h1D->SetYTitle("Events");
    h1D->GetYaxis()->SetTitleOffset(1.4);
    h1D->SetStats(kFALSE);
    h1D->Draw("L HIST");
    h1D->SetLineColor(1);
    h1DExp->SetLineColor(2);
    if(h1DExp2) h1DExp2->SetLineColor(4);
    if(h1DExp3) h1DExp3->SetLineColor(8);    
    h1DExp->Draw("same L HIST");
    if(h1DExp2) h1DExp2->Draw("same L HIST");
    if(h1DExp3) h1DExp3->Draw("same L HIST");
    ///
    if(hName1.Contains("GEN")) l.AddEntry(h1D,"tau momentum");
    if(hName1.Contains("RECO")) l.AddEntry(h1D,"PCA with GEN PV and SV");
    if(hName1.Contains("GEN")) l.AddEntry(h1DExp,"PCA with GEN PV and SV");
    if(hName1.Contains("RECO")) l.AddEntry(h1DExp,"PCA with refit reco vx.");
    if(h1DExp2) l.AddEntry(h1DExp2,"PCA with AOD vx.");
    if(h1DExp3) l.AddEntry(h1DExp3,"PCA with GEN vx.");
    l.Draw();
    aLatex.DrawLatex(0.05,0.01,sysType.c_str());
    
    ///
    c->Print(TString::Format("fig_png/Phi_%s.png",sysType.c_str()).Data());
  }
  /*
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
  */
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
  /*
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
  */
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotPCAResolution(const std::string & sysType){

  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.68,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TH1F* h1DAOD = this->get1DHistogram("h1DCosPhi_PCA_AOD"+sysType);
  TH1F* h1DGen = this->get1DHistogram("h1DCosPhi_PCA_Gen"+sysType);
  TH1F* h1DRefit = this->get1DHistogram("h1DCosPhi_PCA_Refit"+sysType);

  //TH1F* h1DGen = this->get1DHistogram("h1DCosPhi_PCA_Gen_A0");
  //TH1F* h1DRefit = this->get1DHistogram("h1DCosPhi_PCA_Gen_h0");

  //TH1F* h1DGen = this->get1DHistogram("h1DCosPhi_PCA_Refit_A0");
  //TH1F* h1DRefit = this->get1DHistogram("h1DCosPhi_PCA_Refit_h0");
  

  if(h1DGen && h1DRefit && h1DAOD){
    h1DGen->SetLineWidth(3);
    h1DAOD->SetLineWidth(3);
    h1DRefit->SetLineWidth(3);
    //
    h1DGen->Scale(1.0/h1DGen->Integral(0,h1DGen->GetNbinsX()+1));
    h1DAOD->Scale(1.0/h1DAOD->Integral(0,h1DAOD->GetNbinsX()+1));
    h1DRefit->Scale(1.0/h1DRefit->Integral(0,h1DRefit->GetNbinsX()+1));
    ///
    h1DGen->SetLineColor(1);
    h1DAOD->SetLineColor(2);
    h1DRefit->SetLineColor(4);
    ///
    h1DGen->SetYTitle("Events");
    h1DGen->SetXTitle("#hat{n}_{GEN}^{#pi^{+}} #bullet #hat{n}_{RECO}^{#pi^{+}}");
    h1DGen->GetYaxis()->SetTitleOffset(1.4);
    h1DGen->SetStats(kFALSE);
    float max = h1DGen->GetMaximum();
    if(h1DAOD->GetMaximum()>max) max = h1DAOD->GetMaximum();
    if(h1DRefit->GetMaximum()>max) max = h1DRefit->GetMaximum();
    h1DGen->SetMaximum(1.05*max);
    //h1DGen->SetMaximum(0.72);
    
    h1DGen->Draw();
    h1DAOD->Draw("same");
    h1DRefit->Draw("same");

    l.AddEntry(h1DGen,"Generator PV");
    l.AddEntry(h1DAOD,"AOD PV");
    l.AddEntry(h1DRefit,"Refitted PV");
    //l.AddEntry(h1DRefit,"Generator PV h0");
    l.Draw();
    //c->SetLogy();
    c->Print(TString::Format("fig_png/%s.png",("CosPhi_PCA"+sysType).c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotProfiles(const std::string & sysType){

  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.55,0.15,0.75,0.35,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TProfile* h1DAOD = this->getProfile("hProfPhiVsMag_AOD"+sysType);
  TProfile* h1DGen = this->getProfile("hProfPhiVsMag_Gen"+sysType);
  TProfile* h1DRefit = this->getProfile("hProfPhiVsMag_Refit"+sysType);

  if(h1DGen && h1DRefit && h1DAOD){
    h1DGen->SetLineWidth(3);
    h1DAOD->SetLineWidth(3);
    h1DRefit->SetLineWidth(3);
    //
    h1DGen->SetLineColor(1);
    h1DAOD->SetLineColor(2);
    h1DRefit->SetLineColor(4);
    ///
    h1DGen->SetYTitle("<#hat{n}_{GEN}^{#pi^{+}} #bullet #hat{n}_{RECO}^{#pi^{+}}>");
    h1DGen->SetXTitle("|n_{RECO}^{#pi^{+}}|");
    h1DGen->GetYaxis()->SetTitleOffset(1.4);
    h1DGen->SetStats(kFALSE);
    
    h1DGen->Draw();
    h1DAOD->Draw("same");
    h1DRefit->Draw("same");

    l.AddEntry(h1DGen,"Generator PV");
    l.AddEntry(h1DAOD,"AOD PV");
    l.AddEntry(h1DRefit,"Refitted PV");
    l.Draw();
    c->Print(TString::Format("fig_png/%s.png",("CosPhiVsMag_PCA"+sysType).c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void CPHistograms::plotVerticesPulls(const std::string & hName){
  
   TCanvas* c = new TCanvas("Vertices","Vertices resolutions",			   
			   460,500);

  TLegend l(0.15,0.68,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(hName.find("2D")!=std::string::npos){
    TProfile* hProfile_AOD = this->get2DHistogram((hName+"_AOD").c_str())->ProfileX();
    TProfile* hProfile_Refit = this->get2DHistogram((hName+"_RefitBS").c_str())->ProfileX();
    TProfile* hProfile_RefitNoBS = this->get2DHistogram((hName+"_RefitNoBS").c_str())->ProfileX();

    hProfile_AOD->SetLineWidth(3);
    hProfile_Refit->SetLineWidth(3);
    hProfile_RefitNoBS->SetLineWidth(3);

    hProfile_AOD->SetLineColor(1);
    hProfile_Refit->SetLineColor(4);
    hProfile_RefitNoBS->SetLineColor(8);
    
    hProfile_AOD->Draw();
    hProfile_Refit->Draw("same");
    hProfile_RefitNoBS->Draw("same");

    l.AddEntry(hProfile_AOD,"AOD weights");
    l.AddEntry(hProfile_Refit,"refitted with BS");
    l.AddEntry(hProfile_RefitNoBS,"refitted w/o #tau & BS");
    l.Draw();
    
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());

    return;
  }

  TH1F* h1D_AOD = this->get1DHistogram((hName+"_AOD").c_str());
  TH1F* h1D_PF = this->get1DHistogram((hName+"_PF").c_str());
  TH1F* h1D_Refit = this->get1DHistogram((hName+"_RefitBS").c_str());
  TH1F* h1D_RefitNoBS = this->get1DHistogram((hName+"_RefitNoBS").c_str());


  
  
  if(h1D_AOD && h1D_PF && h1D_Refit && h1D_RefitNoBS){
    
    h1D_AOD->SetLineWidth(3);
    h1D_PF->SetLineWidth(3);
    h1D_Refit->SetLineWidth(3);
    h1D_RefitNoBS->SetLineWidth(3);
    ///
    h1D_AOD->SetLineColor(1);
    h1D_PF->SetLineColor(2);
    h1D_Refit->SetLineColor(4);
    h1D_RefitNoBS->SetLineColor(8);
    ///
    h1D_AOD->Scale(1.0/h1D_AOD->Integral(0,h1D_AOD->GetNbinsX()+1));
    h1D_PF->Scale(1.0/h1D_PF->Integral(0,h1D_PF->GetNbinsX()+1));
    h1D_Refit->Scale(1.0/h1D_Refit->Integral(0,h1D_Refit->GetNbinsX()+1));
    h1D_RefitNoBS->Scale(1.0/h1D_RefitNoBS->Integral(0,h1D_RefitNoBS->GetNbinsX()+1));
    ///
    h1D_Refit->Fit("gaus");
    gStyle->SetOptFit(0001);
    gStyle->SetOptStat(0);
    ///
    h1D_AOD->SetYTitle("Events");
    h1D_AOD->SetXTitle("coordinate GEN - RECO [cm]");
    h1D_AOD->GetYaxis()->SetTitleOffset(1.4);
    h1D_AOD->SetStats(kFALSE);    
    ///
    float max =     h1D_AOD->GetMaximum();
    if(h1D_PF->GetMaximum()>max) max = h1D_PF->GetMaximum();
    if(h1D_Refit->GetMaximum()>max) max = h1D_Refit->GetMaximum();
    h1D_AOD->SetMaximum(1.05*max);
    h1D_AOD->Draw();
    //h1D_PF->Draw("same");
    h1D_Refit->Draw("same");
    //h1D_RefitNoBS->Draw("same");

    //TEST
    //h1D_Refit->Draw();
    //h1D_RefitNoBS->Draw("same");
    //TEST

    l.AddEntry(h1D_AOD,"AOD weights");
    //l.AddEntry(h1D_PF,"PF weights");
    l.AddEntry(h1D_Refit,"refitted with BS");
    //l.AddEntry(h1D_RefitNoBS,"refitted w/o #tau & BS");
    l.Draw();
    
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
    if(hName.find("DeltaR")!=std::string::npos) c->SetLogy();
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
