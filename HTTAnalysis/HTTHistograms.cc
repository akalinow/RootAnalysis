
#include <iostream>
#include <cmath>

#include "commonUtils.h"
#include "HTTHistograms.h"
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
float HTTHistograms::getLumi(){

  return 5579820.829*1E-6;//pb-1

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
float HTTHistograms::getSampleNormalisation(const std::string & sampleName){

  float genPresEff = 1.0;
  float recoPresEff = 1.0;
  float presEff = genPresEff*recoPresEff;
  float kFactor = 1.0;

  std::string hName = "h1DStats"+sampleName;
  TH1F *hStats = get1DHistogram(hName.c_str());
  
  float crossSection = 1.0;//TEST FIXME
  int nEventsAnalysed = hStats->GetBinContent(1);

  ///FIXME stupid if
  if(sampleName=="DY"){
    //xsection for 3xZ->mu mu in [pb]
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 3*2008.4; 
  }
  
  float weight = crossSection*presEff/nEventsAnalysed;
  if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;

  //std::cout<<"Sample name: "<<file->GetName()<<std::endl;
  std::cout<<"Mean cross section: "<<crossSection<<" [pb] "<<std::endl;
  std::cout<<"Number of events analyzed: "<<nEventsAnalysed<<std::endl;
  std::cout<<"Gen preselection efficiency: "<<genPresEff<<std::endl;
  std::cout<<"Reco preselection efficiency: "<<recoPresEff<<std::endl;
  std::cout<<"External scaling: "<<kFactor<<std::endl;
  std::cout<<"Final weight: "<<weight<<std::endl;

  return weight;  

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(std::string fileName, int opt){

  AnalysisHistograms::init(fileName);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TFileDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TFileDirectory *myDir, const std::vector<std::string> & flavours){
 selectionFlavours_ = flavours;

AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::~HTTHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
bool HTTHistograms::fill1DHistogram(const std::string& name, float val, float weight){
  
  std::string hTemplateName = "";
  if(!AnalysisHistograms::fill1DHistogram(name,val,weight)){
    if(name.find("h1DNPV")!=std::string::npos) hTemplateName = "h1DNPVTemplate";
    if(name.find("h1DSVfit")!=std::string::npos) hTemplateName = "h1DSVFitTemplate";
    if(name.find("h1DStats")!=std::string::npos) hTemplateName = "h1DStatsTemplate";
    std::cout<<"Adding histogram: "<<name<<" "<<file_<<" "<<file_->fullPath()<<std::endl;
    this->add1DHistogram(name,"",
			 this->get1DHistogram(hTemplateName,true)->GetNbinsX(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmin(),
			 this->get1DHistogram(hTemplateName,true)->GetXaxis()->GetXmax(),
			 file_);
    return AnalysisHistograms::fill1DHistogram(name,val,weight);
  }
  return true;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   //Make template histos
   std::cout<<"Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;
   add1DHistogram("h1DStatsTemplate","",10,0.5,10.5,file_);
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   add1DHistogram("h1DSVFitTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_);
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(int nRuns, float weight){

  plotStack("NPV",0);
  plotStack("SVfit",0);
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
THStack*  HTTHistograms::plotStack(std::string varName, int selType){

  std::string hName = "h1D"+varName;
  TH1F *hWJets = get1DHistogram((hName+"Data").c_str());//TEST should be WJets
  TH1F *hDYJets = get1DHistogram((hName+"Data").c_str());//TEST should be DY
  TH1F *hSoup = get1DHistogram((hName+"Data").c_str());

  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "Data";//TEST should be DY
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJets->Scale(scale);
  hWJets->Scale(0.0);
  //////////////////////////////////////////////////////
  
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);

  hWJets->SetFillColor(41);
  hDYJets->SetFillColor(28);

  hSoup->SetLineWidth(3);
  int rebinFactor = 1;
  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  
  THStack *hs = new THStack("hs","Stacked histograms");      
  /////////
  hs->Add(hWJets,"hist");
  hs->Add(hDYJets,"hist");
  ////////
  TH1F *hMCSum = (TH1F*)hWJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJets);
  hMCSum->Add(hWJets);

  std::cout<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC: "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;  
  std::cout<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->ll: "<<hDYJets->Integral(0,hDYJets->GetNbinsX()+1)<<std::endl;  

  TCanvas *c1 = getDefaultCanvas();
  c1->SetName("c1");
  c1->SetTitle("HTauTau analysis");
  c1->Divide(2);

  TPad *pad1 = (TPad*)c1->GetPad(1);
  TPad *pad2 = (TPad*)c1->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.29);
  pad1->SetRightMargin(0.22);
  pad2->SetRightMargin(0.22);
  pad2->SetFillStyle(4000);
  ///
  pad1->Draw();
  pad1->cd();

  hs->SetTitle("");
  hs->SetMaximum(4400);
  hs->Draw("hist");
  hMCSum->SetFillColor(5);
  /////////
  float highEnd = 150;
  float lowEnd = -150;
  
  int binHigh = hs->GetXaxis()->FindBin(highEnd);  
  int binLow = hs->GetXaxis()->FindBin(lowEnd);
  hs->GetXaxis()->SetRange(binLow,binHigh);
 
  char yTitle[200];
  sprintf(yTitle,"Events/%2.1f",hSoup->GetXaxis()->GetBinWidth(1));
  hs->GetYaxis()->SetTitle(yTitle);
  hs->SetTitle("");

  if(hName.find("svfit")!=std::string::npos)
    hs->GetXaxis()->SetTitle("SVFit mass [GeV/c^{2}]");
  
  float max = hs->GetMaximum();
  if(hSoup->GetMaximum()>max) max = hSoup->GetMaximum();

  hs->GetHistogram()->SetTitleOffset(1.0);
  hs->SetMaximum(1.1*max);
  hs->SetMinimum(0.1);

  hSoup->Draw("same");
  TH1F *hEmpty = new TH1F("hEmpty","",1,0,1);
  hEmpty->SetLineColor(10);
  hEmpty->SetFillColor(10);

  TLegend *leg = new TLegend(0.79,0.32,0.99,0.82,NULL,"brNDC");
  setupLegend(leg);
  leg->AddEntry(hSoup,"Data","lep");
  leg->AddEntry(hDYJets,"Z#rightarrow ll","f");
  leg->AddEntry(hWJets,"W#rightarrow l #nu","f");
  leg->SetHeader(Form("#int L = %.3f pb^{-1}",lumi));
  leg->Draw();

  float x = 0.6*(hs->GetXaxis()->GetXmax() - 
		 hs->GetXaxis()->GetXmin()) +
    hs->GetXaxis()->GetXmin(); 

  float y = 0.8*(max - 
		 hs->GetMinimum()) +
                 hs->GetMinimum(); 
  c1->cd();
  pad2->Draw();
  pad2->cd();

  hMCSum->GetXaxis()->SetRange(binLow,binHigh);
  hMCSum->SetTitle("");
  hMCSum->SetXTitle("");
  //hMCSum->SetYTitle("#frac{N_{obs} - N_{exp}}{#sqrt{N_{obs}}}");
  hMCSum->SetYTitle("#frac{N_{obs} - N_{exp}}{N_{obs}}");//TEST
  hMCSum->GetXaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetTitleSize(0.09);
  hMCSum->GetYaxis()->SetTitleOffset(0.5);
  hMCSum->Add(hSoup,-1);
  for(int i=0;i<hMCSum->GetNbinsX()+1;++i){
    //if(hSoup->GetBinContent(i)>0) hMCSum->SetBinContent(i,-hMCSum->GetBinContent(i)/sqrt(hSoup->GetBinContent(i)));
    if(hSoup->GetBinContent(i)>0) hMCSum->SetBinContent(i,-hMCSum->GetBinContent(i)/hSoup->GetBinContent(i)); //TEST
    else  hMCSum->SetBinContent(i,0);
    hMCSum->SetBinError(i,0);
  }
  hMCSum->SetLineWidth(3);
  hMCSum->SetMinimum(-2);
  hMCSum->SetMaximum(2);
  hMCSum->SetStats(kFALSE);
  hMCSum->Draw();
  TLine *aLine = new TLine(hMCSum->GetXaxis()->GetXmin(),0.0,highEnd,0.0);
  aLine->SetLineColor(1);
  aLine->SetLineWidth(2);
  aLine->Draw();


  string plotName;
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + ".png";    
  else plotName = "fig_png/hTree_"+hName+Form("_Sel%d",selType)+".png";
  c1->Print(plotName.c_str());

  if(hName.find_last_of("/")<string::npos) plotName = "fig_C/" + hName.substr(hName.find_last_of("/")) + ".C";    
  else plotName = "fig_C/hTree_"+hName+Form("_Sel%d",selType)+".C";
  c1->Print(plotName.c_str()); 

  pad1->SetLogy(1);
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + "_LogY.png";    
  else plotName = "fig_png/hTree_"+hName+Form("_Sel%d",selType)+"_LogY.png";
  c1->Print(plotName.c_str()); 

  return hs;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
