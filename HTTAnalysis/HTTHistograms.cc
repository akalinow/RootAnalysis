
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

  //brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i lumiSummary.json
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 19    | 36   | 10256 | 10256 | 573248145.792     | 552672886.226    |

  return 552672886.226e-6;//pb-1 data for NTUPLES_23_11_2015

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
  
  float crossSection = 1.0;
  int nEventsAnalysed = hStats->GetBinContent(1);

  ///FIXME stupid if
  ///Cross sections taken from
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
  if(sampleName=="DYJets"){
    //xsection for 3xZ->mu mu M50 in [pb]  
    crossSection = 3*2008.4; 
  }
  if(sampleName=="WJets"){
    //xsection for 3xW->mu nu in [pb]
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 3*20508.9;
  }
  if(sampleName=="TTbar"){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 831.76; 
  }
  
  float weight = crossSection*presEff/nEventsAnalysed;
  if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;

  std::cout<<"Sample name: "<<sampleName<<std::endl;
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
    if(name.find("h1DMass")!=std::string::npos) hTemplateName = "h1DMassTemplate";
    if(name.find("h1DStats")!=std::string::npos) hTemplateName = "h1DStatsTemplate";
    if(name.find("h1DPt")!=std::string::npos) hTemplateName = "h1DPtTemplate";
    if(name.find("h1DEta")!=std::string::npos) hTemplateName = "h1DEtaTemplate";
    if(name.find("h1DIso")!=std::string::npos) hTemplateName = "h1DIsoTemplate";
    if(name.find("h1DPhi")!=std::string::npos) hTemplateName = "h1DPhiTemplate";
    if(name.find("h1DMt")!=std::string::npos) hTemplateName = "h1DMtTemplate";
    std::cout<<"fill1DHistogram Adding histogram: "<<name<<" "<<file_<<" "<<file_->fullPath()<<std::endl;
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
   std::cout<<"defineHistograms Adding histogram: "<<file_<<" "<<file_->fullPath()<<std::endl;
   add1DHistogram("h1DStatsTemplate","",10,0.5,10.5,file_);
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_);
   add1DHistogram("h1DPtTemplate",";p_{T}; Events",20,0,100,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",20,0,2,file_);
   add1DHistogram("h1DPhiTemplate",";#phi; Events",30,-3,3,file_);
   add1DHistogram("h1DMtTemplate",";m_T; Events",50,0,200,file_);
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(int nRuns, float weight){

  //getWNormalisation("wselOS");
  //getWNormalisation("wselSS");
  //plotAnyHistogram("h1DMtTauWJets");
  //return;

  plotStack("NPV",0);

  plotStack("MassSV",0);
  plotStack("MassVis",0);  
  plotStack("MassTrans",0);
  

  plotStack("PtMuon",0);
  plotStack("EtaMuon",0);
  plotStack("IsoMuon",0);
  
  plotStack("PtTau",0);  
  plotStack("EtaTau",0);
 
  plotStack("PhiMuon",0);
  plotStack("PhiTau",0);
  plotStack("MtTau",0);


}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getQCDOStoSS(int selType){
  std::string hName = "h1DIso";
  // SS selection
  TH1F *hWJetsQ = get1DHistogram((hName+"WJets"+"qcdselSS").c_str());
  TH1F *hDYJetsQ = get1DHistogram((hName+"DY"+"qcdselSS").c_str());
  TH1F *hTTQ = get1DHistogram((hName+"TT"+"qcdselSS").c_str());
  TH1F *hOtherQ = get1DHistogram((hName+"Other"+"qcdselSS").c_str());
  TH1F *hSoupQ = get1DHistogram((hName+"Data"+"qcdselSS").c_str());
  // OS selection
  TH1F *hWJets = get1DHistogram((hName+"WJets"+"qcdselOS").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"DY"+"qcdselOS").c_str());
  TH1F *hTT = get1DHistogram((hName+"TT"+"qcdselOS").c_str());
  TH1F *hOther = get1DHistogram((hName+"Other"+"qcdselOS").c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselOS").c_str());
  // normalisation:: czy to jest poprawne? lumi jest dla wszystkich tych probek taka sama?
  float lumi = getLumi()/1000.0;
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DY";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJets->Scale(scale);
  hDYJetsQ->Scale(scale);
  sampleName = "WJets";
  scale = getSampleNormalisation(sampleName);
  hWJets->Scale(scale);
  hWJetsQ->Scale(scale);
  sampleName = "TT";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTT->Scale(scale);
  hOther->Scale(scale);
  hTTQ->Scale(scale);
  hOtherQ->Scale(scale);
  // OS and SS without background
  // NA RAZIE ZROBIC BEZ ODEJMOWANIA!!!
  TH1F* dataIsoSS= new TH1F("dataIsoSS","; Iso; Events",50,0,0.5);
  TH1F* dataIsoOS= new TH1F("dataIsoOS","; Iso",50,0,0.5);
  TH1F* dataIsoSSb= new TH1F("dataIsoSS","; Iso; Events",50,0,0.5);
  TH1F* dataIsoOSb= new TH1F("dataIsoOS","; Iso",50,0,0.5);
  dataIsoSSb->Add(hSoupQ,1);
  dataIsoOSb->Add(hSoup,1);
  dataIsoSS->Add(hSoupQ,1);
  dataIsoSS->Add(hWJetsQ,-1);
  dataIsoSS->Add(hDYJetsQ,-1);
  dataIsoSS->Add(hTTQ,-1);
  dataIsoSS->Add(hOtherQ,-1);
  dataIsoOS->Add(hSoup,1);
  dataIsoOS->Add(hWJets,-1);
  dataIsoOS->Add(hDYJets,-1);
  dataIsoOS->Add(hTT,-1);
  dataIsoOS->Add(hOther,-1);
  TH1F* diffSS= new TH1F("diffSS","; Iso; Events",50,0,0.5);
  TH1F* diffOS= new TH1F("diffOS","; Iso",50,0,0.5);
  diffSS->Add(dataIsoSSb,1);
  diffSS->Add(dataIsoSS,-1);
  diffOS->Add(dataIsoOSb,1);
  diffOS->Add(dataIsoOS,-1);
  // roznice w zliczeniach:
  cout<<"!!!! roznica w zliczeniach miedzy histogramami; diffSS "<<diffSS->Integral(30,50)<<" diffOS "<<diffOS->Integral(30,50)<<endl;
  cout<<"!!!! same histogramy; dataIsoSSb "<<dataIsoSSb->Integral(30,50)<<" dataIsoSS "<<dataIsoSS->Integral(30,50)<<" dataIsoOSb "<<dataIsoOSb->Integral(30,50)<<" dataIsoOS "<<dataIsoOS->Integral(30,50)<<endl;
  dataIsoOSb->Divide(dataIsoSSb);
  dataIsoOS->Divide(dataIsoSS);
  //funtion fitting
  TF1 *line2=new TF1("line2","[0]",0.3,0.5);
  line2->SetParameter(0,1);
  dataIsoOSb->Fit("line2","","",0.3,0.5);
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",460,500);
  dataIsoOSb->Draw();
  line2->SetLineColor(kGreen);
  line2->Draw("same");
  c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  // fitting do odjetych
  hName=hName+"BezTla";
  TF1 *line=new TF1("line","[0]",0.3,0.5);
  line->SetParameter(0,1);
  dataIsoOS->Fit("line","","",0.3,0.5);
  TCanvas* c1 = new TCanvas("AnyHistogram","AnyHistogram",460,500);
  dataIsoOS->Draw();
  line->SetLineColor(kGreen);
  line->Draw("same");
  c1->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  // OS/SS
  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);
  return std::make_pair(param, dparam);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1* HTTHistograms::getQCDbackground(std::string varName, int selType){

  std::string hName = "h1D" + varName;

  // SS selection
  TH1F *hWJets = get1DHistogram((hName+"WJets"+"qcdselSS").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"DYJets"+"qcdselSS").c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"qcdselSS").c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselSS").c_str());

  std::cout << "Data SS integral: " <<  hSoup->Integral() << std::endl;
  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "WJets";
  float dataToMCScale = getWNormalisation("wselSS").first;
  scale = getSampleNormalisation(sampleName)*lumi*dataToMCScale;
  hWJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);

  // OS and SS without background
  hSoup->SetName(("h1DQCDEstimate"+varName).c_str());
  hSoup->Add(hWJets,-1);
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);

  // scale background FIXME: has to be computed.
  scale = 1.06;
  scale= 1.0;
  hSoup->Scale(scale);

  return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
THStack*  HTTHistograms::plotStack(std::string varName, int selType){

  std::string hName = "h1D"+varName;
  TH1F *hWJets = get1DHistogram((hName+"WJets").c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar").c_str());
  TH1F *hDYJets = get1DHistogram((hName+"DYJets").c_str());   
  TH1F *hSoup = get1DHistogram((hName+"Data").c_str());
  TH1F *hQCD = (TH1F*)getQCDbackground(varName,0);        

  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "WJets";
  std::string WselType = "wselOS";
  float weight = getSampleNormalisation(sampleName);
  float dataToMCScale = getWNormalisation(WselType).first;
  float scale = weight*lumi*dataToMCScale;
  hWJets->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);
  //////////////////////////////////////////////////////
  
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);
  hSoup->SetMarkerStyle(20);
  //hSoup->SetMarkerSize(3);

  hWJets->SetFillColor(kRed+2);
  hTTbar->SetFillColor(kBlue+2);
  hDYJets->SetFillColor(kOrange-4);
  hQCD->SetFillColor(kMagenta-10);

  hSoup->SetLineWidth(1);
  int rebinFactor = 1;
  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJets->Rebin(rebinFactor);
  
  THStack *hs = new THStack("hs","Stacked histograms");      
  /////////
  hs->Add(hQCD,"hist");
  hs->Add(hWJets,"hist");
  hs->Add(hTTbar,"hist");
  hs->Add(hDYJets,"hist");
  ////////
  TH1F *hMCSum = (TH1F*)hWJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJets);
  hMCSum->Add(hWJets);
  hMCSum->Add(hTTbar);
  hMCSum->Add(hQCD);

  std::cout<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC: "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;  
  std::cout<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC TTbar: "<<hTTbar->Integral(0,hTTbar->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->ll: "<<hDYJets->Integral(0,hDYJets->GetNbinsX()+1)<<std::endl;  
  std::cout<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl; 

  TCanvas *c1 = getDefaultCanvas();
  c1->SetName("c1");
  c1->SetTitle("HTauTau analysis");
  c1->Divide(2);

  TPad *pad1 = (TPad*)c1->GetPad(1);
  TPad *pad2 = (TPad*)c1->GetPad(2);
  pad1->SetPad(0.01,0.29,0.99,0.99);
  pad2->SetPad(0.01,0.01,0.99,0.29);
  pad1->SetRightMargin(0.23);
  pad2->SetRightMargin(0.23);
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
  
  if(hs->GetXaxis()->GetXmax()<highEnd) highEnd = hs->GetXaxis()->GetXmax();
  if(hs->GetXaxis()->GetXmin()>lowEnd) lowEnd = hs->GetXaxis()->GetXmin();
 
  int binHigh = hs->GetXaxis()->FindBin(highEnd);  
  int binLow = hs->GetXaxis()->FindBin(lowEnd);
  hs->GetXaxis()->SetRange(binLow,binHigh);
 
  char yTitle[200];
  sprintf(yTitle,"Events/%2.1f",hSoup->GetXaxis()->GetBinWidth(1));
  hs->GetYaxis()->SetTitle(yTitle);
  hs->SetTitle("");

  if(hName.find("MassSV")!=std::string::npos)
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
  leg->AddEntry(hTTbar,"TTbar","f");
  leg->AddEntry(hQCD,"QCD","f");
  leg->SetHeader(Form("#int L = %.2f fb^{-1}",lumi/1000));
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
  hMCSum->SetYTitle("#frac{N_{obs} - N_{exp}}{#sqrt{N_{obs}}}");
  hMCSum->GetXaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetLabelSize(0.09);
  hMCSum->GetYaxis()->SetTitleSize(0.09);
  hMCSum->GetYaxis()->SetTitleOffset(0.5);
  hMCSum->Add(hSoup,-1);
  for(int i=0;i<hMCSum->GetNbinsX()+1;++i){
    if(hSoup->GetBinContent(i)>0) hMCSum->SetBinContent(i,-hMCSum->GetBinContent(i)/sqrt(hSoup->GetBinContent(i)));
    //if(hSoup->GetBinContent(i)>0) hMCSum->SetBinContent(i,-hMCSum->GetBinContent(i)/hSoup->GetBinContent(i)); //TEST
    else  hMCSum->SetBinContent(i,0);
    hMCSum->SetBinError(i,0);
  }
  hMCSum->SetLineWidth(3);
  hMCSum->SetMinimum(-5);
  hMCSum->SetMaximum(5);
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
void HTTHistograms::plotAnyHistogram(const std::string & hName){
  
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
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getWNormalisation(const std::string & selType){

  std::string hName = "h1DMt";
  TH1F *hWJets = get1DHistogram((hName+"WJets"+selType).c_str());
  TH1F *hDYJets = get1DHistogram((hName+"DYJets"+selType).c_str());
  TH1F *hTT = get1DHistogram((hName+"TTbar"+selType).c_str());
  //TH1F *hOther = get1DHistogram((hName+"Other"+selType).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+selType).c_str());
  float lumi = getLumi();
  /* TEST 
  TH1F *hDYJets = get1DHistogram((hName+"WJets"+selType).c_str());
  TH1F *hTT = get1DHistogram((hName+"WJets"+selType).c_str());
  hDYJets->Reset();
  hTT->Reset();
  */
 
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hWJets->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTT->Scale(scale);
  //hOther->Scale(scale);

  // Create a histogram with data minus bacgroungs: DYJets, hTT, Other
  TH1F* datamtlo = (TH1F*)hSoup->Clone("datamtlo");
  datamtlo->Add(hDYJets,-1);
  datamtlo->Add(hTT,-1);
  //datamtlo->Add(hOther,-1);

  float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
  float intdata=datamtlo->Integral(0,datamtlo->GetNbinsX()+1);

  // Calculate weight
  weight=intdata/inthWJets;
  float dweight;
  float inthSoup = hSoup->Integral(0,hSoup->GetNbinsX()+1);
  float inthDYJets = hDYJets->Integral(0,hDYJets->GetNbinsX()+1);
  float inthTT = hTT->Integral(0,hTT->GetNbinsX()+1);
  float inthOther = 0;//hOther->Integral(0,hOther->GetNbinsX()+1);
  dweight=((inthSoup+inthDYJets+inthTT+inthOther)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
  dweight=sqrt(dweight);
  cout<<"Selectino type: "<<selType<<std::endl;
  cout<<"DATA: "<<inthSoup<<" DYJets: "<<inthDYJets<<" TTbar: "<<inthTT<<" Other: "<<inthOther<<endl;
  cout<<"MC WJets "<<inthWJets<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;
  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
