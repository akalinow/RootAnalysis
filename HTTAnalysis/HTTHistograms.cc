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

  //pileupCalc.py -i 15_12_2015.json --inputLumiJSON pileup_latest.txt --calcMode observed --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 MyDataPileupHistogram.root

  //brilcalc lumi --normtag /afs/cern.ch/user/c/cmsbril/public/normtag_json/OfflineNormtagV1.json -i lumiSummary.json
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 19    | 36   | 10256 | 10256 | 573248145.792     | 552672886.226    |

  //return 552672886.226e-6;//pb-1 data for NTUPLES_23_11_2015

  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i 15_12_2015.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 47    | 115  | 30707 | 30707 | 2135679014.929    | 2066764067.818   |
  //+-------+------+-------+-------+-------------------+------------------+

  //return 2066764067.818e-6;//pb-1 data for NTUPLES_15_12_2015


  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i 18_12_2015.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 47    | 116  | 32567 | 32567 | 2282566714.405    | 2207823354.548   |
  //+-------+------+-------+-------+-------------------+------------------+

  //return 2207823354.548e-6;//pb-1 data for NTUPLES_18_12_2015
  
  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i lumiSummary_Run2015C_16Dec2015_v1.json
  //+-------+------+------+------+-------------------+------------------+
  //| nfill | nrun | nls  | ncms | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+------+------+-------------------+------------------+
  //| 6     | 8    | 1099 | 1099 | 17767116.986      | 17225935.728     |
  //+-------+------+------+------+-------------------+------------------+
  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i lumiSummary_Run2015D_16Dec2015_v1.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 41    | 108  | 31592 | 31592 | 2185867468.728    | 2114239169.533   |
  //+-------+------+-------+-------+-------------------+------------------+
  //./.local/bin/brilcalc lumi --normtag  /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i lumiSummary_Run2016B_PromptReco_v2.json
  //+-------+------+-------+-------+-------------------+------------------+
  //| nfill | nrun | nls   | ncms  | totdelivered(/ub) | totrecorded(/ub) |
  //+-------+------+-------+-------+-------------------+------------------+
  //| 22    | 61   | 24919 | 24914 | 2151456680.808    | 2059830201.769   |
  //+-------+------+-------+-------+-------------------+------------------+  
  return (17225935.728 + 2114239169.533 + 2059830201.769)*1e-6;//pb-1 data for NTUPLES_20_06_2016
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
float HTTHistograms::getSampleNormalisation(std::string sampleName){

  std::string hName = "h1DStats"+sampleName;
  TH1F *hStats = 0;
  if(sampleName=="DYJets") hStats = get1D_DY_Histogram(hName.c_str());
  else if(sampleName=="WJets") hStats = get1D_WJet_Histogram(hName.c_str());
  else hStats = get1DHistogram(hName.c_str());

  if(!hStats) return 0;

  float genPresEff = 1.0;
  float recoPresEff = hStats->GetBinContent(3)/hStats->GetBinContent(2);
  float presEff = genPresEff*recoPresEff;
  float kFactor = 1.0;

  float crossSection = 1.0;
  int nEventsAnalysed = hStats->GetBinContent(1);

  ///Cross sections taken from
  if(sampleName=="DYJetsLowM"){
    //https://cmsweb.cern.ch/das/request?input=mcm%20prepid=SMP-RunIISpring15MiniAODv2-00016
    crossSection = 71600;
  }
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
    crossSection = 831.76*0.95; 
  }
  if(sampleName=="H"){
    ///mH = 125, gg fussion only
    //https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
    //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3
    crossSection = 44.14*6.32E-02; 
  }
  
  float weight = crossSection*presEff/nEventsAnalysed;
  if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;

  std::cout<<"Sample name: "<<sampleName<<" ";
  std::cout<<"Xsection: "<<crossSection<<" [pb] "<<" ";
  std::cout<<"Events analyzed: "<<nEventsAnalysed<<" ";
  //std::cout<<"Gen preselection efficiency: "<<genPresEff<<" ";
  std::cout<<"Reco preselection efficiency: "<<recoPresEff<<" ";
  //std::cout<<"External scaling: "<<kFactor<<" ";
  std::cout<<"Final weight: "<<weight<<std::endl;

  return weight;  

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TDirectory *myDir){

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TDirectory *myDir, const std::vector<std::string> & flavours){

  selectionFlavours_ = flavours;

  AnalysisHistograms::init(myDir);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::~HTTHistograms(){ }
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_DY_Histogram(const std::string& name){
  
  TString hName = name;
  hName.ReplaceAll("DYJets","DYJetsMuTau");
  TH1F *hDYJetsMuTau = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsMuMu");
  TH1F *hDYJetsMuMu = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsEE");
  TH1F *hDYJetsEE = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("DYJets","DYJetsOther");
  TH1F *hDYJetsOther = get1DHistogram(hName.Data());

  if(!hDYJetsMuTau) return 0;

  if(name.find("h1DStats")==std::string::npos){
    hDYJetsMuTau->Scale(0.85);//TEST AK
    hDYJetsMuTau->Scale(1.2);//TEST AK
  }
  
  if(hDYJetsMuMu) hDYJetsMuTau->Add(hDYJetsMuMu);
  if(hDYJetsEE) hDYJetsMuTau->Add(hDYJetsEE);
  if(hDYJetsOther) hDYJetsMuTau->Add(hDYJetsOther);
  hDYJetsMuTau->SetName(name.c_str());

  return hDYJetsMuTau;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_WJet_Histogram(const std::string& name){

  TString hName = name;
  hName.ReplaceAll("WJets","WJetsHT0");
  TH1F *hWJets0to100 = get1DHistogram(hName.Data());
  
  hName = name;
  hName.ReplaceAll("WJets","WJetsHT100to200");
  TH1F *hWJets100to200 = get1DHistogram(hName.Data());
    
  hName = name;
  hName.ReplaceAll("WJets","WJetsHT200to400");
  TH1F *hWJets200to400 = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("WJets","WJetsHT400to600");
  TH1F *hWJets400to600 = get1DHistogram(hName.Data());

  hName = name;
  hName.ReplaceAll("WJets","WJetsHT600toInf");
  TH1F *hWJets600toInf = get1DHistogram(hName.Data());

  TH1F *hWJets = 0;
  if(hWJets0to100) hWJets = (TH1F*)hWJets0to100->Clone(name.c_str());
  else if(hWJets100to200) hWJets = (TH1F*)hWJets100to200->Clone(name.c_str());
  else if(hWJets200to400) hWJets = (TH1F*)hWJets200to400->Clone(name.c_str());

  if(!hWJets) return 0;
  
  hWJets->Clear();
  if(hWJets0to100) hWJets->Add(hWJets0to100);
  if(hWJets100to200) hWJets->Add(hWJets100to200);
  if(hWJets200to400) hWJets->Add(hWJets200to400);
  if(hWJets400to600) hWJets->Add(hWJets400to600);
  if(hWJets600toInf) hWJets->Add(hWJets600toInf);
  
  return hWJets;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HTTHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";
  if(name.find("hProf")!=std::string::npos && name.find("VsMag")!=std::string::npos) templateName = "hProfVsMagTemplate";
  if(name.find("hProf")!=std::string::npos && name.find("VsPt")!=std::string::npos) templateName = "hProfVsPtTemplate";
  if(name.find("hProf")!=std::string::npos && name.find("VsCos")!=std::string::npos) templateName = "hProfVsCosTemplate";

  if(name.find("h1DNPV")!=std::string::npos) templateName = "h1DNPVTemplate";
  if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";
  if(name.find("h1DStats")!=std::string::npos) templateName = "h1DStatsTemplate";
  if(name.find("h1DPt")!=std::string::npos) templateName = "h1DPtTemplate";
  if(name.find("h1DEta")!=std::string::npos) templateName = "h1DEtaTemplate";
  if(name.find("h1DIso")!=std::string::npos) templateName = "h1DIsoTemplate";
  if(name.find("h1DPhi")!=std::string::npos) templateName = "h1DPhiTemplate";
  if(name.find("h1DCosPhi")!=std::string::npos) templateName = "h1DCosPhiTemplate";
  if(name.find("h1DCSVBtag")!=std::string::npos) templateName = "h1DCSVBtagTemplate";
  if(name.find("h1DID")!=std::string::npos) templateName = "h1DIDTemplate";
  if(name.find("h1DVxPull")!=std::string::npos) templateName = "h1DVxPullTemplate";
  if(name.find("h1DnPCA")!=std::string::npos) templateName = "h1DnPCATemplate";

  if(name.find("h2DVxPullVsNTrack")!=std::string::npos) templateName = "h2DVxPullVsNTrackTemplate";

  if(name.find("h1DyTau")!=std::string::npos) templateName = "h1DyTauTemplate";

  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){

   add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   add1DHistogram("h1DMassTemplate",";SVFit mass [GeV/c^{2}]; Events",50,0,200,file_);   
   add1DHistogram("h1DPtTemplate",";p_{T}; Events",20,0,100,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   add1DHistogram("h1DPhiTemplate",";#phi; Events",12,-M_PI,2*M_PI,file_);
   add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",10,-1.0,1.0,file_);
   add1DHistogram("h1DCSVBtagTemplate",";CSV btag; Events",20,0,1,file_);
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",10,0,0.5,file_);
   add1DHistogram("h1DIDTemplate",";ID; Events",30,-0.5,15.5,file_);
   add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",11,-0.01,0.01,file_);
   add1DHistogram("h1DyTauTemplate",";yTau; Events",30,-1.5,1.5,file_);
   add1DHistogram("h1DnPCATemplate",";#hat{n}_{RECO}>; Events",10,0,0.015,file_);

   add2DHistogram("h2DVxPullVsNTrackTemplate","",21,-0.5,20.5,11,-0.01,0.01,file_);
   
   addProfile("hProfVsMagTemplate","",10,0,0.015,file_);
   addProfile("hProfVsPtTemplate","",20,15,55,file_);
   addProfile("hProfVsCosTemplate","",20,-1,1,file_);
   
   histosInitialized_ = true;
 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  //plotCPhistograms(nRuns, weight);

  wselOSCorrection =  std::pair<float,float>(1,0);
  wselSSCorrection =  std::pair<float,float>(1,0);
  
  wselOSCorrection = getWNormalisation("wselOS");
  wselSSCorrection = getWNormalisation("wselSS");

  ///Control regions plots
  plotStack("Iso","qcdselOS");
  plotStack("Iso","qcdselSS");
  plotStack("StatsDecayMode","");
  
  plotStack("MassVis","qcdselSS");
  plotStack("StatsNJets30","qcdselSS");
  plotStack("CSVBtagLeadingJet","qcdselSS");
  
  plotStack("MassTrans","wselOS");  
  plotStack("MassTrans","wselSS");

  plotStack("MassTrans","ttselOS");  
  plotStack("MassTrans","ttselSS");
  plotStack("IsoMuon","ttselOS");

  plotStack("PtMET","mumuselOS");
  plotStack("PtMET","mumuselSS");
  
  ///Baseline selection plots
  plotStack("MassSV","");
  plotStack("MassVis","");  
  plotStack("MassTrans","");

  plotStack("PtMuon","");
  plotStack("EtaMuon","");
  plotStack("IsoMuon","");
  
  plotStack("PtTau","");  
  plotStack("EtaTau","");
  plotStack("IDTau","");
  plotStack("StatsDecayMode","");
  plotStack("PtTauLeadingTk","");
 
  plotStack("PhiMuon","");
  plotStack("PhiTau","");

  plotStack("PtMET","");  

  plotStack("StatsNJets30","");
  
  plotStack("PtLeadingJet","");
  plotStack("EtaLeadingJet","");
  plotStack("CSVBtagLeadingJet","");
  
  plotStack("PtLeadingBJet","");
  plotStack("EtaLeadingBJet","");

  plotStack("nPCAMuon","");
  plotStack("nPCATau","");
  plotStack("Phi_nVectors","");
  plotStack("Phi_nVectors","wselOS");

  plotStack("NPV","");
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotCPhistograms(int nRuns, float weight){

  plotProfiles("hProfMagVsPt_","H");
  plotProfiles("hProfMagVsPt_","A");
  plotProfiles("hProfMagVsPt_","DYJetsMuTau");

  plotProfiles("hProfMagVsCos_","H");
  plotProfiles("hProfMagVsCos_","A");
  plotProfiles("hProfMagVsCos_","DYJetsMuTau");

  plotProfiles("hProfPtVsMag_","H");
  plotProfiles("hProfPtVsMag_","A");
  plotProfiles("hProfPtVsMag_","DYJetsMuTau");

  plotProfiles("hProfRecoVsMagGen_","H");
  plotProfiles("hProfRecoVsMagGen_","A");
  plotProfiles("hProfRecoVsMagGen_","DYJetsMuTau");

  plotProfiles("hProfPhiVsMag_","H");
  plotProfiles("hProfPhiVsMag_","A");
  plotProfiles("hProfPhiVsMag_","DYJetsMuTau");

  plotVerticesPulls("h1DVxPullX_H");
  plotVerticesPulls("h1DVxPullY_H");
  plotVerticesPulls("h1DVxPullZ_H");

  plotVerticesPulls("h2DVxPullVsNTrackTrans_H");
  plotVerticesPulls("h2DVxPullVsNTrackLong_H");
  
  plotPhiDecayPlanes("Phi_nVectorsData");
  plotPhiDecayPlanes("Phi_nVectorsH");
  plotPhiDecayPlanes("Phi_nVectorsA");  
  plotPhiDecayPlanes("Phi_nVectorsDYJetsMuTau");
  plotPhiDecayPlanes("Phi_nVectorsWJetsHT0");
  
  plotPhiDecayPlanes("Phi_nVecIP_yTauPosData");
  plotPhiDecayPlanes("Phi_nVecIP_yTauNegData");
  plotPhiDecayPlanes("Phi_nVecIP_yTauPosH");
  plotPhiDecayPlanes("Phi_nVecIP_yTauNegH");
  plotPhiDecayPlanes("Phi_nVecIP_yTauPosA");
  plotPhiDecayPlanes("Phi_nVecIP_yTauNegA");
  plotPhiDecayPlanes("Phi_nVecIP_yTauPosDYJetsMuTau");
  plotPhiDecayPlanes("Phi_nVecIP_yTauNegDYJetsMuTau");
  plotPhiDecayPlanes("Phi_nVecIP_yTauPosWJetsHT0");
  plotPhiDecayPlanes("Phi_nVecIP_yTauNegWJetsHT0");

  plotSingleHistogram("yTauH");
  plotSingleHistogram("yTauA");
  plotSingleHistogram("yTauDYJetsMuTau");
  plotSingleHistogram("yTauWJetsHT0");

  plotPhiDecayPlanes("CosPhi_CosPositiveH");
  plotPhiDecayPlanes("CosPhi_CosNegativeH");

  plotPhiDecayPlanes("CosPhiNN_H");
  plotPhiDecayPlanes("CosPhiNN_A");
  plotPhiDecayPlanes("CosPhiNN_DYJetsMuTau");
  plotPhiDecayPlanes("CosPhiNN_WJetsHT0");

  plotPhiDecayPlanes("CosPhi_CosPositiveA");
  plotPhiDecayPlanes("CosPhi_CosNegativeA");

  plotnPCA("H");
  plotnPCA("A");
  plotnPCA("DYJetsMuTau");
  plotnPCA("WJetsHT0");
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotnPCA(const std::string & type){

  TH1F* h1DTau = get1DHistogram("h1DnPCATau"+type);
  TH1F* h1DMuon = get1DHistogram("h1DnPCAMuon"+type);
  if(!h1DTau || !h1DMuon) return;
  
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.12,0.35,0.22,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
   
  h1DTau->SetLineWidth(3);
  h1DTau->Scale(1.0/h1DTau->Integral(0,h1DTau->GetNbinsX()+1));

  h1DMuon->SetLineWidth(3);
  h1DMuon->SetLineColor(2);
  h1DMuon->Scale(1.0/h1DMuon->Integral(0,h1DMuon->GetNbinsX()+1));
  h1DMuon->GetYaxis()->SetTitleOffset(1.5);
  h1DMuon->SetStats(kFALSE);
  h1DMuon->SetYTitle("Events");
  h1DMuon->SetXTitle("|n_{RECO}|");

  h1DMuon->Draw();
  h1DTau->Draw("same");

  l.AddEntry(h1DTau,"hadronic tau");
  l.AddEntry(h1DMuon,"leptonic tau");
  l.Draw();
  
  c->Print(TString::Format("fig_png/nPCA_length_%s.png",type.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotVerticesPulls(const std::string & hName){
  
   TCanvas* c = new TCanvas("Vertices","Vertices resolutions",			   
			   460,500);

  TLegend l(0.1,0.7,0.25,0.85,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(hName.find("2D")!=std::string::npos){
    TProfile* hProfile_AOD = this->get2DHistogram((hName+"AODPV").c_str())->ProfileX();
    TProfile* hProfile_Refit = this->get2DHistogram((hName+"RefitPV").c_str())->ProfileX();

    if(!hProfile_AOD || !hProfile_Refit) return;

    hProfile_AOD->SetLineWidth(3);
    hProfile_Refit->SetLineWidth(3);

    hProfile_AOD->SetLineColor(1);
    hProfile_Refit->SetLineColor(4);
    
    hProfile_AOD->Draw();
    hProfile_Refit->Draw("same");

    float min = hProfile_AOD->GetMinimum();
    if(hProfile_Refit->GetMinimum()<min) min = hProfile_Refit->GetMinimum();
    hProfile_AOD->SetMinimum(0.95*min);

    l.AddEntry(hProfile_AOD,"from AOD");
    l.AddEntry(hProfile_Refit,"#splitline{refitted from}{mAOD, with BS}");
    l.Draw();
    
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
    return;
  }

  TH1F* h1D_AOD = this->get1DHistogram((hName+"AODPV").c_str());
  TH1F* h1D_Refit = this->get1DHistogram((hName+"RefitPV").c_str());
  
  if(h1D_AOD && h1D_Refit){
    
    h1D_AOD->SetLineWidth(3);
    h1D_Refit->SetLineWidth(3);
    ///
    h1D_AOD->SetLineColor(1);
    h1D_Refit->SetLineColor(4);
    ///
    h1D_AOD->Scale(1.0/h1D_AOD->Integral(0,h1D_AOD->GetNbinsX()+1));
    h1D_Refit->Scale(1.0/h1D_Refit->Integral(0,h1D_Refit->GetNbinsX()+1));
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
    if(h1D_Refit->GetMaximum()>max) max = h1D_Refit->GetMaximum();
    h1D_AOD->SetMaximum(1.05*max);
    h1D_AOD->Draw();
    h1D_Refit->Draw("same");

    l.AddEntry(h1D_AOD,"from AOD");
    l.AddEntry(h1D_Refit,"#splitline{refitted}{with BS}");
    l.Draw();
    
    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotProfiles(const std::string & hName,
				 const std::string & sysType){

  TProfile* h1DAOD = this->getProfile(hName+sysType+"AODPV");
  //TProfile* h1DGen = this->getProfile(hName+sysType+"GenNoOfflineSel");
  TProfile* h1DGen = this->getProfile(hName+sysType+"GenPV");  
  TProfile* h1DRefit = this->getProfile(hName+sysType+"RefitPV");

  if(!h1DGen || !h1DRefit || !h1DAOD) return;

  TCanvas c("AnyHistogram","AnyHistogram",460,500);			   
  c.SetLeftMargin(0.13);

  TLegend l(0.55,0.15,0.75,0.35,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(h1DGen && h1DRefit && h1DAOD){
    h1DGen->SetLineWidth(3);
    h1DAOD->SetLineWidth(3);
    h1DRefit->SetLineWidth(3);
    //
    h1DGen->SetLineColor(1);
    h1DAOD->SetLineColor(2);
    h1DRefit->SetLineColor(4);
    ///
    h1DGen->SetYTitle("<#hat{n}_{GEN} #bullet #hat{n}_{RECO}>");

    h1DGen->SetXTitle("|n_{GEN}|");
    h1DGen->GetYaxis()->SetTitleOffset(1.9);
    h1DGen->SetStats(kFALSE);

    if(hName.find("RecoVsMagGen")!=std::string::npos){
      h1DGen->SetYTitle("<|n_{RECO}|>");
      h1DGen->SetMinimum(0);
    }

    if(hName.find("MagVsPt")!=std::string::npos){
      h1DGen->SetYTitle("<|n_{GEN}|>");
      h1DGen->SetXTitle("p_{T}^{leading tk.}");
      h1DGen->SetMinimum(0);
    }
    if(hName.find("MagVsCos")!=std::string::npos){
      h1DGen->SetYTitle("<|n_{GEN}|>");
      h1DGen->SetXTitle("cos(#phi)");
    }
    
    h1DGen->Draw();
    h1DAOD->Draw("same");
    h1DRefit->Draw("same");

    if(hName.find("RecoVsMagGen")!=std::string::npos){
      TF1 *line=new TF1("line","x",0,0.014);
      line->Draw("same");
    }

    l.AddEntry(h1DGen,"Generator PV");
    l.AddEntry(h1DAOD,"AOD PV");
    l.AddEntry(h1DRefit,"Refitted PV");
    if(hName.find("RecoVsMagGen")!=std::string::npos) l.Draw();
    
    c.Print(TString::Format("fig_png/%s.png",(hName+sysType).c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotPhiDecayPlanes(const std::string & name){

  TCanvas aCanvas(TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  460,500);
  
  TLegend l(0.15,0.65,0.35,0.85,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TString hName = "h1D"+name+"RefitPV";
  TH1F* h1DRefitPV = get1DHistogram(hName.Data());

  hName = "h1D"+name+"AODPV";
  TH1F* h1DAODPV = get1DHistogram(hName.Data());
  
  hName = "h1D"+name+"GenPV";
  TH1F* h1DGenPV = get1DHistogram(hName.Data());

  hName = "h1D"+name+"GenNoOfflineSel";
  TH1F* h1DGen = get1DHistogram(hName.Data());

  if(h1DGen){
    h1DGen->SetLineWidth(4);
    h1DGen->Scale(1.0/h1DGen->Integral(0,h1DGen->GetNbinsX()+1));
    h1DGen->SetLineColor(1);
  }

  if(h1DAODPV){
    h1DAODPV->SetLineWidth(3);
    h1DAODPV->Scale(1.0/h1DAODPV->Integral(0,h1DAODPV->GetNbinsX()+1));
    h1DAODPV->SetLineColor(2);
  }

  if(h1DGenPV){
    h1DGenPV->SetLineWidth(3);
    h1DGenPV->Scale(1.0/h1DGenPV->Integral(0,h1DGenPV->GetNbinsX()+1));
    h1DGenPV->SetLineColor(3);
  }
  
  if(h1DRefitPV){
    h1DRefitPV->SetLineWidth(3);
    h1DRefitPV->SetLineColor(4);    
    h1DRefitPV->Scale(1.0/h1DRefitPV->Integral(0,h1DRefitPV->GetNbinsX()+1));
    h1DRefitPV->SetXTitle("#phi^{*}");
    h1DRefitPV->SetYTitle("Events");
    h1DRefitPV->SetTitle(name.c_str());
    h1DRefitPV->GetYaxis()->SetTitleOffset(1.4);
    h1DRefitPV->SetStats(kFALSE);
    //h1DRefitPV->GetXaxis()->SetRangeUser(0,M_PI);

    //h1DRefitPV->SetMaximum(1);
    h1DRefitPV->SetMinimum(0);
    h1DRefitPV->Draw("HISTO");

    h1DGenPV->Print();
    h1DGenPV->Print("all");
    
    l.AddEntry(h1DRefitPV,"nPCA with refit. PV");
    if(h1DGenPV){
      h1DGenPV->Draw("HISTO same");
      l.AddEntry(h1DGenPV,"nPCA with gen. PV");
    }    
    if(h1DAODPV){
      h1DAODPV->Draw("HISTO same");
      l.AddEntry(h1DAODPV,"nPCA with AOD PV");
    }
    if(h1DGen){
      h1DGen->Draw("HISTO same");
      l.AddEntry(h1DGen,"#splitline{PCA with gen. particles}{no offline selection}");
    }

    h1DRefitPV->Draw("HISTO same");
    l.Draw();
    aCanvas.Print(TString::Format("fig_png/%s.png",name.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
THStack*  HTTHistograms::plotStack(std::string varName, std::string selName){

  std::cout<<"--- Drawing THStack for variable: "<<varName
	   <<" selection: "<<selName<<std::endl;

  std::string hName = "h1D"+varName;
  TH1F *hHiggs = get1DHistogram((hName+"H"+selName).c_str());
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName).c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+selName).c_str());

  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+selName).c_str());
  TH1F *hDYJetsOther = get1DHistogram((hName+"DYJetsOther"+selName).c_str());
  TH1F *hDYJetsMuMu = get1DHistogram((hName+"DYJetsMuMu"+selName).c_str());
  TH1F *hDYJetsEE = get1DHistogram((hName+"DYJetsEE"+selName).c_str());
  TH1F *hDYJetsMuTau = get1DHistogram((hName+"DYJetsMuTau"+selName).c_str());  
  TH1F *hSoup = get1DHistogram((hName+"Data"+selName).c_str());
  pair<float,float> qcdOStoSS = getQCDOStoSS(selName, wselOSCorrection, wselSSCorrection);
  TH1F *hQCD = (TH1F*)getQCDbackground(varName,selName);
  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.

  if(!hSoup) return 0;
  
  if(!hQCD){
    hQCD = (TH1F*)hSoup->Clone((hName+"QCD"+selName).c_str()); hQCD->Reset();
  }
  if(!hWJets){
    hWJets = (TH1F*)hSoup->Clone((hName+"WJets"+selName).c_str()); hWJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+selName).c_str()); hDYJetsLowM->Reset();
  }  
  if(!hDYJetsOther){
    hDYJetsOther = (TH1F*)hSoup->Clone((hName+"hDYJetsOther"+selName).c_str()); hDYJetsOther->Reset();
  }
  if(!hDYJetsMuTau){
    hDYJetsMuTau = (TH1F*)hSoup->Clone((hName+"hDYJetsMuTau"+selName).c_str()); hDYJetsMuTau->Reset();
  }
  if(!hDYJetsMuMu){
    hDYJetsMuMu = (TH1F*)hSoup->Clone((hName+"hDYJetsMuMu"+selName).c_str()); hDYJetsMuMu->Reset();
  }
  if(!hDYJetsEE){
    hDYJetsEE = (TH1F*)hSoup->Clone((hName+"hDYJetsEE"+selName).c_str()); hDYJetsEE->Reset();
  }  
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+selName).c_str()); hTTbar->Reset();
  }
  if(!hHiggs){
    hHiggs = (TH1F*)hSoup->Clone((hName+"hH"+selName).c_str()); hHiggs->Reset();
  } 
  /////////////////////////////////////////////////////////////////

  float lumi = getLumi(); 
  std::string sampleName = "WJets";
  std::string WselType = "wselOS";
  if(selName.find("SS")!=std::string::npos) WselType = "wselSS";
  float weight = getSampleNormalisation(sampleName);
  ///TEST
  //pair<float,float> dataToMCScale = getWNormalisation(WselType);

  pair<float,float> dataToMCScale = wselOSCorrection;
  /////
  float scale = weight*lumi*dataToMCScale.first;
  hWJets->Scale(scale);

  sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM->Scale(scale);
  
  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsOther->Scale(scale);
  hDYJetsMuMu->Scale(1.2*scale);//TEST AK
  hDYJetsEE->Scale(scale);
  hDYJetsMuTau->Scale(0.85*scale);//TEST AK

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);

  sampleName = "H";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hHiggs->Scale(scale);
  //////////////////////////////////////////////////////
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);
  hSoup->SetMarkerStyle(20);

  hWJets->SetFillColor(kRed+2);
  hTTbar->SetFillColor(kBlue+2);
  hDYJetsOther->SetFillColor(kOrange-1);
  hDYJetsMuMu->SetFillColor(kOrange-3);
  hDYJetsEE->SetFillColor(kOrange-6);
  hDYJetsMuTau->SetFillColor(kOrange-9);
  hDYJetsLowM->SetFillColor(kOrange-7);

  hQCD->SetFillColor(kMagenta-10);
  hHiggs->SetFillColor(kCyan+4);

  hSoup->SetLineWidth(1);
  int rebinFactor = 1;  
  hSoup->Rebin(rebinFactor);
  hWJets->Rebin(rebinFactor);
  hTTbar->Rebin(rebinFactor);
  hDYJetsOther->Rebin(rebinFactor);
  hDYJetsMuMu->Rebin(rebinFactor);
  hDYJetsEE->Rebin(rebinFactor);
  hDYJetsMuTau->Rebin(rebinFactor);
  hDYJetsLowM->Rebin(rebinFactor);
  hHiggs->Rebin(rebinFactor);

  THStack *hs = new THStack("hs","Stacked histograms");      
  /////////
  hs->Add(hHiggs,"hist");    
  hs->Add(hQCD,"hist");
  hs->Add(hTTbar,"hist");
  hs->Add(hWJets,"hist");
  hs->Add(hDYJetsLowM,"hist");
  hs->Add(hDYJetsOther,"hist");
  hs->Add(hDYJetsMuMu,"hist");
  hs->Add(hDYJetsEE,"hist");
  hs->Add(hDYJetsMuTau,"hist");
  ////////
  TH1F *hMCSum = (TH1F*)hWJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJetsLowM);
  hMCSum->Add(hDYJetsMuTau);
  hMCSum->Add(hDYJetsMuMu);
  hMCSum->Add(hDYJetsEE);
  hMCSum->Add(hDYJetsOther);  
  hMCSum->Add(hWJets);
  hMCSum->Add(hTTbar);
  hMCSum->Add(hQCD);
  hMCSum->Add(hHiggs);

  if(!selName.size()) selName = "baseline";
  cout<<"Event count summary for selecion name: "<<selName<<std::endl;
  std::cout<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC: "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;  
  std::cout<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC TTbar: "<<hTTbar->Integral(0,hTTbar->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->mu tau: "<<hDYJetsMuTau->Integral(0,hDYJetsMuTau->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->mu mu: "<<hDYJetsMuMu->Integral(0,hDYJetsMuMu->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->e e: "<<hDYJetsEE->Integral(0,hDYJetsEE->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->ll(m<50): "<<hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC Z->other: "<<hDYJetsOther->Integral(0,hDYJetsOther->GetNbinsX()+1)<<std::endl;
  std::cout<<"MC H->tau tau: "<<hHiggs->Integral(0,hHiggs->GetNbinsX()+1)<<std::endl;  
  std::cout<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl; 
  std::cout<<"Correction factors:"<<std::endl;
  std::cout<<"QCD SS to OS: "<<qcdOStoSS.first<<" +- "<<qcdOStoSS.second<<std::endl;
  std::cout<<"W DATA to MC: "<<dataToMCScale.first<<" +- "<<dataToMCScale.second<<std::endl;
  std::cout<<"----------------------------------------"<<std::endl;

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

  if(!selName.size()) selName = "baseline";
  hs->SetTitle(("Variable: "+varName+" selection: "+selName).c_str());
  hs->SetMaximum(4400);
  hs->Draw("hist");
  hs->GetXaxis()->SetTitle(varName.c_str());
  hs->GetYaxis()->SetTitleOffset(1.4);
  hMCSum->SetFillColor(5);
  /////////
  float highEnd = 150;
  float lowEnd = -150;

  if(varName.find("Phi_")!=std::string::npos) lowEnd = 0;
    
  int binHigh = hs->GetXaxis()->FindBin(highEnd);  
  int binLow = hs->GetXaxis()->FindBin(lowEnd);

  if(hs->GetXaxis()->GetXmax()<highEnd || hs->GetXaxis()->GetXmax()>300) binHigh = hs->GetXaxis()->GetNbins();
  if(hs->GetXaxis()->GetXmin()>lowEnd) lowEnd = 1;

  hs->GetXaxis()->SetRange(binLow,binHigh);
  highEnd =  hs->GetXaxis()->GetBinUpEdge(binHigh);

  char yTitle[200];
  sprintf(yTitle,"Events/%2.1f",hSoup->GetXaxis()->GetBinWidth(1));
  hs->GetYaxis()->SetTitle(yTitle);
  
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
  leg->AddEntry(hDYJetsMuTau,"Z#rightarrow #mu #tau_{h}","f");
  leg->AddEntry(hDYJetsMuMu,"Z#rightarrow #mu #mu","f");
  leg->AddEntry(hDYJetsEE,"Z#rightarrow e e","f");
  leg->AddEntry(hDYJetsLowM,"Z#rightarrow ll(m<50)","f");
  leg->AddEntry(hDYJetsOther,"Z#rightarrow other","f");
  leg->AddEntry(hWJets,"W#rightarrow l #nu","f");
  leg->AddEntry(hTTbar,"TTbar","f");
  leg->AddEntry(hQCD,"QCD","f");
  leg->AddEntry(hHiggs,"H#rightarrow #tau #tau","f");
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
    else  hMCSum->SetBinContent(i,0);
    hMCSum->SetBinError(i,0);
  }
  hMCSum->SetLineWidth(3);
  hMCSum->SetMinimum(-5);
  hMCSum->SetMaximum(5);
  hMCSum->SetStats(kFALSE);
  hMCSum->Draw("hist");
  TLine *aLine = new TLine(hMCSum->GetXaxis()->GetXmin(),0.0,highEnd,0.0);
  aLine->SetLineColor(1);
  aLine->SetLineWidth(2);
  aLine->Draw();

  string plotName;
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + ".png";    
  else plotName = "fig_png/hTree_"+hName+Form("_%s",selName.c_str())+".png";
  c1->Print(plotName.c_str());

  if(hName.find_last_of("/")<string::npos) plotName = "fig_C/" + hName.substr(hName.find_last_of("/")) + ".C";    
  else plotName = "fig_C/hTree_"+hName+Form("_%s",selName.c_str())+".C";
  c1->Print(plotName.c_str()); 

  pad1->SetLogy(1);
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + "_LogY.png";    
  else plotName = "fig_png/hTree_"+hName+Form("_%s",selName.c_str())+"_LogY.png";
  c1->Print(plotName.c_str()); 

  std::cout<<"-------------------------------------------------------------"<<std::endl;

  return hs;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotSingleHistogram(std::string hName){

  TH1F* h1D = get1DHistogram(hName.c_str());
  if(!h1D) return;
  
  TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",			   
			   460,500);

  TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);
   
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
std::pair<float,float> HTTHistograms::getQCDOStoSS(std::string selName,
						   std::pair<float,float> wselOSCorrection,
						   std::pair<float,float> wselSSCorrection){

  std::cout<<"Calling method: "<<__func__<<std::endl;
  if(selName.find("SS")!=std::string::npos) return std::make_pair(1.0,0.0);

  std::string hName = "h1DIso";

  // SS selection
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS").c_str());
  TH1F *hDYJetsLowMSS = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS").c_str());
  TH1F *hDYJetsSS = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS").c_str());  
  TH1F *hTTSS = get1DHistogram((hName+"TTbar"+"qcdselSS").c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+"qcdselSS").c_str());
  TH1F *hSoupSSb = get1DHistogram((hName+"Data"+"qcdselSS").c_str());
  // OS selection
  TH1F *hWJetsOS = get1D_WJet_Histogram((hName+"WJets"+"qcdselOS").c_str());
  TH1F *hDYJetsLowMOS = get1DHistogram((hName+"DYJetsLowM"+"qcdselOS").c_str());
  TH1F *hDYJetsOS = get1D_DY_Histogram((hName+"DYJets"+"qcdselOS").c_str());
  TH1F *hTTOS = get1DHistogram((hName+"TTbar"+"qcdselOS").c_str());
  TH1F *hSoupOS = get1DHistogram((hName+"Data"+"qcdselOS").c_str());
  TH1F *hSoupOSb = get1DHistogram((hName+"Data"+"qcdselOS").c_str());

  if(!hSoupSS){
    std::cout<<"No data histograms for SS category!"<<std::endl;
    return  std::make_pair(1.0,0.0);
  }
  
  TH1F *hEmpty = (TH1F*)hSoupSS->Clone("hEmpty");
  hEmpty->Clear();
  
  if(!hSoupSS || !hSoupOS) return std::make_pair(1.0,0.0);
  if(!hDYJetsLowMSS) hDYJetsLowMSS = (TH1F*)hEmpty->Clone("hDYJetsLowMqcdselSS");
  if(!hDYJetsLowMOS) hDYJetsLowMOS = (TH1F*)hEmpty->Clone("hDYJetsLowMqcdselOS");

  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowMOS->Scale(scale);
  hDYJetsLowMSS->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsOS->Scale(scale);
  hDYJetsSS->Scale(scale);
  
  sampleName = "WJets";
  scale = getSampleNormalisation(sampleName);
  hWJetsOS->Scale(scale*wselOSCorrection.first);
  hWJetsSS->Scale(scale*wselSSCorrection.first);
  
  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTOS->Scale(scale);
  hTTSS->Scale(scale);
 
  ///Subtract backgrounds other than QCD using MC
  hSoupSS->Add(hWJetsSS,-1);
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);
  hSoupSS->Add(hTTSS,-1);
  
  hSoupOS->Add(hWJetsOS,-1);
  hSoupOS->Add(hDYJetsLowMOS,-1);
  hSoupOS->Add(hDYJetsOS,-1);
  hSoupOS->Add(hTTOS,-1);
  
  hSoupOS->Divide(hSoupSS);

  //funtion fitting
  TF1 *line=new TF1("line","[0]",0,2);
  line->SetParameter(0,1);
  TCanvas* c = new TCanvas("QCD_OStoSS","QCD_OStoSS",460,500);
  hSoupOS->SetLineWidth(3);
  hSoupOS->GetYaxis()->SetTitleOffset(1.4);
  hSoupOS->GetYaxis()->SetTitle("OS/SS");
  hSoupOS->GetXaxis()->SetTitle("muon relative isolation");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(11);
  hSoupOS->Draw();
  hSoupOS->Fit("line","","",0.2,0.3);
  c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
  c->Print(TString::Format("fig_C/%s.C",hName.c_str()).Data());

  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);

  std::cout<<"QCD OS/SS ratio: "<<param<<" +- "<<dparam<<std::endl;

  return std::make_pair(param, dparam);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::getQCDbackground(std::string varName, std::string selName,
				      std::pair<float,float> wselOSCorrection,
				      std::pair<float,float> wselSSCorrection){
				      

  std::cout<<"Calling method: "<<__func__<<std::endl;

  float qcdScale = getQCDOStoSS(selName, wselOSCorrection, wselSSCorrection).first;
  
  ///Not very clear and elegant. AK
  ///Need this to avoid resursive control region labels like
  ///qcdselSSqcdselOS
  if(selName.find("qcdsel")!=std::string::npos) selName = "";

  std::string hName = "h1D" + varName;
  // SS selection
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+"qcdselSS"+selName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+"qcdselSS"+selName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+"qcdselSS"+selName).c_str());
  TH1F *hTTbar = get1DHistogram((hName+"TTbar"+"qcdselSS"+selName).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+"qcdselSS"+selName).c_str());

  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.
  if(!hSoup) return 0;
  if(!hWJets){
    hWJets = (TH1F*)hSoup->Clone((hName+"WJets"+"qcdselSS").c_str()); hWJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hSoup->Clone((hName+"hDYJetsLowM"+"qcdselSS").c_str()); hDYJetsLowM->Reset();
  }
  if(!hDYJets){
    hDYJets = (TH1F*)hSoup->Clone((hName+"hDYJets"+"qcdselSS").c_str()); hDYJets->Reset();
  }
  if(!hTTbar){
    hTTbar = (TH1F*)hSoup->Clone((hName+"hTTbar"+"qcdselSS").c_str()); hTTbar->Reset();
  }
  //////////////////////////////////////////////////////////////////////
  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYJetsLowM";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "WJets";
  scale = getSampleNormalisation(sampleName)*lumi*wselSSCorrection.first;
  hWJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTTbar->Scale(scale);

  hSoup->SetName(("h1DQCDEstimate"+varName).c_str());
  hSoup->Add(hWJets,-1);
  hSoup->Add(hDYJetsLowM,-1);
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);

  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
  for(unsigned int iBinX=0;iBinX<=hSoup->GetNbinsX();++iBinX){
    if(hSoup->GetBinContent(iBinX)<3.0) hSoup->SetBinContent(iBinX,0);
  }
  
  hSoup->Scale(qcdScale);

  return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getWNormalisation(std::string selName){

  std::cout<<"Calling method: "<<__func__<<std::endl;

  std::string hName = "h1DMassTrans";
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+selName).c_str());
  TH1F *hDYJets = get1D_DY_Histogram((hName+"DYJets"+selName).c_str());
  TH1F *hDYJetsLowM = get1DHistogram((hName+"DYJetsLowM"+selName).c_str());
  TH1F *hTT = get1DHistogram((hName+"TTbar"+selName).c_str());
  TH1F *hQCD = (TH1F*)getQCDbackground("MassTrans",selName, wselOSCorrection, wselSSCorrection);
  TH1F *hSoup = get1DHistogram((hName+"Data"+selName).c_str());
  float lumi = getLumi();

  if(!hDYJets){
    hDYJets = (TH1F*)hWJets->Clone((hName+"hDYJets"+selName).c_str()); hDYJets->Reset();
  }
  if(!hDYJetsLowM){
    hDYJetsLowM = (TH1F*)hWJets->Clone((hName+"hDYJetsLowM"+selName).c_str()); hDYJetsLowM->Reset();
  }
  if(!hTT){
    hTT = (TH1F*)hWJets->Clone((hName+"hTTbar"+selName).c_str()); hTT->Reset();
  } 
		 
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hWJets->Scale(scale);

  sampleName = "DYJetsLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM->Scale(scale);
  
  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "TTbar";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hTT->Scale(scale);

  // Create a histogram with data minus backgrounds: DYJets, hTT, Other
  TH1F* datamtlo = (TH1F*)hSoup->Clone("datamtlo");
  datamtlo->Add(hDYJets,-1);
  datamtlo->Add(hDYJetsLowM,-1);
  datamtlo->Add(hTT,-1);
  datamtlo->Add(hQCD,-1);

  float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
  float intdata=datamtlo->Integral(0,datamtlo->GetNbinsX()+1);

  // Calculate weight
  weight=intdata/inthWJets;
  float dweight;
  float inthSoup = hSoup->Integral(0,hSoup->GetNbinsX()+1);
  float inthDYJetsLowM = hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1);
  float inthDYJets = hDYJets->Integral(0,hDYJets->GetNbinsX()+1);
  float inthTT = hTT->Integral(0,hTT->GetNbinsX()+1);
  float inthQCD = hQCD ? hQCD->Integral(0,hQCD->GetNbinsX()+1) : 0;
  float inthOther = 0;//hOther->Integral(0,hOther->GetNbinsX()+1);
  dweight=((inthSoup+inthDYJets+inthTT+inthOther)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
  dweight=sqrt(dweight);
  cout<<"Selecion name: "<<selName<<std::endl;
  cout<<"DATA: "<<inthSoup<<" DATA - MC(!WJets): "<<intdata<<" MC WJets "<<inthWJets
      <<" DYJets: "<<inthDYJets<<" DYJetsLowM: "<<inthDYJetsLowM
      <<" TTbar: "<<inthTT<<" QCD: "<<inthQCD
      <<" Other: "<<inthOther<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;
  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
