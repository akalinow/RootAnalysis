#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "HTTHistograms.h"
#include "HTTAnalyzer.h"
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

  //./.local/bin/brilcalc lumi --normtag ~lumipro/public/normtag_file/OfflineNormtagV2.json -i lumiSummary_Run2015C_16Dec2015_v1.json
  float run2015C = 17225935.728*1E-6;
  float run2015D = 2114239169.533*1E-6;

  //./.local/bin/brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i processedLumis_SingleMuon.json

  float run2016BPromptReco = 5923961370.727;
  float run2016BReReco = 5933308579.501;

  float run2016CPromptReco = 2645968083.093;
  float run2016CReReco = 2645968083.093;

  float run2016DPromptReco = 4353448810.554;
  float run2016DReReco = 4353448810.554;

  float run2016EReReco = 4049255306.406;
  float run2016FReReco = 3160088401.247;
  float run2016GReReco = 7554453635.136;

  float run2016HPromptReco_v2 = 8545039541.475;
  float run2016HPromptReco_v3 = 216782873.203;

  float run2016 = run2016BReReco + run2016CReReco +
    run2016DReReco + run2016EReReco +
    run2016FReReco + run2016GReReco +
    run2016HPromptReco_v2 + run2016HPromptReco_v3;

  return run2016*1E-6;//pb-1 data for NTUPLES_05_12_2016
  //return 36446609816.686*1E-6;//pb-1 data for NTUPLES_05_12_2016
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
float HTTHistograms::getSampleNormalisation(std::string sampleName){

  std::string hName = "h1DStats"+sampleName;
  TH1F *hStats = new TH1F("h1DStats","",11,-0.5,10.5);
  hStats->SetBinContent(1,1);
  hStats->SetBinContent(2,1);
  hStats->SetBinContent(3,1);

  if(sampleName.find("DYJetsMatch")!=std::string::npos || sampleName=="DYJets") {/*DYJets50 are normalised for analysed events and preselection for each nJets sample*/}
  else if(sampleName.find("WJets")!=std::string::npos) {/*WJets are normalised for analysed events and preselection for each nJets sample*/ }
  else if(sampleName=="ST") {/*WJets are normalised for analysed events and preselection for each nJets sample*/ }
  else hStats = get1DHistogram(hName.c_str());

  if(!hStats) hStats = get1D_TauMatchJetSum(hName,true,false);

  if(!hStats) return 0;

  float genPresEff = 1.0;
  float recoPresEff = hStats->GetBinContent(3)/hStats->GetBinContent(2);
  float presEff = genPresEff*recoPresEff;
  float kFactor = 1.0;

  float crossSection = 1.0;
  int nEventsAnalysed = hStats->GetBinContent(1);

  //test HACK fixme!!!
  std::vector<std::string> sampleNames = {"TTbar", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};
  for(auto sampleNameTmp:sampleNames){
        if(sampleName.find(sampleNameTmp+"Match")!=std::string::npos) {hStats = get1D_TauMatchJetSum(hName, true, false);nEventsAnalysed=hStats->GetBinContent(1);}
        }
  //test

  ///Cross sections taken from
  if(sampleName=="DYLowM"){
    //https://cmsweb.cern.ch/das/request?input=mcm%20prepid=SMP-RunIISpring15MiniAODv2-00016
    crossSection = 71600;
  }
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
  if(sampleName.find("DYJetsMatch")!=std::string::npos || sampleName=="DYJets"){
    //xsection for 3xZ->mu mu M50 in [pb]
    crossSection = 3*1921.8;
  }
  if(sampleName.find("WJets")!=std::string::npos){
    //xsection for 3xW->mu nu in [pb]
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 3*20508.9;
  }
  if(sampleName.find("TTbar")!=std::string::npos){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
    crossSection = 831.76*ttScale;
  }

  //https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
  //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2014#Higgs_2_fermions
  //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014
  //Xsection for mass!=125 are calculated using luminosity ratio, and cross section for 8 TeV

  if(sampleName=="ggH120") crossSection = 5.222E+01*6.981E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
  if(sampleName=="ggH125") crossSection = 4.858E+01*6.272E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
  if(sampleName=="ggH130") crossSection = 4.531E+01*5.411E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR

  if(sampleName=="qqH120") crossSection = 1.676E+00*6.981E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
  if(sampleName=="qqH125") crossSection = 1.601E+00*6.272E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
  if(sampleName=="qqH130") crossSection = 1.531E+00*5.411E-02;//CERNYellowReportPageAt13TeV*CERNYellowReportPageBR

  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
  if(sampleName=="WplusH120") crossSection = 1.565*0.0698*0.5;
  if(sampleName=="WplusH125") crossSection = 1.373*0.0627*0.5;
  if(sampleName=="WplusH130") crossSection = 1.209*0.0541*0.5;

  if(sampleName=="WminusH120") crossSection = 1.565*0.0698*0.5;
  if(sampleName=="WminusH125") crossSection = 1.373*0.0627*0.5;
  if(sampleName=="WminusH130") crossSection = 1.209*0.0541*0.5;

  if(sampleName=="ZH120") crossSection = 0.994*0.0698;
  if(sampleName=="ZH125") crossSection = 0.884*0.0627;
  if(sampleName=="ZH130") crossSection = 0.790*0.0541;

  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
  if(sampleName.find("ZZTo2L2Q")!=std::string::npos) crossSection = 3.22;
  if(sampleName.find("ZZTo4L")!=std::string::npos) crossSection = 1.212;
  if(sampleName.find("WZTo1L3Nu")!=std::string::npos) crossSection = 3.05;
  if(sampleName.find("WZJToLLLNu")!=std::string::npos) crossSection = 4.708;
  if(sampleName.find("WWTo1L1Nu2Q")!=std::string::npos) crossSection = 1.212;
  if(sampleName.find("WZTo1L1Nu2Q")!=std::string::npos) crossSection = 10.71;
  if(sampleName.find("VVTo2L2Nu")!=std::string::npos) crossSection = 11.95;
  if(sampleName.find("WZTo2L2Q")!=std::string::npos) crossSection = 5.595;

  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
  if(sampleName=="Wantitop") crossSection = 35.6;
  if(sampleName=="Wtop") crossSection = 35.6;
  if(sampleName=="t-channel_top") crossSection = 136.02;
  if(sampleName=="t-channel_antitop") crossSection = 80.95;
  if(sampleName=="EWKWMinus") crossSection = 20.25;
  if(sampleName=="EWKWPlus") crossSection = 25.62;
  if(sampleName=="EWKZ2JetsZToLL") crossSection = 3.987;
  if(sampleName=="EWKZ2JetsZToNuNu") crossSection = 10.01;


  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
  ///matching eff. = 0.00042
  if(sampleName=="QCD_MC") crossSection = 720648000*0.00042;

  float weight = crossSection*presEff/nEventsAnalysed;
  if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;

  std::cout<<"Sample name: "<<sampleName<<" ";
  std::cout<<"Xsection: "<<crossSection<<" [pb] "<<" ";
  std::cout<<"Events analyzed: "<<nEventsAnalysed<<" ";
  std::cout<<"Reco preselection efficiency: "<<recoPresEff<<" ";
  std::cout<<"Final weight: "<<weight<<std::endl;

  return weight;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
HTTHistograms::HTTHistograms(TDirectory *myDir){ AnalysisHistograms::init(myDir); }
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
TH1F *HTTHistograms::get1D_EWK2JetsSum(const std::string& name){

  std::vector<std::string> ewkSamples = {"EWKWMinus", "EWKWPlus", "EWKZ2JetsZToLL", "EWKZ2JetsZToNuNu"};
  return get1D_SumPattern_Histogram(name, "EWK2Jets", ewkSamples);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_TauMatchJetSum(const std::string& name, bool sumDecayModes, bool sumJetBins){

  std::vector<std::string> decayNames = {"MatchL", "MatchJ", "MatchT"};

  TString hName = name;
  TH1F *hSum = 0;

  if(sumDecayModes){
    ///Strip decay name if any.
    for(auto decayName:decayNames){hName.ReplaceAll(decayName.c_str(), "");}
    for(auto decayName:decayNames){
      TString hNameTmp = hName;
      if(hNameTmp.First("_")>0) hNameTmp.Replace(hNameTmp.First("_"),1,(decayName+"_").c_str()); else hNameTmp.Append(decayName.c_str());
      TH1F *hDecayMode = 0;
      if(sumJetBins) hDecayMode = get1D_VJetSum(hNameTmp.Data());
      else hDecayMode = get1DHistogram(hNameTmp.Data());
      if(!hSum && hDecayMode){
	hSum = (TH1F*)hDecayMode->Clone(name.c_str());
	hSum->Reset();
      }
      if(hDecayMode) hSum->Add(hDecayMode);
    }
  }
  else if(sumJetBins) hSum = get1D_VJetSum(name);

  if(hSum) hSum->SetName(name.c_str());
  return hSum;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_WJet_Histogram(const std::string& name){return get1D_VJetSum(name);}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_DYJet_Histogram(const std::string& name){

  bool sumDecayModes = true;
  bool sumJetBins = true;
  TH1F *histo = get1D_TauMatchJetSum(name, sumDecayModes, sumJetBins);

  return histo;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::getNormalised_NJet_Histogram(const std::string& hName){

  TH1F *hNJets = get1DHistogram(hName);
  if(!hNJets) return hNJets;

  std::string sampleName = "";
  TH1F *hNJetsStats = 0;
  if(hName.find("W")!=std::string::npos){
    sampleName = hName.substr(hName.find("W"));
    std::string selName = sampleName.substr(sampleName.find("_"));
    sampleName = sampleName.substr(0,sampleName.size()-selName.size());
    hNJetsStats = get1DHistogram("h1DStats"+sampleName);
  }
  if(hName.find("DY")!=std::string::npos){
    sampleName = hName.substr(hName.find("DY"));
    std::string selName = sampleName.substr(sampleName.find("_"));
    sampleName = sampleName.substr(0,sampleName.size()-selName.size());
    bool sumDecayModes = true;
    bool sumJetBins = false;
    hNJetsStats = get1D_TauMatchJetSum("h1DStats"+sampleName, sumDecayModes, sumJetBins);
  }

  if(!hNJetsStats) return hNJets;

  float recoPresEff = hNJetsStats->GetBinContent(3)/hNJetsStats->GetBinContent(2);
  int nEventsAnalysed = hNJetsStats->GetBinContent(1);

  if(sampleName.find("0Jets")!=std::string::npos ||
     sampleName.find("AllJets")!=std::string::npos){
    TString allJetsName = "h1DStats"+sampleName;
    if(sampleName.find("0Jets")!=std::string::npos) allJetsName.ReplaceAll("0Jets","AllJets");
    if(sampleName.find("AllJets")!=std::string::npos) allJetsName.ReplaceAll("AllJets","0Jets");
    TH1F *hAllJetsStats = 0;
    bool sumDecayModes = true;
    bool sumJetBins = false;
    if(hName.find("W")!=std::string::npos) hAllJetsStats = get1DHistogram(allJetsName.Data());
    if(hName.find("DY")!=std::string::npos) hAllJetsStats = get1D_TauMatchJetSum(allJetsName.Data(), sumDecayModes, sumJetBins);
    recoPresEff =  (hNJetsStats->GetBinContent(3) + hAllJetsStats->GetBinContent(3));
    recoPresEff /= (hNJetsStats->GetBinContent(2) + hAllJetsStats->GetBinContent(2));
    nEventsAnalysed = hNJetsStats->GetBinContent(1) + hAllJetsStats->GetBinContent(1);
  }

  if(hNJets) hNJets->Scale(recoPresEff/nEventsAnalysed);
  return hNJets;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1D_VJetSum(const std::string& name){

  TString hName = name;

  ///TEST
  //if(name.find("Jets")==std::string::npos) return get1DHistogram(name.c_str());
  //return getNormalised_NJet_Histogram(name.c_str());
  ////////////////////

  if(name.find("AllJets")!=std::string::npos) return getNormalised_NJet_Histogram(name.c_str());
  if(name.find("0Jets")!=std::string::npos) return getNormalised_NJet_Histogram(name.c_str());
  if(name.find("Jets")==std::string::npos) return get1DHistogram(name.c_str());

  hName.ReplaceAll("Jets","0Jets");
  TH1F *h0Jets = getNormalised_NJet_Histogram(hName.Data());

  hName = name;
  hName.ReplaceAll("Jets","1Jets");
  TH1F *h1Jets = getNormalised_NJet_Histogram(hName.Data());

  hName = name;
  hName.ReplaceAll("Jets","2Jets");
  TH1F *h2Jets = getNormalised_NJet_Histogram(hName.Data());

  hName = name;
  hName.ReplaceAll("Jets","3Jets");
  TH1F *h3Jets = getNormalised_NJet_Histogram(hName.Data());

  hName = name;
  hName.ReplaceAll("Jets","4Jets");
  TH1F *h4Jets = getNormalised_NJet_Histogram(hName.Data());

  TH1F *hJets = 0;
  if(h0Jets) hJets = (TH1F*)h0Jets->Clone(name.c_str());
  else if(h1Jets) hJets = (TH1F*)h1Jets->Clone(name.c_str());
  else if(h2Jets) hJets = (TH1F*)h2Jets->Clone(name.c_str());

  if(!hJets) return 0;
  if(!h1Jets && !h2Jets && !h3Jets && !h4Jets) return getNormalised_NJet_Histogram(name.c_str());

  std::vector<float> jetsLOSigma(5);
  if(name.find("W")!=std::string::npos) jetsLOSigma = {50380, 9644.5, 3144.5, 954.8, 485.6};
  if(name.find("DY")!=std::string::npos) jetsLOSigma = {4954.0, 1012.5, 332.8, 101.8, 54.8};

  hJets->Reset();
  if(h0Jets) hJets->Add(h0Jets, jetsLOSigma[0]/jetsLOSigma[0]);
  if(h1Jets) hJets->Add(h1Jets, jetsLOSigma[1]/jetsLOSigma[0]);
  if(h2Jets) hJets->Add(h2Jets, jetsLOSigma[2]/jetsLOSigma[0]);
  if(h3Jets) hJets->Add(h3Jets, jetsLOSigma[3]/jetsLOSigma[0]);
  if(h4Jets) hJets->Add(h4Jets, jetsLOSigma[4]/jetsLOSigma[0]);

  return hJets;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::get1D_VV_Histogram(const std::string& name, std::string tauMatchSuffix){

  std::vector<std::string> sampleNamesVV = {"ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};
  return get1D_SumPattern_Histogram(name, "DiBoson", sampleNamesVV, tauMatchSuffix);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::get1D_SumPattern_Histogram(const std::string& name, std::string pattern, std::vector<std::string> sampleNames, std::string tauMatchSuffix){

  TString hName = name;
  TH1F *hSum = 0;

  for(auto sampleName:sampleNames){
    TString hNameTmp = hName;
    hNameTmp.ReplaceAll(pattern,(sampleName+tauMatchSuffix).c_str());
    TH1F *histo = get1DHistogram(hNameTmp.Data());

    if(tauMatchSuffix == "" && !histo) histo = get1D_TauMatchJetSum(hNameTmp.Data(), true, false);

    if(!hSum && histo){
      hSum = (TH1F*)histo->Clone(name.c_str());
      hSum->Reset();
    }
    if(histo && hSum){
      float scale = getSampleNormalisation(sampleName+tauMatchSuffix);
      hSum->Add(histo, scale);
    }
  }

  hName.ReplaceAll(pattern, (pattern+tauMatchSuffix).c_str());
  if(hSum) hSum->SetName(hName.Data());

  return hSum;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::get1D_ST_Histogram(const std::string& name){

  std::vector<std::string> sampleNamesST = {"Wtop", "Wantitop","t-channel_top","t-channel_antitop"};
  return get1D_SumPattern_Histogram(name, "ST", sampleNamesST);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::get1D_TT_Histogram(const std::string& name, std::string tauMatchSuffix){

  std::vector<std::string> sampleNamesTT = {"TTbar"};
  return get1D_SumPattern_Histogram(name, "TTbar", sampleNamesTT, tauMatchSuffix);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::string HTTHistograms::getTemplateName(const std::string& name){

  std::string templateName = "";
  if(name.find("hProf")!=std::string::npos && name.find("VsMag")!=std::string::npos) templateName = "hProfVsMagTemplate";
  if(name.find("hProf")!=std::string::npos && name.find("VsPt")!=std::string::npos) templateName = "hProfVsPtTemplate";
  if(name.find("hProf")!=std::string::npos && name.find("VsCos")!=std::string::npos) templateName = "hProfVsCosTemplate";

  if(name.find("h1DNPV")!=std::string::npos) templateName = "h1DNPVTemplate";
  if(name.find("h1DNPU")!=std::string::npos) templateName = "h1DNPUTemplate";
  if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";
  if(name.find("h1DWideMass")!=std::string::npos) templateName = "h1DWideMassTemplate";
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
  if(name.find("h1DyTau")!=std::string::npos) templateName = "h1DyTauTemplate";
  if(name.find("h1DNPartons")!=std::string::npos) templateName = "h1DStatsTemplate";

  if(name.find("h2DVxPullVsNTrack")!=std::string::npos) templateName = "h2DVxPullVsNTrackTemplate";

  if(name.find("h1DUnRollTauPtMassVis")!=std::string::npos) templateName = "h1DUnRollTauPtMassVisTemplate";
  if(name.find("h2DRollTauPtMassVis")!=std::string::npos) templateName = "h2DRollTauPtMassVisTemplate";

  if(name.find("h1DUnRollHiggsPtMassSV")!=std::string::npos) templateName = "h1DUnRollHiggsPtMassSVTemplate";
  if(name.find("h2DRollHiggsPtMassSV")!=std::string::npos) templateName = "h2DRollHiggsPtMassSVTemplate";

  if(name.find("h1DUnRollMjjMassSV")!=std::string::npos) templateName = "h1DUnRollMjjMassSVTemplate";
  if(name.find("h2DRollMjjMassSV")!=std::string::npos) templateName = "h2DRollMjjMassSVTemplate";

  if(name.find("h1DUnRollMassSVPhiCP")!=std::string::npos) templateName = "h1DUnRollMassSVPhiCPTemplate";
  if(name.find("h2DRollMassSVPhiCP")!=std::string::npos) templateName = "h2DRollMassSVPhiCPTemplate";

  if(name.find("h1DUnRollMassSVYCP")!=std::string::npos) templateName = "h1DUnRollMassSVYCPTemplate";
  if(name.find("h2DRollMassSVYCP")!=std::string::npos) templateName = "h2DRollMassSVYCPTemplate";

  return templateName;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::defineHistograms(){

 using namespace std;

 if(!histosInitialized_){
   add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
   add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
   add1DHistogram("h1DNPUTemplate",";Number of PV; Events",600,0,60,file_);
   add1DHistogram("h1DMassTemplate",";mass [GeV/c^{2}]; Events",35,0,350,file_);
   add1DHistogram("h1DWideMassTemplate",";mass [GeV/c^{2}]; Events",25,0,1500,file_);
   add1DHistogram("h1DPtTemplate",";p_{T}; Events",20,0,100,file_);
   add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
   add1DHistogram("h1DPhiTemplate",";#phi; Events",12,0,2*M_PI,file_);
   add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",10,-1.0,1.0,file_);
   add1DHistogram("h1DCSVBtagTemplate",";CSV btag; Events",20,0,1,file_);
   add1DHistogram("h1DIsoTemplate",";Isolation; Events",20,0,0.3,file_);
   add1DHistogram("h1DIDTemplate",";ID; Events",20,0.8,1,file_);
   add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",11,-0.01,0.01,file_);
   add1DHistogram("h1DyTauTemplate",";yTau; Events",15,-1,1,file_);
   add1DHistogram("h1DnPCATemplate",";#hat{n}_{RECO}>; Events",10,0,0.015,file_);

   // 0jet: unroll in tau_Pt using bins of {30,35,40,45,50,55,300} and in visible mass using bins of {0,60,65,70,75,80,85,90,95,100,105,110,400}
   std::vector<double> massVisBins = {0,60,65,70,75,80,85,90,95,100,105,110,400};
   std::vector<double> tauPtBins = {30,35,40,45,50,55,300};
   addRollHistogram("h1DUnRollTauPtMassVisTemplate","Visible Mass vs tau pt; Events",massVisBins, tauPtBins, file_);
   //Boosted: unroll in Higgs_Pt using Higgs_Pt bins of {0,100,150,200,250,300,5000}, and in SV mass using bins of {0,80,90,100,110,120,130,140,150,160,300}
   std::vector<double> higgsPtBins = {0,100,150,200,250,300,5000};
   std::vector<double> svMassBins = {0,80,90,100,110,120,130,140,150,160,300};
   addRollHistogram("h1DUnRollHiggsPtMassSVTemplate","SV Mass vs Higgs Pt; Events",svMassBins, higgsPtBins, file_);
   //VBF: unroll in mjj using bins of {300,700,1100,1500,10000}, and in SV mass using bins of {0,95,115,135,155,400}
   vector<double> mjjBins = {300,700,1100,1500,10000};
   vector<double> svMassBinsVBF =  {0,95,115,135,155,400};
   addRollHistogram("h1DUnRollMjjMassSVTemplate","SV Mass vs Mjj; Events",svMassBinsVBF, mjjBins, file_);
   ///2D CP histograms
   vector<double> phiBins;
   for(unsigned int iBin=0;iBin<=12;++iBin) phiBins.push_back(iBin*2.0*M_PI/12);
   addRollHistogram("h1DUnRollMassSVPhiCPTemplate","#phi_{IP,IP} CP vs SV Mass; Events;#phi_{IP,IP} CP",phiBins, svMassBins, file_);
   addRollHistogram("h1DUnRollMassSVYCPTemplate","#phi_{IP,#rho} CP vs SV Mass; Events;#phi_{IP,#rho} CP", phiBins, svMassBins, file_);

   addProfile("hProfVsMagTemplate","",10,0,0.015,file_);
   addProfile("hProfVsPtTemplate","",20,15,55,file_);
   addProfile("hProfVsCosTemplate","",20,-1,1,file_);


 }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(int nRuns, float weight){

  AnalysisHistograms::finalizeHistograms();

  ////Code below tests W+n jets normalisation
  ///Samples split into jet multiplicity are compared to
  ///inclusive sample.
  /*
  TH1F *hW = get1D_WJet_Histogram("h1DNPartonsWJets");
  TH1F *hWAllJets = get1D_WJet_Histogram("h1DNPartonsWAllJets");
  hWAllJets->Print("all");
  hW->Print("all");
  hWAllJets->Divide(hW);
  hWAllJets->Print("all");
  return;
  */
  /*
  TH1F *hDY = get1D_VJetSum("h1DNPartonsDYMuTauJets");
  TH1F *hDYAllJets = get1D_VJetSum("h1DNPartonsDYMuTauAllJets");
  hDYAllJets->Print("all");
  hDY->Print("all");
  hDYAllJets->Divide(hDY);
  hDYAllJets->Print("all");
  */
  //TH1F *hDY =  get1D_VJetSum("h1DMassVisDYMuTauAllJets");
  //hDY->Add(get1D_VJetSum("h1DMassVisDYMuTau0Jets"));
  /*
  TH1F *hDYSum = get1D_DYJet_Histogram("h1DMassVisDYJets");
  hDYSum->Print("all");
  return;
  */
  //////////////
  ///Control regions plots
  ttScale = 1.0;

return;

  for(unsigned int iCategory = (int)HTTAnalyzer::jet0;
      iCategory<(int)HTTAnalyzer::boosted;++iCategory){

    plotCPhistograms(iCategory);

    wselOSCorrection =  std::pair<float,float>(1.0,0);
    wselSSCorrection =  std::pair<float,float>(1.0,0);

    wselOSCorrection = getWNormalisation(iCategory, "OS");
    wselSSCorrection = getWNormalisation(iCategory, "SS");

    //plotStack(iCategory, "MassVis");
    plotStack(iCategory, "UnRollTauPtMassVis");
    //plotStack(iCategory, "UnRollHiggsPtMassSV");
    //plotStack(iCategory, "UnRollMjjMassSV");
    //plotStack(iCategory, "UnRollMassSVPhiCP");
    //plotStack(iCategory, "UnRollMassSVYCP");
    //plotStack(iCategory, "MassTrans");

    continue;

    plotStack(iCategory, "PtMuon");
    plotStack(iCategory, "EtaMuon");
    plotStack(iCategory, "IsoMuon");

    plotStack(iCategory, "PtTau");
    plotStack(iCategory, "EtaTau");
    plotStack(iCategory, "IDTau");
    plotStack(iCategory, "StatsDecayMode");
    plotStack(iCategory, "PtTauLeadingTk");

    plotStack(iCategory, "PhiMuon");
    plotStack(iCategory, "PhiTau");

    plotStack(iCategory, "PtMET");
    plotStack(iCategory, "PtMuTauMET");

    plotStack(iCategory, "StatsNJ30");
    plotStack(iCategory, "PtLeadingJet");
    plotStack(iCategory, "EtaLeadingJet");
    plotStack(iCategory, "CSVBtagLeadingJet");
    plotStack(iCategory, "PtLeadingBJet");
    plotStack(iCategory, "EtaLeadingBJet");
    plotStack(iCategory, "WideMass2J");

    plotStack(iCategory, "nPCAMuon");
    plotStack(iCategory, "nPCATau");
    plotStack(iCategory, "Phi-nVectors");
    plotStack(iCategory, "Phi-nVecIP");
    plotStack(iCategory, "NPV");
    }

    ///Make the ststematic effect histos.

    for(unsigned int iSystEffect = (unsigned int)sysEffects::NOMINAL_SVFIT;
	iSystEffect<(unsigned int)sysEffects::DUMMY;++iSystEffect){

      for(unsigned int iCategory = (int)HTTAnalyzer::jet0;
	  iCategory<(int)HTTAnalyzer::boosted;++iCategory){

	wselOSCorrection =  std::pair<float,float>(1.0,0);
	wselSSCorrection =  std::pair<float,float>(1.0,0);

	wselOSCorrection = getWNormalisation(iCategory, "OS", iSystEffect);
	wselSSCorrection = getWNormalisation(iCategory, "SS", iSystEffect);

	plotStack(iCategory, "MassSV", "OS", iSystEffect);
	plotStack(iCategory, "UnRollTauPtMassVis", "OS", iSystEffect);
	plotStack(iCategory, "UnRollHiggsPtMassSV", "OS", iSystEffect);
	plotStack(iCategory, "UnRollMjjMassSV", "OS", iSystEffect);
	plotStack(iCategory, "UnRollMassSVPhiCP", "OS", iSystEffect);
	plotStack(iCategory, "UnRollMassSVYCP", "OS", iSystEffect);
      }
    }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotCPhistograms(unsigned int iCategory){

  std::string hNameSuffix = "";
  plot_HAZ_Histograms("Phi-nVecIP-yTauNeg",hNameSuffix+"_GenNoOfflineSel");
  plot_HAZ_Histograms("Phi-nVecIP-yTauPos",hNameSuffix+"_GenNoOfflineSel");
  plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_GenNoOfflineSel");
  plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_GenNoOfflineSel");

  hNameSuffix =  "_OS_"+std::to_string(iCategory);

  //plot_HAZ_Histograms("Phi-nVecIP-yTauNeg",hNameSuffix+"_Gen");
  //plot_HAZ_Histograms("Phi-nVecIP-yTauPos",hNameSuffix+"_Gen");
  plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_Gen");
  plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_Gen");

  plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_RefitPV");
  plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_RefitPV");

  plotnPCA("ggH125"+hNameSuffix);
  plotnPCA("A"+hNameSuffix);
  plotnPCA("DYJets"+hNameSuffix);
  plotnPCA("WJets"+hNameSuffix);

  plotPhiDecayPlanes("Phi-nVectorsggH125"+hNameSuffix);
  plotPhiDecayPlanes("Phi-nVectorsA"+hNameSuffix);
  plotPhiDecayPlanes("Phi-nVectorsDYJets"+hNameSuffix);

  plotPhiDecayPlanes("Phi-nVecIPggH125"+hNameSuffix);
  plotPhiDecayPlanes("Phi-nVecIPA"+hNameSuffix);
  plotPhiDecayPlanes("Phi-nVecIPDYJets"+hNameSuffix);

  plotProfiles("hProfRecoVsMagGen","ggH125"+hNameSuffix);
  plotProfiles("hProfPhiVsMag","ggH125"+hNameSuffix);

  plotVerticesPulls("h1DVxPullXggH125"+hNameSuffix);
  plotVerticesPulls("h1DVxPullYggH125"+hNameSuffix);
  plotVerticesPulls("h1DVxPullZggH125"+hNameSuffix);

  plotPhiDecayPlanes("Phi-nVecIP-yTauPosggH125"+hNameSuffix);
  plotPhiDecayPlanes("Phi-nVecIP-yTauNegggH125"+hNameSuffix);

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotnPCA(const std::string & type){

  TH1F* h1DTau = 0;
  TH1F* h1DMuon = 0;

  if(type.find("DYJets")!=std::string::npos){
    std::string typeztt="DYJetsMatchT"+type.substr(type.find("DYJets")+6);
    h1DTau = get1D_TauMatchJetSum(("h1DnPCATau"+typeztt).c_str(),false,true);
    h1DMuon = get1D_TauMatchJetSum(("h1DnPCAMuon"+typeztt).c_str(),false,true);
  }
  else if(type.find("WJets")!=std::string::npos){
    h1DTau = get1D_WJet_Histogram("h1DnPCATau"+type);
    h1DMuon = get1D_WJet_Histogram("h1DnPCAMuon"+type);
  }
  else{
    h1DTau = get1DHistogram("h1DnPCATau"+type);
    h1DMuon = get1DHistogram("h1DnPCAMuon"+type);
  }

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
   c->SetLeftMargin(0.15);

  TLegend l(0.15,0.7,0.3,0.85,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  if(hName.find("2D")!=std::string::npos){
    TProfile* hProfile_AOD = this->get2DHistogram((hName+"_AODPV").c_str())->ProfileX();
    TProfile* hProfile_Refit = this->get2DHistogram((hName+"_RefitPV").c_str())->ProfileX();

    if(!hProfile_AOD || !hProfile_Refit) return;

    hProfile_AOD->SetLineWidth(3);
    hProfile_Refit->SetLineWidth(3);

    hProfile_AOD->SetLineColor(1);
    hProfile_Refit->SetLineColor(4);

    hProfile_AOD->Draw();
    hProfile_Refit->Draw("same");

    hProfile_AOD->GetXaxis()->SetTitle("number of tracks in PV");
    hProfile_AOD->GetYaxis()->SetTitle("#sigma(PV^{RECO} - PV^{GEN})");
    hProfile_AOD->GetYaxis()->SetTitleOffset(2.2);

    float min = hProfile_AOD->GetMinimum();
    if(hProfile_Refit->GetMinimum()<min) min = hProfile_Refit->GetMinimum();
    hProfile_AOD->SetMinimum(0.95*min);

    l.AddEntry(hProfile_AOD,"from AOD");
    l.AddEntry(hProfile_Refit,"#splitline{refitted from}{mAOD, with BS}");
    l.Draw();

    c->Print(TString::Format("fig_png/%s.png",hName.c_str()).Data());
    return;
  }

  TH1F* h1D_AOD = this->get1DHistogram((hName+"_AODPV").c_str());
  TH1F* h1D_Refit = this->get1DHistogram((hName+"_RefitPV").c_str());

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

  TProfile* h1DAOD = this->getProfile(hName+sysType+"_AODPV");
  //TProfile* h1DGen = this->getProfile(hName+sysType+"_GenNoOfflineSel");
  TProfile* h1DGen = this->getProfile(hName+sysType+"_GenPV");
  TProfile* h1DRefit = this->getProfile(hName+sysType+"_RefitPV");

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
    if(hName.find("PtVsMag")!=std::string::npos){
      h1DGen->SetYTitle("p_{T}^{leading tk.}");
      h1DGen->SetXTitle("<|n_{GEN}|>");
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
    if(hName.find("PhiVsMag")!=std::string::npos) l.Draw();

    c.Print(TString::Format("fig_png/%s.png",(hName+sysType).c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotPhiDecayPlanes(const std::string & name){

  TCanvas aCanvas(TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  TString::Format("PhiDecayPlanes_%s",name.c_str()),
		  460,500);

  TLegend l(0.15,0.15,0.35,0.4,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TString hName = "h1D"+name+"_RefitPV";
  TH1F* h1DRefitPV = 0;
  if(name.find("DYJets")!=std::string::npos) {
    std::string n2=name;
    n2.insert(n2.find("DYJets")+6,"MatchT");
    hName = "h1D"+n2+"_RefitPV";
    h1DRefitPV = get1D_TauMatchJetSum(hName.Data(),false,true);
  }
  else h1DRefitPV = get1DHistogram(hName.Data());

  hName = "h1D"+name+"_AODPV";
  TH1F* h1DAODPV = 0;
  if(name.find("DYJets")!=std::string::npos) {
    std::string n2=name;
    n2.insert(n2.find("DYJets")+6,"MatchT");
    hName = "h1D"+n2+"_AODPV";
    h1DAODPV = get1D_TauMatchJetSum(hName.Data(),false,true);
  }
  else h1DAODPV = get1DHistogram(hName.Data());

  hName = "h1D"+name+"_GenPV";
  TH1F* h1DGenPV = 0;
  if(name.find("DYJets")!=std::string::npos) {
    std::string n2=name;
    n2.insert(n2.find("DYJets")+6,"MatchT");
    hName = "h1D"+n2+"_GenPV";
    h1DGenPV = get1D_TauMatchJetSum(hName.Data(),false,true);
  }
  else h1DGenPV = get1DHistogram(hName.Data());

  //hName = "h1D"+name+"_GenNoOfflineSel";
  TH1F* h1DGen = 0;
  //if(name.find("DYJets")!=std::string::npos) {
  //  std::string n2=name;
  //  n2.insert(n2.find("DYJets")+6,"MatchT");
  //  hName = "h1D"+n2+"_GenNoOfflineSel";
  //  h1DGen = get1D_TauMatchJetSum(hName.Data(),false,true);
  //}
  //else h1DGen = get1DHistogram(hName.Data());

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

    if(h1DGenPV && h1DGenPV->GetMaximum()> h1DRefitPV->GetMaximum()) h1DRefitPV->SetMaximum(1.02*h1DGenPV->GetMaximum());
    h1DRefitPV->SetMinimum(0);
    h1DRefitPV->Draw("HISTO");

    h1DGenPV->Print();

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
      l.AddEntry(h1DGen,"PCA gen. particles, no sel.");
    }

    h1DRefitPV->Draw("HISTO same");
    l.Draw();
    aCanvas.Print(TString::Format("fig_png/%s.png",name.c_str()).Data());
  }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plot_HAZ_Histograms(const std::string & hName,
					const std::string & sysType){

  TCanvas* c = new TCanvas(TString::Format("%s_%s",hName.c_str(), sysType.c_str()),
			   TString::Format("%s_%s",hName.c_str(), sysType.c_str()),
			   460,500);

  TLegend l(0.35,0.15,0.55,0.35,NULL,"brNDC");
  l.SetTextSize(0.05);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  l.SetFillColor(10);

  TLatex aLatex(0,0,"");

  TString name = "h1D"+hName+"ggH125"+sysType;
  TH1F* h_h = this->get1DHistogram(name.Data());
  name = "h1D"+hName+"A"+sysType;
  TH1F* h_A = this->get1DHistogram(name.Data());
  name = "h1D"+hName+"DYJetsMatchT"+sysType;
  TH1F* h_Z = get1D_TauMatchJetSum(name.Data(),false,true);

  if(!h_h || !h_A || !h_Z) return;

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
  l.AddEntry(h_h,"SM ggH125");
  l.AddEntry(h_A,"MSSM A");
  l.AddEntry(h_Z,"SM Z");
  l.Draw();
  //aLatex.DrawLatex(0.05,0.02,sysType.c_str());
  ///
  c->Print(TString::Format("fig_png/%s_h_A_Z_%s.png",hName.c_str(), sysType.c_str()).Data());
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
THStack*  HTTHistograms::plotStack(unsigned int iCategory,
				   std::string varName, std::string selName,
				   unsigned int iSystEffect){

  std::cout<<"--- Drawing THStack for variable: "<<varName
	   <<" category: "<<iCategory<<std::endl;

  std::string hName = "h1D"+varName;
  std::string categoryName = HTTAnalyzer::categoryName(iCategory);
  std::string systEffectName = HTTAnalyzer::systEffectName(iSystEffect);

  if(systEffectName.find("CAT")!=std::string::npos){
    systEffectName.replace(systEffectName.find("CAT"),3,categoryName);
  }

  std::string hNameSuffix =  "_"+selName+"_"+std::to_string(iCategory);
  hNameSuffix+=systEffectName;

  TH1F *hggHiggs120 = get1DHistogram((hName+"ggH120"+hNameSuffix).c_str());
  TH1F *hggHiggs125 = get1DHistogram((hName+"ggH125"+hNameSuffix).c_str());
  TH1F *hggHiggs130 = get1DHistogram((hName+"ggH130"+hNameSuffix).c_str());

  TH1F *hqqHiggs120 = get1DHistogram((hName+"qqH120"+hNameSuffix).c_str());
  TH1F *hqqHiggs125 = get1DHistogram((hName+"qqH125"+hNameSuffix).c_str());
  TH1F *hqqHiggs130 = get1DHistogram((hName+"qqH130"+hNameSuffix).c_str());

  TH1F *hZHiggs120 = get1DHistogram((hName+"ZH120"+hNameSuffix).c_str());
  TH1F *hZHiggs125 = get1DHistogram((hName+"ZH125"+hNameSuffix).c_str());
  TH1F *hZHiggs130 = get1DHistogram((hName+"ZH130"+hNameSuffix).c_str());

  TH1F *hWplusHiggs120 = get1DHistogram((hName+"WplusH120"+hNameSuffix).c_str());
  TH1F *hWplusHiggs125 = get1DHistogram((hName+"WplusH125"+hNameSuffix).c_str());
  TH1F *hWplusHiggs130 = get1DHistogram((hName+"WplusH130"+hNameSuffix).c_str());

  TH1F *hWminusHiggs120 = get1DHistogram((hName+"WminusH120"+hNameSuffix).c_str());
  TH1F *hWminusHiggs125 = get1DHistogram((hName+"WminusH125"+hNameSuffix).c_str());
  TH1F *hWminusHiggs130 = get1DHistogram((hName+"WminusH130"+hNameSuffix).c_str());

  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix).c_str());
  TH1F *hTTbarJ = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str(),"MatchJ");
  TH1F *hTTbarT = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str(),"MatchT");
  TH1F *hST = get1D_ST_Histogram((hName+"ST"+hNameSuffix).c_str());
  TH1F *hVVJ = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str(), "MatchJ");
  TH1F *hVVT = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str(), "MatchT");
  TH1F *hDYJetsLowM = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix).c_str());

  bool sumDecayModes = false;
  bool sumJetBins = true;
  TH1F *hDYJetsZJ = get1D_TauMatchJetSum((hName+"DYJetsMatchJ"+hNameSuffix).c_str(), sumDecayModes, sumJetBins);
  TH1F *hDYJetsZL = get1D_TauMatchJetSum((hName+"DYJetsMatchL"+hNameSuffix).c_str(), sumDecayModes, sumJetBins);
  TH1F *hDYJetsZTT = get1D_TauMatchJetSum((hName+"DYJetsMatchT"+hNameSuffix).c_str(), sumDecayModes, sumJetBins);

  TH1F *hEWK2Jets = get1D_EWK2JetsSum(hName+"EWK2Jets"+hNameSuffix);

  TH1F *hSoup = get1DHistogram((hName+"Data"+hNameSuffix).c_str(),true);
  pair<float,float> qcdOStoSS = getQCDOStoSS(iCategory, wselOSCorrection, wselSSCorrection, iSystEffect);
  TH1F *hQCD = (TH1F*)getQCDbackground(iCategory,varName,wselOSCorrection, wselSSCorrection, iSystEffect);
  TH1F *hQCD_MC =  get1DHistogram((hName+"QCD_MC"+hNameSuffix).c_str());

  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.
  if(!hSoup) return 0;

  TH1F *hEmpty = (TH1F*)hSoup->Clone("hEmpty");
  hEmpty->Reset();
  if(!hQCD) hQCD = (TH1F*)hEmpty->Clone((hName+"QCDEstimate_"+std::to_string(iCategory)).c_str());
  if(!hQCD_MC) hQCD_MC = (TH1F*)hEmpty->Clone((hName+"QCD_MC"+hNameSuffix).c_str());
  if(!hWJets) hWJets = (TH1F*)hEmpty->Clone((hName+"WJets"+hNameSuffix).c_str());
  if(!hDYJetsLowM) hDYJetsLowM = (TH1F*)hEmpty->Clone((hName+"DYLowM"+hNameSuffix).c_str());
  if(!hDYJetsZJ) hDYJetsZJ = (TH1F*)hEmpty->Clone((hName+"DYJetsMatchJ"+hNameSuffix).c_str());
  if(!hDYJetsZL) hDYJetsZL = (TH1F*)hEmpty->Clone((hName+"DYJetsMatchL"+hNameSuffix).c_str());
  if(!hDYJetsZTT) hDYJetsZTT = (TH1F*)hEmpty->Clone((hName+"DYJetsMatchT"+hNameSuffix).c_str());
  if(!hTTbarJ) hTTbarJ = (TH1F*)hEmpty->Clone((hName+"TTbarMatchJ"+hNameSuffix).c_str());
  if(!hTTbarT) hTTbarT = (TH1F*)hEmpty->Clone((hName+"TTbarMatchT"+hNameSuffix).c_str());
  if(!hST) hST = (TH1F*)hEmpty->Clone((hName+"ST"+hNameSuffix).c_str());
  if(!hVVJ) hVVJ = (TH1F*)hEmpty->Clone((hName+"DiBosonMatchJ"+hNameSuffix).c_str());
  if(!hVVT) hVVT = (TH1F*)hEmpty->Clone((hName+"DiBosonMatchT"+hNameSuffix).c_str());
  if(!hEWK2Jets) hEWK2Jets = (TH1F*)hEmpty->Clone((hName+"EWK2Jets"+hNameSuffix).c_str());

  if(!hggHiggs120) hggHiggs120 = (TH1F*)hEmpty->Clone((hName+"ggH120"+hNameSuffix).c_str());
  if(!hggHiggs125) hggHiggs125 = (TH1F*)hEmpty->Clone((hName+"ggH125"+hNameSuffix).c_str());
  if(!hggHiggs130) hggHiggs130 = (TH1F*)hEmpty->Clone((hName+"ggH130"+hNameSuffix).c_str());

  if(!hqqHiggs120) hqqHiggs120 = (TH1F*)hEmpty->Clone((hName+"qqH120"+hNameSuffix).c_str());
  if(!hqqHiggs125) hqqHiggs125 = (TH1F*)hEmpty->Clone((hName+"qqH125"+hNameSuffix).c_str());
  if(!hqqHiggs130) hqqHiggs130 = (TH1F*)hEmpty->Clone((hName+"qqH130"+hNameSuffix).c_str());

  if(!hZHiggs120) hZHiggs120 = (TH1F*)hEmpty->Clone((hName+"ZH120"+hNameSuffix).c_str());
  if(!hZHiggs125) hZHiggs125 = (TH1F*)hEmpty->Clone((hName+"ZH125"+hNameSuffix).c_str());
  if(!hZHiggs130) hZHiggs130 = (TH1F*)hEmpty->Clone((hName+"ZH130"+hNameSuffix).c_str());

  if(!hWplusHiggs120) hWplusHiggs120 = (TH1F*)hEmpty->Clone((hName+"WplusH120"+hNameSuffix).c_str());
  if(!hWplusHiggs125) hWplusHiggs125 = (TH1F*)hEmpty->Clone((hName+"WplusH125"+hNameSuffix).c_str());
  if(!hWplusHiggs130) hWplusHiggs130 = (TH1F*)hEmpty->Clone((hName+"WplusH130"+hNameSuffix).c_str());

  if(!hWminusHiggs120) hWminusHiggs120 = (TH1F*)hEmpty->Clone((hName+"WminusH120"+hNameSuffix).c_str());
  if(!hWminusHiggs125) hWminusHiggs125 = (TH1F*)hEmpty->Clone((hName+"WminusH125"+hNameSuffix).c_str());
  if(!hWminusHiggs130) hWminusHiggs130 = (TH1F*)hEmpty->Clone((hName+"WminusH130"+hNameSuffix).c_str());

  hggHiggs125->Print();

  ///Set histograms directory, so the histograms are saved
  if(hQCD) hQCD->SetDirectory(hSoup->GetDirectory());
  if(hQCD_MC) hQCD_MC->SetDirectory(hSoup->GetDirectory());
  if(hWJets) hWJets->SetDirectory(hSoup->GetDirectory());
  if(hDYJetsLowM) hDYJetsLowM->SetDirectory(hSoup->GetDirectory());
  if(hDYJetsZJ) hDYJetsZJ->SetDirectory(hSoup->GetDirectory());
  if(hDYJetsZL) hDYJetsZL->SetDirectory(hSoup->GetDirectory());
  if(hDYJetsZTT) hDYJetsZTT->SetDirectory(hSoup->GetDirectory());
  if(hTTbarJ) hTTbarJ->SetDirectory(hSoup->GetDirectory());
  if(hTTbarT) hTTbarT->SetDirectory(hSoup->GetDirectory());
  if(hST) hST->SetDirectory(hSoup->GetDirectory());
  if(hVVJ) hVVJ->SetDirectory(hSoup->GetDirectory());
  if(hVVT) hVVT->SetDirectory(hSoup->GetDirectory());
  if(hEWK2Jets) hEWK2Jets->SetDirectory(hSoup->GetDirectory());

  if(hggHiggs120) hggHiggs120->SetDirectory(hSoup->GetDirectory());
  if(hggHiggs125) hggHiggs125->SetDirectory(hSoup->GetDirectory());
  if(hggHiggs130) hggHiggs130->SetDirectory(hSoup->GetDirectory());

  if(hqqHiggs120) hqqHiggs120->SetDirectory(hSoup->GetDirectory());
  if(hqqHiggs125) hqqHiggs125->SetDirectory(hSoup->GetDirectory());
  if(hqqHiggs130) hqqHiggs130->SetDirectory(hSoup->GetDirectory());

  if(hZHiggs120) hZHiggs120->SetDirectory(hSoup->GetDirectory());
  if(hZHiggs125) hZHiggs125->SetDirectory(hSoup->GetDirectory());
  if(hZHiggs130) hZHiggs130->SetDirectory(hSoup->GetDirectory());

  if(hWplusHiggs120) hWplusHiggs120->SetDirectory(hSoup->GetDirectory());
  if(hWplusHiggs125) hWplusHiggs125->SetDirectory(hSoup->GetDirectory());
  if(hWplusHiggs130) hWplusHiggs130->SetDirectory(hSoup->GetDirectory());

  if(hWminusHiggs120) hWminusHiggs120->SetDirectory(hSoup->GetDirectory());
  if(hWminusHiggs125) hWminusHiggs125->SetDirectory(hSoup->GetDirectory());
  if(hWminusHiggs130) hWminusHiggs130->SetDirectory(hSoup->GetDirectory());


  TH1F *hHiggs = (TH1F*)hggHiggs125->Clone("hHiggs");
  hHiggs->Reset();

  float lumi = getLumi();
  std::string sampleName = "WJets";
  std::string WselType = "wselOS";
  pair<float,float> dataToMCScale = wselOSCorrection;

  if(hNameSuffix.find("SS")!=std::string::npos){
    WselType = "wselSS";
    dataToMCScale = wselSSCorrection;
  }
  float weight = getSampleNormalisation(sampleName);

  /////
  float scale = weight*lumi*dataToMCScale.first;
  hWJets->Scale(scale);

  sampleName = "DYLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJetsMatchJ";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsZJ->Scale(scale);

  sampleName = "DYJetsMatchL";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsZL->Scale(scale);

  sampleName = "DYJetsMatchT";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsZTT->Scale(scale);

  //TT samples scaled for preselection during stiching step
  sampleName = "TTbarMatchJ";
  scale=lumi;
  hTTbarJ->Scale(scale);

  //TT samples scaled for preselection during stiching step
  sampleName = "TTbarMatchT";
  scale=lumi;
  hTTbarT->Scale(scale);

  sampleName = "QCD_MC";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hQCD_MC->Scale(scale);

  //single top samples scaled for preselection during stiching step
  sampleName = "ST";
  scale = lumi;
  hST->Scale(scale);

  //VV samples scaled for preselection during stiching step
  sampleName = "DiBosonMatchJ";
  scale = lumi;
  hVVJ->Scale(scale);

  //VV samples scaled for preselection during stiching step
  sampleName = "DiBosonMatchT";
  scale = lumi;
  hVVT->Scale(scale);

  ///EWK 2Jets scaled for preselection during stiching step
  sampleName = "EWK2Jets";
  scale = lumi;
  hEWK2Jets->Scale(scale);

  sampleName = "ggH120";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hggHiggs120->Scale(scale);

  sampleName = "qqH120";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hqqHiggs120->Scale(scale);

  sampleName = "ggH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hggHiggs125->Scale(scale);

  sampleName = "qqH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hqqHiggs125->Scale(scale);

  sampleName = "ggH130";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hggHiggs130->Scale(scale);

  sampleName = "qqH130";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hqqHiggs130->Scale(scale);

  sampleName = "ZH120";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hZHiggs120->Scale(scale);

  sampleName = "ZH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hZHiggs125->Scale(scale);

  sampleName = "ZH130";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hZHiggs130->Scale(scale);

  sampleName = "WminusH120";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWminusHiggs120->Scale(scale);

  sampleName = "WminusH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWminusHiggs125->Scale(scale);

  sampleName = "WminusH130";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWminusHiggs130->Scale(scale);

  sampleName = "WplusH120";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWplusHiggs120->Scale(scale);

  sampleName = "WplusH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWplusHiggs125->Scale(scale);

  sampleName = "WplusH130";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hWplusHiggs130->Scale(scale);
  
  //join Wplus and Wminus processes////////////////////
  TH1F *hWHiggs120 = hWplusHiggs120;
  hWHiggs120->Add(hWminusHiggs120);
  hWHiggs120->SetName((hName+"WH120"+hNameSuffix).c_str());
  TH1F *hWHiggs125 = hWplusHiggs125;
  hWHiggs125->Add(hWminusHiggs125);
  hWHiggs125->SetName((hName+"WH125"+hNameSuffix).c_str());
  TH1F *hWHiggs130 = hWplusHiggs130;
  hWHiggs130->Add(hWminusHiggs130);
  hWHiggs130->SetName((hName+"WH130"+hNameSuffix).c_str());
  hWHiggs120->SetDirectory(hSoup->GetDirectory());
  hWHiggs125->SetDirectory(hSoup->GetDirectory());
  hWHiggs130->SetDirectory(hSoup->GetDirectory());
  /////////////////////////////////////////////////////

  hHiggs->Add(hggHiggs125);
  hHiggs->Add(hqqHiggs125);
  hHiggs->Add(hZHiggs125);
  hHiggs->Add(hWplusHiggs125);
  hHiggs->Add(hWminusHiggs125);
  //////////////////////////////////////////////////////
  hSoup->SetLineColor(1);
  hSoup->SetFillColor(1);
  hSoup->SetMarkerStyle(20);

  hWJets->SetFillColor(kRed+2);
  hTTbarJ->SetFillColor(kBlue+2);
  hTTbarT->SetFillColor(kBlue+10);
  hST->SetFillColor(kYellow-10);
  hVVJ->SetFillColor(kRed-10);
  hVVT->SetFillColor(kBlue-10);
  hDYJetsZL->SetFillColor(kOrange-3);
  hDYJetsZJ->SetFillColor(kOrange-6);
  hDYJetsZTT->SetFillColor(kOrange-9);
  hDYJetsLowM->SetFillColor(kOrange-7);
  hQCD->SetFillColor(kMagenta-10);
  hHiggs->SetFillColor(kCyan+4);

  hSoup->SetLineWidth(1);

  THStack *hs = new THStack("hs","Stacked histograms");
  /////////
  hs->Add(hHiggs,"hist");

  hs->Add(hEWK2Jets,"hist");
  hs->Add(hQCD,"hist");
  hs->Add(hTTbarJ,"hist");
  hs->Add(hTTbarT,"hist");
  hs->Add(hST,"hist");
  hs->Add(hVVJ,"hist");
  hs->Add(hVVT,"hist");
  hs->Add(hWJets,"hist");
  hs->Add(hDYJetsLowM,"hist");
  hs->Add(hDYJetsZJ,"hist");
  hs->Add(hDYJetsZL,"hist");
  hs->Add(hDYJetsZTT,"hist");
  ////////
  TH1F *hMCSum = (TH1F*)hWJets->Clone("hMCSum");
  hMCSum->Reset();
  hMCSum->Add(hDYJetsLowM);
  hMCSum->Add(hDYJetsZTT);
  hMCSum->Add(hDYJetsZL);
  hMCSum->Add(hDYJetsZJ);
  hMCSum->Add(hWJets);
  hMCSum->Add(hTTbarJ);
  hMCSum->Add(hTTbarT);
  hMCSum->Add(hST);
  hMCSum->Add(hVVJ);
  hMCSum->Add(hVVT);
  hMCSum->Add(hQCD);
  hMCSum->Add(hEWK2Jets);

  hNameSuffix =  "_"+selName+"_"+std::to_string(iCategory);
  hNameSuffix+="_"+categoryName;
  hNameSuffix+=systEffectName;

  std::stringstream outputStream;

  outputStream<<"Event count summary for selecion name: "<<hNameSuffix<<std::endl;
  outputStream<<"Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC (no Higgs): "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC TTbarJ: "<<hTTbarJ->Integral(0,hTTbarJ->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC TTbarT: "<<hTTbarT->Integral(0,hTTbarT->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC single T: "<<hST->Integral(0,hST->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC DiBosonJ: "<<hVVJ->Integral(0,hVVJ->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC DiBosonT: "<<hVVT->Integral(0,hVVT->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC ZTT: "<<hDYJetsZTT->Integral(0,hDYJetsZTT->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC ZL: "<<hDYJetsZL->Integral(0,hDYJetsZL->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC ZJ: "<<hDYJetsZJ->Integral(0,hDYJetsZJ->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC Z->ll(m<50): "<<hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1)<<std::endl;
  outputStream<<"MC H(125)->tau tau: "<<hHiggs->Integral(0,hHiggs->GetNbinsX()+1)<<std::endl;
  outputStream<<"qqH(125)->tau tau: "<<hqqHiggs125->Integral(0,hqqHiggs125->GetNbinsX()+1)<<std::endl;
  outputStream<<"ggH(125)->tau tau: "<<hggHiggs125->Integral(0,hggHiggs125->GetNbinsX()+1)<<std::endl;
  outputStream<<"ZH(125)->tau tau: "<<hZHiggs125->Integral(0,hZHiggs125->GetNbinsX()+1)<<std::endl;
  outputStream<<"WplusH(125)->tau tau: "<<hWplusHiggs125->Integral(0,hWplusHiggs125->GetNbinsX()+1)<<std::endl;
  outputStream<<"WminusH(125)->tau tau: "<<hWminusHiggs125->Integral(0,hWminusHiggs125->GetNbinsX()+1)<<std::endl;
  outputStream<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl;
  outputStream<<"QCD_MC: "<<hQCD_MC->Integral(0,hQCD_MC->GetNbinsX()+1)<<std::endl;
  outputStream<<"EWK 2Jets: "<<hEWK2Jets->Integral(0,hEWK2Jets->GetNbinsX()+1)<<std::endl;
  outputStream<<"Correction factors:"<<std::endl;
  outputStream<<"QCD SS to OS: "<<qcdOStoSS.first<<" +- "<<qcdOStoSS.second<<std::endl;
  outputStream<<"W MC to DATA: "<<dataToMCScale.first<<" +- "<<dataToMCScale.second<<std::endl;
  outputStream<<"----------------------------------------"<<std::endl;

  std::cout<<outputStream.str();
  ofstream eventCountFile("eventCount.txt",ios::out | ios::app);
  eventCountFile<<"HTTHistograms compilation time: "<<__TIMESTAMP__<<std::endl;
  eventCountFile<<outputStream.str();
  eventCountFile.close();

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

  if(!hNameSuffix.size()) hNameSuffix = "baseline";
  hs->SetTitle(("Variable: "+varName+" selection: "+hNameSuffix).c_str());
  hs->SetMaximum(4400);
  hs->Draw("hist");
  hs->GetXaxis()->SetTitle(varName.c_str());
  hs->GetYaxis()->SetTitleOffset(1.4);
  hMCSum->SetFillColor(5);
  /////////
  float highEnd = 170;
  float lowEnd = -150;

  if(varName.find("Phi-")!=std::string::npos) lowEnd = 0;

  int binHigh = hs->GetXaxis()->FindBin(highEnd);
  int binLow = hs->GetXaxis()->FindBin(lowEnd);

  if(hs->GetXaxis()->GetXmax()<highEnd || hs->GetXaxis()->GetXmax()>300) binHigh = hs->GetXaxis()->GetNbins();
  if(hs->GetXaxis()->GetXmin()>lowEnd) lowEnd = 1;

  if(varName.find("Roll")==std::string::npos) hs->GetXaxis()->SetRange(binLow,binHigh);
  highEnd =  hs->GetXaxis()->GetBinUpEdge(binHigh);

  char yTitle[200];
  sprintf(yTitle,"Events/%2.1f",hSoup->GetXaxis()->GetBinWidth(1));
  hs->GetYaxis()->SetTitle(yTitle);

  float max = hs->GetMaximum();
  if(hSoup->GetMaximum()>max) max = hSoup->GetMaximum();

  hs->GetHistogram()->SetTitleOffset(1.0);
  hs->SetMaximum(1.1*max);
  hs->SetMinimum(0.1);

  hSoup->DrawCopy("same");

  TLegend *leg = new TLegend(0.79,0.12,0.99,0.82,NULL,"brNDC");
  setupLegend(leg);
  leg->AddEntry(hSoup,"Data","lep");
  leg->AddEntry(hDYJetsZTT,"Z#rightarrow#mu#tau_{h}","f");
  leg->AddEntry(hDYJetsZL,"Z#rightarrow#mu#mu","f");
  leg->AddEntry(hDYJetsZJ,"Z, j#rightarrow#tau_{h}","f");
  leg->AddEntry(hDYJetsLowM,"Z#rightarrow#it{ll}(m<50)","f");
  leg->AddEntry(hWJets,"W#rightarrow#it{l}#nu","f");
  leg->AddEntry(hTTbarJ,"t#bar{t}, j#rightarrow#tau_{h}","f");
  leg->AddEntry(hTTbarT,"t#bar{t}, #tau_{h}","f");
  leg->AddEntry(hST,"single-t","f");
  leg->AddEntry(hVVJ,"VV, j#rightarrow#tau_{h}","f");
  leg->AddEntry(hVVT,"VV, tau_{h}","f");
  leg->AddEntry(hQCD,"QCD","f");
  leg->AddEntry(hEWK2Jets,"EWK","f");
  leg->AddEntry(hHiggs,"H(125)#rightarrow#tau#tau","f");
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

  hSoup = (TH1F*)hSoup->Clone("hDataMCRatio");
  hSoup->SetDirectory(0);
  hSoup->GetXaxis()->SetRange(binLow,binHigh);
  hSoup->SetTitle("");
  hSoup->SetXTitle("");
  //hSoup->SetYTitle("#frac{N_{obs} - N_{exp}}{#sqrt{N_{obs}}}");
  hSoup->SetYTitle("#frac{N_{obs}}{N_{exp}}");
  hSoup->GetXaxis()->SetLabelSize(0.09);
  hSoup->GetYaxis()->SetLabelSize(0.09);
  hSoup->GetYaxis()->SetTitleSize(0.09);
  hSoup->GetYaxis()->SetTitleOffset(0.5);
  hSoup->Divide(hMCSum);
  hSoup->SetLineWidth(3);
  hSoup->SetMinimum(0.55);
  hSoup->SetMaximum(1.55);
  hSoup->SetStats(kFALSE);
  hSoup->SetFillStyle(0);
  hSoup->Draw("E1");
  TLine *aLine = new TLine(hSoup->GetXaxis()->GetXmin(),1.0,highEnd,1.0);
  aLine->SetLineColor(1);
  aLine->SetLineWidth(2);
  aLine->Draw();

  string plotName;
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + ".png";
  else plotName = "fig_png/hTree_"+hName+Form("_%s",hNameSuffix.c_str())+".png";
  c1->Print(plotName.c_str());

  if(hName.find_last_of("/")<string::npos) plotName = "fig_C/" + hName.substr(hName.find_last_of("/")) + ".C";
  else plotName = "fig_C/hTree_"+hName+Form("_%s",hNameSuffix.c_str())+".C";
  c1->Print(plotName.c_str());

  pad1->SetLogy(1);
  if(hName.find_last_of("/")<string::npos) plotName = "fig_png/" + hName.substr(hName.find_last_of("/")) + "_LogY.png";
  else plotName = "fig_png/hTree_"+hName+Form("_%s",hNameSuffix.c_str())+"_LogY.png";
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
std::pair<float,float> HTTHistograms::getQCDOStoSS(unsigned int iCategory,
               std::pair<float,float> wselOSCorrection,
						   std::pair<float,float> wselSSCorrection,
               unsigned int iSystEffect){

                 ///Return fixed values according to SMH2016 TWiki:
                 ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMTauTau2016#QCD_background_estimation_in_Lta
std::pair<float,float> result(1,0);
if(iCategory == (unsigned int)HTTAnalyzer::jet0) result = std::pair<float,float>(1.0,0.15);
else if(iCategory == (unsigned int)HTTAnalyzer::boosted) result = std::pair<float,float>(1.15,0.15*1.15);
else if(iCategory == (unsigned int)HTTAnalyzer::vbf) result = std::pair<float,float>(1.2,0.30*1.2);
else result = std::pair<float,float>(1.0,0.15);

if(iSystEffect==(unsigned int)sysEffects::QCDSFUp) result.first+=result.second;
if(iSystEffect==(unsigned int)sysEffects::QCDSFDown) result.first-=result.second;

return result;

  std::string hName = "h1DIso";
  std::string hNameSuffix;

  // OS selection
  hNameSuffix =  "_OSnoMuIso_"+std::to_string(iCategory);
  TH1F *hWJetsOS = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix).c_str());
  TH1F *hDYJetsLowMOS = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix).c_str());
  TH1F *hDYJetsOS = get1D_DYJet_Histogram((hName+"DYJets"+hNameSuffix).c_str());
  TH1F *hTTOS = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str());
  TH1F *hSTOS = get1D_ST_Histogram((hName+"ST"+hNameSuffix).c_str());
  TH1F *hVVOS = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str());
  TH1F *hSoupOS = get1DHistogram((hName+"Data"+hNameSuffix).c_str());
  TH1F *hSoupOSb = get1DHistogram((hName+"Data"+hNameSuffix).c_str());

  if(!hSoupOS){
    std::cout<<"No data histograms for OS category!"<<std::endl;
    return  std::make_pair(1.0,0.0);
  }

  TH1F *hEmpty = (TH1F*)hSoupOS->Clone("hEmpty");
  hEmpty->Reset();

  if(!hTTOS) hTTOS = (TH1F*)hEmpty->Clone((hName+"TTbar"+hNameSuffix).c_str());
  if(!hSTOS) hSTOS = (TH1F*)hEmpty->Clone((hName+"ST"+hNameSuffix).c_str());
  if(!hVVOS) hVVOS = (TH1F*)hEmpty->Clone((hName+"DiBoson"+hNameSuffix).c_str());
  if(!hWJetsOS) hWJetsOS = (TH1F*)hEmpty->Clone((hName+"WJets"+hNameSuffix).c_str());
  if(!hDYJetsOS) hDYJetsOS = (TH1F*)hEmpty->Clone((hName+"DYJets"+hNameSuffix).c_str());
  if(!hDYJetsLowMOS) hDYJetsLowMOS = (TH1F*)hEmpty->Clone(("hDYLowM"+hNameSuffix).c_str());
  // SS selection
  hNameSuffix =  "_SSnoMuIso_"+std::to_string(iCategory);
  TH1F *hWJetsSS = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix).c_str());
  TH1F *hDYJetsLowMSS = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix).c_str());
  TH1F *hDYJetsSS = get1D_DYJet_Histogram((hName+"DYJets"+hNameSuffix).c_str());
  TH1F *hTTSS = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str());
  TH1F *hSTSS = get1D_ST_Histogram((hName+"ST"+hNameSuffix).c_str());
  TH1F *hVVSS = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str());
  TH1F *hSoupSS = get1DHistogram((hName+"Data"+hNameSuffix).c_str());
  TH1F *hSoupSSb = get1DHistogram((hName+"Data"+hNameSuffix).c_str());

  if(!hSoupSS){
    std::cout<<"No data histograms for SS category!"<<std::endl;
    return  std::make_pair(1.0,0.0);
  }

  if(!hTTSS) hTTSS = (TH1F*)hEmpty->Clone((hName+"TTbar"+hNameSuffix).c_str());
  if(!hSTSS) hSTSS = (TH1F*)hEmpty->Clone((hName+"ST"+hNameSuffix).c_str());
  if(!hVVSS) hVVSS = (TH1F*)hEmpty->Clone((hName+"DiBoson"+hNameSuffix).c_str());
  if(!hWJetsSS) hWJetsSS = (TH1F*)hEmpty->Clone((hName+"WJets"+hNameSuffix).c_str());
  if(!hDYJetsSS) hDYJetsSS = (TH1F*)hEmpty->Clone((hName+"DYJets"+hNameSuffix).c_str());
  if(!hDYJetsLowMSS) hDYJetsLowMSS = (TH1F*)hEmpty->Clone(("hDYLowM"+hNameSuffix).c_str());

  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYLowM";
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
  scale = lumi;
  hTTOS->Scale(scale);
  hTTSS->Scale(scale);

  sampleName = "ST";
  scale = lumi;
  hSTOS->Scale(scale);
  hSTSS->Scale(scale);

  sampleName = "DiBoson";
  scale = lumi;
  hVVOS->Scale(scale);
  hVVSS->Scale(scale);

  ///Subtract backgrounds other than QCD using MC
  hSoupSS->Add(hWJetsSS,-1);
  hSoupSS->Add(hDYJetsLowMSS,-1);
  hSoupSS->Add(hDYJetsSS,-1);
  hSoupSS->Add(hTTSS,-1);
  hSoupSS->Add(hSTSS,-1);
  hSoupSS->Add(hVVSS,-1);

  hSoupOS->Add(hWJetsOS,-1);
  hSoupOS->Add(hDYJetsLowMOS,-1);
  hSoupOS->Add(hDYJetsOS,-1);
  hSoupOS->Add(hTTOS,-1);
  hSoupOS->Add(hSTOS,-1);
  hSoupOS->Add(hVVOS,-1);

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

  std::string categoryName = HTTAnalyzer::categoryName(iCategory);
  c->Print(TString::Format("fig_png/%s_%s.png",hName.c_str(), categoryName.c_str()).Data());
  c->Print(TString::Format("fig_C/%s_%s.C",hName.c_str(), categoryName.c_str()).Data());

  float param, dparam;
  param=line->GetParameter(0);
  dparam=line->GetParError(0);

  std::cout<<"QCD OS/SS ratio: "<<param<<" +- "<<dparam<<std::endl;

  return std::make_pair(param, dparam);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::getQCDbackground(unsigned int iCategory,
				      std::string varName,
				      std::pair<float,float> wselOSCorrection,
				      std::pair<float,float> wselSSCorrection,
               unsigned int iSystEffect){

  float qcdScale = getQCDOStoSS(iCategory, wselOSCorrection, wselSSCorrection, iSystEffect).first;

  std::string hName = "h1D" + varName;
  std::string hNameSuffix =  "_SS_"+std::to_string(iCategory);
  if(HTTAnalyzer::categoryName(iCategory).find("antiiso")!=std::string::npos) hNameSuffix.replace(hNameSuffix.find("_SS_"),4,"_SSnoMuIso_");
  // SS selection
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix).c_str());
  TH1F *hDYJetsLowM = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix).c_str());
  TH1F *hDYJets = get1D_DYJet_Histogram((hName+"DYJets"+hNameSuffix).c_str());
  TH1F *hTTbar = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str());
  TH1F *hST = get1D_ST_Histogram((hName+"ST"+hNameSuffix).c_str());
  TH1F *hVV = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str());
  TH1F *hggH125 = get1DHistogram((hName+"ggH125"+hNameSuffix).c_str());
  TH1F *hqqH125 = get1DHistogram((hName+"qqH125"+hNameSuffix).c_str());
  TH1F *hSoup = get1DHistogram((hName+"Data"+hNameSuffix).c_str());

  ///Protection against null pointers
  ///Null pointers happen when sample was not read, or there were no
  ///events passing particular selection.
  if(!hSoup) return 0;

  TH1F *hEmpty = (TH1F*)hSoup->Clone("hEmpty");
  hEmpty->Reset();

  if(!hWJets) hWJets = (TH1F*)hEmpty->Clone((hName+"WJets"+hNameSuffix).c_str());
  if(!hDYJetsLowM) hDYJetsLowM = (TH1F*)hEmpty->Clone((hName+"hDYLowM"+hNameSuffix).c_str());
  if(!hDYJets) hDYJets = (TH1F*)hEmpty->Clone((hName+"hDYJets"+hNameSuffix).c_str());
  if(!hTTbar) hTTbar = (TH1F*)hEmpty->Clone((hName+"hTTbar"+hNameSuffix).c_str());
  if(!hST) hST = (TH1F*)hEmpty->Clone((hName+"hST"+hNameSuffix).c_str());
  if(!hVV) hVV = (TH1F*)hEmpty->Clone((hName+"hDiBoson"+hNameSuffix).c_str());
  if(!hqqH125) hqqH125 = (TH1F*)hEmpty->Clone((hName+"qqH125"+hNameSuffix).c_str());
  if(!hggH125) hggH125 = (TH1F*)hEmpty->Clone((hName+"ggH125"+hNameSuffix).c_str());
  //////////////////////////////////////////////////////////////////////
  float lumi = getLumi();
  ///Normalise MC histograms according to cross sections
  std::string sampleName = "DYLowM";
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
  scale = lumi;
  hTTbar->Scale(scale);

  sampleName = "qqH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hqqH125->Scale(scale);

  sampleName = "ggH125";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hggH125->Scale(scale);

  sampleName = "ST";
  scale = lumi;
  hST->Scale(scale);

  sampleName = "DiBoson";
  scale = lumi;
  hVV->Scale(scale);

  hSoup->SetName(("h1D"+varName+"QCDEstimate_"+std::to_string(iCategory)).c_str());
  hSoup->Add(hWJets,-1);
  hSoup->Add(hDYJetsLowM,-1);
  hSoup->Add(hDYJets,-1);
  hSoup->Add(hTTbar,-1);
  hSoup->Add(hST,-1);
  hSoup->Add(hVV,-1);
  //hSoup->Add(hqqH125,-1);
  //hSoup->Add(hggH125,-1);

  ///Clean up the QCD shape, and remove fluctuations around 0 counts.
  for(unsigned int iBinX=0;iBinX<=hSoup->GetNbinsX();++iBinX){
    if(hSoup->GetBinContent(iBinX)<3.0) hSoup->SetBinContent(iBinX,0);
  }

  hSoup->Scale(qcdScale);

  return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getWNormalisation(unsigned int iCategory,
      std::string selName,
    unsigned int iSystEffect){

  if(iCategory==(unsigned int)(HTTAnalyzer::jet0)) iCategory = HTTAnalyzer::wjets_jet0;
  else if(iCategory==(unsigned int)(HTTAnalyzer::boosted)) iCategory = HTTAnalyzer::wjets_boosted;
  else if(iCategory==(unsigned int)(HTTAnalyzer::vbf)) iCategory = HTTAnalyzer::wjets_vbf;
  else iCategory = HTTAnalyzer::W;

  std::string hNameSuffix =  "_"+selName+"_"+std::to_string(iCategory);

  std::string hName = "h1DMassTrans";
  TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix).c_str());
  TH1F *hDYJets = get1D_DYJet_Histogram((hName+"DYJets"+hNameSuffix).c_str());
  TH1F *hDYJetsLowM = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix).c_str());
  TH1F *hTT = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix).c_str());
  TH1F *hST = get1D_ST_Histogram((hName+"ST"+hNameSuffix).c_str());
  TH1F *hVV = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix).c_str());
  TH1F *hQCD = (TH1F*)getQCDbackground(iCategory,"MassTrans",wselOSCorrection, wselSSCorrection, iSystEffect);
  TH1F *hSoup = get1DHistogram((hName+"Data"+hNameSuffix).c_str());
  float lumi = getLumi();

  TH1F *hEmpty = (TH1F*)hWJets->Clone("hEmpty");
  hEmpty->Reset();

  if(!hDYJets) hDYJets = (TH1F*)hEmpty->Clone((hName+"hDYJets"+hNameSuffix).c_str());
  if(!hDYJetsLowM) hDYJetsLowM = (TH1F*)hEmpty->Clone((hName+"hDYLowM"+hNameSuffix).c_str());
  if(!hTT) hTT = (TH1F*)hEmpty->Clone((hName+"hTTbar"+hNameSuffix).c_str());
  if(!hST) hST = (TH1F*)hEmpty->Clone((hName+"hST"+hNameSuffix).c_str());
  if(!hVV) hVV = (TH1F*)hEmpty->Clone((hName+"hDiBoson"+hNameSuffix).c_str());
  if(!hQCD) hQCD = (TH1F*)hEmpty->Clone((hName+"hQCD"+hNameSuffix).c_str());

  ///Normalise MC histograms according to cross sections
  std::string sampleName = "WJets";
  float weight = getSampleNormalisation(sampleName);
  float scale = weight*lumi;
  hWJets->Scale(scale);

  sampleName = "DYLowM";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJetsLowM->Scale(scale);

  sampleName = "DYJets";
  weight = getSampleNormalisation(sampleName);
  scale = weight*lumi;
  hDYJets->Scale(scale);

  sampleName = "TTbar";
  scale = lumi;
  hTT->Scale(scale);

  sampleName = "ST";
  scale = lumi;
  hST->Scale(scale);

  sampleName = "DiBoson";
  scale = lumi;
  hVV->Scale(scale);

  // Create a histogram with data minus backgrounds: DYJets, hTT, Other
  TH1F* datamtlo = (TH1F*)hSoup->Clone("datamtlo");
  datamtlo->Add(hDYJets,-1);
  datamtlo->Add(hDYJetsLowM,-1);
  datamtlo->Add(hTT,-1);
  datamtlo->Add(hST,-1);
  datamtlo->Add(hVV,-1);
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
  float inthST = hST->Integral(0,hST->GetNbinsX()+1);
  float inthVV = hVV->Integral(0,hVV->GetNbinsX()+1);
  float inthQCD = hQCD ? hQCD->Integral(0,hQCD->GetNbinsX()+1) : 0;
  dweight=((inthSoup+inthDYJets+inthTT+inthST+inthVV)/inthWJets/inthWJets+intdata*intdata/(inthWJets*inthWJets*inthWJets));
  dweight=sqrt(dweight);
  cout<<"Selecion name: "<<hNameSuffix<<std::endl;
  cout<<" DATA:              "<<inthSoup<<endl
      <<" DATA - MC(!WJets): "<<intdata<<endl
      <<" MC WJets           "<<inthWJets<<endl
      <<" DYJets:            "<<inthDYJets<<endl
      <<" DYLowM:            "<<inthDYJetsLowM<<endl
      <<" TTbar:             "<<inthTT<<endl
      <<" single T:          "<<inthST<<endl
      <<" DiBoson:           "<<inthVV<<endl
      <<" QCD:               "<<inthQCD<<endl;
  cout<<"WJets scale:"<<weight<<" dweight "<<dweight<<endl;

float WSFUncertainty = 0.0;
  if(iCategory==(unsigned int)HTTAnalyzer::jet0) WSFUncertainty = 0.1;
  if(iCategory==(unsigned int)HTTAnalyzer::boosted) WSFUncertainty = 0.1;
  if(iCategory==(unsigned int)HTTAnalyzer::vbf) WSFUncertainty = 0.3;

  if(iSystEffect==(unsigned int)sysEffects::WSFUp){
    weight*=(1+WSFUncertainty);
    weight*=(1+WSFUncertainty);
  }

  if(iSystEffect==(unsigned int)sysEffects::WSFDown){
    weight*=(1-WSFUncertainty);
    weight*=(1-WSFUncertainty);
  }
std::cout << "WSFUncertainty: " << WSFUncertainty <<std::endl;
cout<<"WJets scale with uncertainty: "<<weight<<" dweight "<<dweight<<endl;

  return std::make_pair(weight, dweight);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
