#include <iostream>
#include <sstream>
#include <cmath>

#include "commonUtils.h"
#include "HTTHistograms.h"
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
float HTTHistograms::getSampleNormalisation(std::string sampleName){

        std::string hName = "h1DStats"+sampleName;
        TH1F *hStats = new TH1F("h1DStats","",11,-0.5,10.5);
        hStats->SetBinContent(1,1);
        hStats->SetBinContent(2,1);
        hStats->SetBinContent(3,1);

        if(sampleName.find("DYJetsMatch")!=std::string::npos || sampleName=="DYJets") { /*DYJets50 are normalised for analysed events and preselection for each nJets sample*/}
        else if(sampleName.find("WJets")!=std::string::npos) { /*WJets are normalised for analysed events and preselection for each nJets sample*/ }
        else if(sampleName=="ST") { /*WJets are normalised for analysed events and preselection for each nJets sample*/ }
        else hStats = get1DHistogram(hName);

        if(!hStats) hStats = get1D_TauMatchJetSum(hName,true,false);

        if(!hStats) return 0;

        float genPresEff = 1.0;
        float recoPresEff = hStats->GetBinContent(3)/hStats->GetBinContent(2);
        float presEff = genPresEff*recoPresEff;
        float kFactor = 1.0;

        float crossSection = HTTAnalysis::getCrossSection(sampleName);
        int nEventsAnalysed = hStats->GetBinContent(1);

        //test HACK fixme!!!
        std::vector<std::string> sampleNames = {"TTbar", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};
        for(auto sampleNameTmp : sampleNames) {
                if(sampleName.find(sampleNameTmp+"Match")!=std::string::npos) {hStats = get1D_TauMatchJetSum(hName, true, false); nEventsAnalysed=hStats->GetBinContent(1); }
        }
        //test

        float weight = crossSection*presEff/nEventsAnalysed;
        if(presEff<0 || fabs(fabs(crossSection)-1.0)<1e-5) weight = 1.0;
        /*
           outputStream<<"Sample name: "<<sampleName<<" ";
           outputStream<<"Xsection: "<<crossSection<<" [pb] "<<" ";
           outputStream<<"Events analyzed: "<<nEventsAnalysed<<" ";
           outputStream<<"Reco preselection efficiency: "<<recoPresEff<<" ";
           outputStream<<"Final weight: "<<weight<<std::endl;
         */
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
HTTHistograms::~HTTHistograms(){
}
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

        if(sumDecayModes) {
                ///Strip decay name if any.
                for(auto decayName : decayNames) {hName.ReplaceAll(decayName.c_str(), ""); }
                for(auto decayName : decayNames) {
                        TString hNameTmp = hName;
                        if(hNameTmp.First("_")>0) hNameTmp.Replace(hNameTmp.First("_"),1,(decayName+"_").c_str()); else hNameTmp.Append(decayName.c_str());
                        TH1F *hDecayMode = 0;
                        if(sumJetBins) hDecayMode = get1D_VJetSum(hNameTmp.Data());
                        else hDecayMode = get1DHistogram(hNameTmp.Data());
                        if(!hSum && hDecayMode) {
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
TH1F *HTTHistograms::get1D_WJet_Histogram(const std::string& name){
        return get1D_VJetSum(name);
}
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

        if(hName.find("W")!=std::string::npos &&
	   hName.find("DY")==std::string::npos) {
                sampleName = hName.substr(hName.find("W"));
                std::string selName = sampleName.substr(sampleName.find("_"));
                sampleName = sampleName.substr(0,sampleName.size()-selName.size());
                hNJetsStats = get1DHistogram("h1DStats"+sampleName);
        }
        if(hName.find("DY")!=std::string::npos) {
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
           sampleName.find("AllJets")!=std::string::npos) {
                TString allJetsName = "h1DStats"+sampleName;
                if(sampleName.find("0Jets")!=std::string::npos) allJetsName.ReplaceAll("0Jets","AllJets");
                if(sampleName.find("AllJets")!=std::string::npos) allJetsName.ReplaceAll("AllJets","0Jets");
                TH1F *hAllJetsStats = 0;
                bool sumDecayModes = true;
                bool sumJetBins = false;
                if(hName.find("W")!=std::string::npos && hName.find("DY")==std::string::npos) hAllJetsStats = get1DHistogram(allJetsName.Data());
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
        if(name.find("W")!=std::string::npos && name.find("DY")==std::string::npos) jetsLOSigma = {50380, 9644.5, 3144.5, 954.8, 485.6};
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
TH1F* HTTHistograms::get1D_SumPattern_Histogram(const std::string& name, std::string pattern, std::vector<std::string> sampleNames, std::string tauMatchSuffix){

        TString hName = name;
        TH1F *hSum = 0;

        for(auto sampleName : sampleNames) {
                TString hNameTmp = hName;
                hNameTmp.ReplaceAll(pattern,(sampleName+tauMatchSuffix).c_str());
                TH1F *histo = get1DHistogram(hNameTmp.Data());

                if(tauMatchSuffix == "" && !histo) histo = get1D_TauMatchJetSum(hNameTmp.Data(), true, false);

                if(!hSum && histo) {
                        hSum = (TH1F*)histo->Clone(name.c_str());
                        hSum->Reset();
                }
                if(histo && hSum) {
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
TH1F* HTTHistograms::get1D_VV_Histogram(const std::string& name, std::string tauMatchSuffix){

        std::vector<std::string> sampleNamesVV = {"ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};
        return get1D_SumPattern_Histogram(name, "DiBoson", sampleNamesVV, tauMatchSuffix);

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

        std::string templateName = "TemplateNotFound: "+name;
        if(name.find("hProf")!=std::string::npos && name.find("VsMag")!=std::string::npos) templateName = "hProfVsMagTemplate";
        if(name.find("hProf")!=std::string::npos && name.find("VsPt")!=std::string::npos) templateName = "hProfVsPtTemplate";
        if(name.find("hProf")!=std::string::npos && name.find("VsCos")!=std::string::npos) templateName = "hProfVsCosTemplate";

        if(name.find("h1DNPV")!=std::string::npos) templateName = "h1DNPVTemplate";
        if(name.find("h1DNPU")!=std::string::npos) templateName = "h1DNPUTemplate";
        if(name.find("h1DMass")!=std::string::npos) templateName = "h1DMassTemplate";
        if(name.find("h1DBigMass")!=std::string::npos) templateName = "h1DBigMassTemplate";
        if(name.find("h1DStats")!=std::string::npos) templateName = "h1DStatsTemplate";
        if(name.find("h1DPt")!=std::string::npos) templateName = "h1DPtTemplate";
        if(name.find("h1DEta")!=std::string::npos && name.find("Jet")==std::string::npos) templateName = "h1DEtaTemplate";
        if(name.find("h1DEta")!=std::string::npos && name.find("Jet")!=std::string::npos) templateName = "h1DJetEtaTemplate";
        if(name.find("h1DDeltaEta")!=std::string::npos) templateName = "h1DDeltaEtaTemplate";
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

        if(!histosInitialized_) {
                add1DHistogram("h1DStatsTemplate","",21,-0.5,20.5,file_);
                add1DHistogram("h1DNPVTemplate",";Number of PV; Events",61,-0.5,60.5,file_);
                add1DHistogram("h1DNPUTemplate",";Number of PV; Events",600,0,60,file_);
                add1DHistogram("h1DMassTemplate",";mass [GeV/c^{2}]; Events",35,0,350,file_);
                add1DHistogram("h1DBigMassTemplate",";mass [GeV/c^{2}]; Events",25,0,1500,file_);
                add1DHistogram("h1DPtTemplate",";p_{T}; Events",20,0,100,file_);
                add1DHistogram("h1DEtaTemplate",";#eta; Events",24,-2.4,2.4,file_);
                add1DHistogram("h1DJetEtaTemplate",";#eta; Events",50,-5,5,file_);
                add1DHistogram("h1DDeltaEtaTemplate",";#Delta#eta; Events",50,0,10,file_);
                add1DHistogram("h1DPhiTemplate",";#phi; Events",2*9,-M_PI,2*M_PI,file_);
                add1DHistogram("h1DCosPhiTemplate",";cos(#phi); Events",10,-1.0,1.0,file_);
                add1DHistogram("h1DCSVBtagTemplate",";CSV btag; Events",20,0,1,file_);
                add1DHistogram("h1DIsoTemplate",";Isolation; Events",20,0,0.3,file_);
                add1DHistogram("h1DIDTemplate",";ID; Events",20,0.8,1,file_);
                add1DHistogram("h1DVxPullTemplate",";#phi^{*} [rad]; Events",11,-0.01,0.01,file_);
                add1DHistogram("h1DyTauTemplate",";yTau; Events",15,-1,1,file_);
                add1DHistogram("h1DnPCATemplate",";#hat{n}_{RECO}>; Events",20,0,0.02,file_);

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
                for(unsigned int iBin=0; iBin<=12; ++iBin) phiBins.push_back(iBin*2.0*M_PI/12);
                addRollHistogram("h1DUnRollMassSVPhiCPTemplate","#phi_{IP,IP} CP vs SV Mass; Events;#phi_{IP,IP} CP",phiBins, svMassBins, file_);
                addRollHistogram("h1DUnRollMassSVYCPTemplate","#phi_{IP,#rho} CP vs SV Mass; Events;#phi_{IP,#rho} CP", phiBins, svMassBins, file_);

                addProfile("hProfVsMagTemplate","",10,0,0.015,file_);
                addProfile("hProfVsPtTemplate","",20,15,55,file_);
                addProfile("hProfVsCosTemplate","",20,-1,1,file_);
        }
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

        std::cout<<"HTTHistograms::finalizeHistograms() START"<<std::endl;

        AnalysisHistograms::finalizeHistograms();
        myCategoryRejester  = aCategoryRejester;
        unsigned int myNumberOfCategories = myCategoryRejester.size();

        ///
        std::vector<std::string> mainCategoryNames = {"0jet","boosted", "vbf",
                                                      "antiIso_0jet", "antiIso_boosted", "antiIso_vbf",
                                                      "antiIso_0jet_QCD","0jet_W_QCD",
                                                      "0jet_QCD", "boosted_QCD", "vbf_QCD",
                                                      "0jet_W", "boosted_W", "vbf_W"
						      ,"inclusive","btag","nobtag"
	};

        std::vector <unsigned int> mainCategoriesRejester;
        for(unsigned int iCategory=0; iCategory<myCategoryRejester.size(); ++iCategory) {
                for(auto nameIt : mainCategoryNames) {
                        if(myCategoryRejester[iCategory]->name()==nameIt) mainCategoriesRejester.push_back(iCategory);
                }
        }

        gErrorIgnoreLevel = kBreak;
        //////////////
        ///Control regions plots
        for(auto iCategory: mainCategoriesRejester) {

                //plotCPhistograms(iCategory);

                plotStack(iCategory, "MassSV");
                plotStack(iCategory, "MassVis");
                plotStack(iCategory, "MassTrans");
                plotStack(iCategory, "UnRollTauPtMassVis");
                plotStack(iCategory, "UnRollHiggsPtMassSV");
                plotStack(iCategory, "UnRollMjjMassSV");
                plotStack(iCategory, "UnRollMassSVPhiCP");
                plotStack(iCategory, "UnRollMassSVYCP");

                plotStack(iCategory, "PtLeg1");
                plotStack(iCategory, "EtaLeg1");
                plotStack(iCategory, "IsoLeg1");

                plotStack(iCategory, "PtLeg2");
                plotStack(iCategory, "EtaLeg2");
                plotStack(iCategory, "IDLeg2");
                plotStack(iCategory, "StatsLeg2DecayMode");
                plotStack(iCategory, "PtLeg2LeadingTk");

                plotStack(iCategory, "PhiLeg1");
                plotStack(iCategory, "PhiLeg2");

                plotStack(iCategory, "PtMET");
                plotStack(iCategory, "PtLeg1Leg2MET");

                plotStack(iCategory, "StatsNJ30");
                plotStack(iCategory, "PtLeadingJet");
                plotStack(iCategory, "EtaLeadingJet");
                plotStack(iCategory, "StatsNBTag");
                plotStack(iCategory, "CSVBtagLeadingJet");
                plotStack(iCategory, "PtLeadingBJet");
                plotStack(iCategory, "EtaLeadingBJet");
                plotStack(iCategory, "BigMass2Jet");
                plotStack(iCategory, "DeltaEta2Jet");

                plotStack(iCategory, "nPCALeg1");
                plotStack(iCategory, "nPCALeg2");
                plotStack(iCategory, "Phi-nVectors");
                plotStack(iCategory, "Phi-nVecIP");
                plotStack(iCategory, "NPV");
        }

        ///Make systematic effect histos.
        for(unsigned int iSystEffect = (unsigned int)HTTAnalysis::NOMINAL;
            iSystEffect<=(unsigned int)HTTAnalysis::ZmumuDown; ++iSystEffect) {
              if(iSystEffect==(unsigned int)HTTAnalysis::DUMMY_SYS) continue;
                for(auto iCategory: mainCategoriesRejester) {
                        plotStack(iCategory, "MassSV", iSystEffect);
                        plotStack(iCategory, "MassTrans", iSystEffect);
                        plotStack(iCategory, "MassVis", iSystEffect);
                        plotStack(iCategory, "UnRollTauPtMassVis", iSystEffect);
                        plotStack(iCategory, "UnRollHiggsPtMassSV", iSystEffect);
                        plotStack(iCategory, "UnRollMjjMassSV", iSystEffect);
                        plotStack(iCategory, "UnRollMassSVPhiCP", iSystEffect);
                        plotStack(iCategory, "UnRollMassSVYCP", iSystEffect);
                }
        }
        ofstream eventCountFile("eventCount.txt",ios::out | ios::app);
        outputStream<<"HTTHistograms compilation time: "<<__TIMESTAMP__<<std::endl;
        eventCountFile<<outputStream.str();
        eventCountFile.close();

        std::cout<<"HTTHistograms::finalizeHistograms() END"<<std::endl;

}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotCPhistograms(unsigned int iCategory){

        std::string hNameSuffix = "";
        plot_HAZ_Histograms("Phi-nVecIP-yTauNeg",hNameSuffix+"_GenNoOfflineSel");
        plot_HAZ_Histograms("Phi-nVecIP-yTauPos",hNameSuffix+"_GenNoOfflineSel");
        plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_GenNoOfflineSel");
        plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_GenNoOfflineSel");
	std::string categoryName = myCategoryRejester[iCategory]->name();
        hNameSuffix =  "_"+categoryName;

        //plot_HAZ_Histograms("Phi-nVecIP-yTauNeg",hNameSuffix+"_Gen");
        //plot_HAZ_Histograms("Phi-nVecIP-yTauPos",hNameSuffix+"_Gen");
        plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_Gen");
        plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_Gen");

        plot_HAZ_Histograms("Phi-nVectors",hNameSuffix+"_RefitPV");
        plot_HAZ_Histograms("Phi-nVecIP",hNameSuffix+"_RefitPV");

        plotnPCA("ggHTT125"+hNameSuffix);
        plotnPCA("ATT"+hNameSuffix);
        plotnPCA("DYJets"+hNameSuffix);
        plotnPCA("WJets"+hNameSuffix);

        plotPhiDecayPlanes("Phi-nVectorsggHTT125"+hNameSuffix);
        plotPhiDecayPlanes("Phi-nVectorsATT"+hNameSuffix);
        plotPhiDecayPlanes("Phi-nVectorsDYJets"+hNameSuffix);

        plotPhiDecayPlanes("Phi-nVecIPggHTT125"+hNameSuffix);
        plotPhiDecayPlanes("Phi-nVecIPATT"+hNameSuffix);
        plotPhiDecayPlanes("Phi-nVecIPDYJets"+hNameSuffix);

        plotProfiles("hProfRecoVsMagGen","ggHTT125"+hNameSuffix);
        plotProfiles("hProfPhiVsMag","ggHTT125"+hNameSuffix);

        plotVerticesPulls("h1DVxPullXggHTT125"+hNameSuffix);
        plotVerticesPulls("h1DVxPullYggHTT125"+hNameSuffix);
        plotVerticesPulls("h1DVxPullZggHTT125"+hNameSuffix);

        plotPhiDecayPlanes("Phi-nVecIP-yTauPosggHTT125"+hNameSuffix);
        plotPhiDecayPlanes("Phi-nVecIP-yTauNegggHTT125"+hNameSuffix);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void HTTHistograms::plotnPCA(const std::string & type){

        TH1F* h1DLeg2 = 0;
        TH1F* h1DLeg1 = 0;

        if(type.find("DYJets")!=std::string::npos) {
                std::string typeztt="DYJetsMatchT"+type.substr(type.find("DYJets")+6);
                h1DLeg2 = get1D_TauMatchJetSum(("h1DnPCALeg2"+typeztt).c_str(),false,true);
                h1DLeg1 = get1D_TauMatchJetSum(("h1DnPCALeg1"+typeztt).c_str(),false,true);
        }
        else if(type.find("WJets")!=std::string::npos) {
                h1DLeg2 = get1D_WJet_Histogram("h1DnPCALeg2"+type);
                h1DLeg1 = get1D_WJet_Histogram("h1DnPCALeg1"+type);
        }
        else{
                h1DLeg2 = get1DHistogram("h1DnPCALeg2"+type);
                h1DLeg1 = get1DHistogram("h1DnPCALeg1"+type);
        }

        if(!h1DLeg2 || !h1DLeg1) return;

        TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);

        TLegend l(0.15,0.12,0.35,0.22,NULL,"brNDC");
        l.SetTextSize(0.05);
        l.SetFillStyle(4000);
        l.SetBorderSize(0);
        l.SetFillColor(10);

        h1DLeg2->SetLineWidth(3);
        h1DLeg2->Scale(1.0/h1DLeg2->Integral(0,h1DLeg2->GetNbinsX()+1));

        h1DLeg1->SetLineWidth(3);
        h1DLeg1->SetLineColor(2);
        h1DLeg1->Scale(1.0/h1DLeg1->Integral(0,h1DLeg1->GetNbinsX()+1));
        h1DLeg1->GetYaxis()->SetTitleOffset(1.5);
        h1DLeg1->SetStats(kFALSE);
        h1DLeg1->SetYTitle("Events");
        h1DLeg1->SetXTitle("|n_{RECO}|");

        h1DLeg1->Draw();
        h1DLeg2->Draw("same");

        l.AddEntry(h1DLeg2,"leg2 (hadronic tau)");
        l.AddEntry(h1DLeg1,"leg1(leptonic/leading had. tau)");
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

        if(hName.find("2D")!=std::string::npos) {
                TProfile* hProfile_AOD = this->get2DHistogram((hName+"_AODPV"))->ProfileX();
                TProfile* hProfile_Refit = this->get2DHistogram((hName+"_RefitPV"))->ProfileX();

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

        TH1F* h1D_AOD = this->get1DHistogram((hName+"_AODPV"));
        TH1F* h1D_Refit = this->get1DHistogram((hName+"_RefitPV"));

        if(h1D_AOD && h1D_Refit) {

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

        if(h1DGen && h1DRefit && h1DAOD) {
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

                if(hName.find("RecoVsMagGen")!=std::string::npos) {
                        h1DGen->SetYTitle("<|n_{RECO}|>");
                        h1DGen->SetMinimum(0);
                }
                if(hName.find("MagVsPt")!=std::string::npos) {
                        h1DGen->SetYTitle("<|n_{GEN}|>");
                        h1DGen->SetXTitle("p_{T}^{leading tk.}");
                        h1DGen->SetMinimum(0);
                }
                if(hName.find("PtVsMag")!=std::string::npos) {
                        h1DGen->SetYTitle("p_{T}^{leading tk.}");
                        h1DGen->SetXTitle("<|n_{GEN}|>");
                        h1DGen->SetMinimum(0);
                }
                if(hName.find("MagVsCos")!=std::string::npos) {
                        h1DGen->SetYTitle("<|n_{GEN}|>");
                        h1DGen->SetXTitle("cos(#phi)");
                }

                h1DGen->Draw();
                h1DAOD->Draw("same");
                h1DRefit->Draw("same");

                if(hName.find("RecoVsMagGen")!=std::string::npos) {
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

        if(h1DGen) {
                h1DGen->SetLineWidth(4);
                h1DGen->Scale(1.0/h1DGen->Integral(0,h1DGen->GetNbinsX()+1));
                h1DGen->SetLineColor(1);
        }

        if(h1DAODPV) {
                h1DAODPV->SetLineWidth(3);
                h1DAODPV->Scale(1.0/h1DAODPV->Integral(0,h1DAODPV->GetNbinsX()+1));
                h1DAODPV->SetLineColor(2);
        }

        if(h1DGenPV) {
                h1DGenPV->SetLineWidth(3);
                h1DGenPV->Scale(1.0/h1DGenPV->Integral(0,h1DGenPV->GetNbinsX()+1));
                h1DGenPV->SetLineColor(3);
        }

        if(h1DRefitPV) {
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

                l.AddEntry(h1DRefitPV,"nPCA with refit. PV");
                if(h1DGenPV) {
                        h1DGenPV->Draw("HISTO same");
                        l.AddEntry(h1DGenPV,"nPCA with gen. PV");
                }
                if(h1DAODPV) {
                        h1DAODPV->Draw("HISTO same");
                        l.AddEntry(h1DAODPV,"nPCA with AOD PV");
                }
                if(h1DGen) {
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

        TString name = "h1D"+hName+"ggHTT125"+sysType;
        TH1F* h_h = this->get1DHistogram(name.Data());
        name = "h1D"+hName+"ATT"+sysType;
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
void HTTHistograms::plotSingleHistogram(std::string hName){

        TH1F* h1D = get1DHistogram(hName);
        if(!h1D) return;

        TCanvas* c = new TCanvas("AnyHistogram","AnyHistogram",
                                 460,500);

        TLegend l(0.15,0.78,0.35,0.87,NULL,"brNDC");
        l.SetTextSize(0.05);
        l.SetFillStyle(4000);
        l.SetBorderSize(0);
        l.SetFillColor(10);

        if(h1D) {
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
THStack*  HTTHistograms::plotStack(unsigned int iCategory,
                                   std::string varName, unsigned int iSystEffect){

        std::string hName = "h1D"+varName;
        std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
        std::string hNameSuffix = "_"+categoryName+systEffectName;

        TH1F *hggHiggs120 = get1DHistogram((hName+"ggHTT120"+hNameSuffix));
        TH1F *hggHiggs125 = get1DHistogram((hName+"ggHTT125"+hNameSuffix));
        TH1F *hggHiggs130 = get1DHistogram((hName+"ggHTT130"+hNameSuffix));

        TH1F *hqqHiggs120 = get1DHistogram((hName+"qqHTT120"+hNameSuffix));
        TH1F *hqqHiggs125 = get1DHistogram((hName+"qqHTT125"+hNameSuffix));
        TH1F *hqqHiggs130 = get1DHistogram((hName+"qqHTT130"+hNameSuffix));

        TH1F *hZHiggs120 = get1DHistogram((hName+"ZHTT120"+hNameSuffix));
        TH1F *hZHiggs125 = get1DHistogram((hName+"ZHTT125"+hNameSuffix));
        TH1F *hZHiggs130 = get1DHistogram((hName+"ZHTT130"+hNameSuffix));

        TH1F *hWplusHiggs120 = get1DHistogram((hName+"WplusHTT120"+hNameSuffix));
        TH1F *hWplusHiggs125 = get1DHistogram((hName+"WplusHTT125"+hNameSuffix));
        TH1F *hWplusHiggs130 = get1DHistogram((hName+"WplusHTT130"+hNameSuffix));

        TH1F *hWminusHiggs120 = get1DHistogram((hName+"WminusHTT120"+hNameSuffix));
        TH1F *hWminusHiggs125 = get1DHistogram((hName+"WminusHTT125"+hNameSuffix));
        TH1F *hWminusHiggs130 = get1DHistogram((hName+"WminusHTT130"+hNameSuffix));

        TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix));
        TH1F *hTTbarJ = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix),"MatchJ");
        TH1F *hTTbarT = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix),"MatchT");
        TH1F *hST = get1D_ST_Histogram((hName+"ST"+hNameSuffix));
        TH1F *hVVJ = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix), "MatchJ");
        TH1F *hVVT = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix), "MatchT");
        TH1F *hDYJetsLowM = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix));

        bool sumDecayModes = false;
        bool sumJetBins = true;
        TH1F *hDYJetsZJ = get1D_TauMatchJetSum((hName+"DYJetsMatchJ"+hNameSuffix), sumDecayModes, sumJetBins);
        TH1F *hDYJetsZL = get1D_TauMatchJetSum((hName+"DYJetsMatchL"+hNameSuffix), sumDecayModes, sumJetBins);
        TH1F *hDYJetsZTT = get1D_TauMatchJetSum((hName+"DYJetsMatchT"+hNameSuffix), sumDecayModes, sumJetBins);

        TH1F *hEWK2Jets = get1D_EWK2JetsSum(hName+"EWK2Jets"+hNameSuffix);

        TH1F *hSoup = get1DHistogram((hName+"Data"+hNameSuffix),true);
        pair<float,float> qcdOStoSS = getQCDControlToSignal(iCategory, iSystEffect);
        TH1F *hQCD = (TH1F*)getQCDbackground(iCategory, varName, iSystEffect);
        TH1F *hQCD_MC =  get1DHistogram((hName+"QCD_MC"+hNameSuffix));

        ///Protection against null pointers
        ///Null pointers happen when sample was not read, or there were no
        ///events passing particular selection.
        if(!hSoup) {
                //std::cout<<"No data events for "<<hName<<"Data"<<hNameSuffix
                  //       <<" ("<<categoryName<<")"
                    //     <<std::endl;
                return 0;
        }

        TH1F *hEmpty = (TH1F*)hSoup->Clone("hEmpty");
        hEmpty->Reset();

        hSoup->SetDirectory(myDirCopy);

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

        if(!hggHiggs120) hggHiggs120 = (TH1F*)hEmpty->Clone((hName+"ggHTT120"+hNameSuffix).c_str());
        if(!hggHiggs125) hggHiggs125 = (TH1F*)hEmpty->Clone((hName+"ggHTT125"+hNameSuffix).c_str());
        if(!hggHiggs130) hggHiggs130 = (TH1F*)hEmpty->Clone((hName+"ggHTT130"+hNameSuffix).c_str());

        if(!hqqHiggs120) hqqHiggs120 = (TH1F*)hEmpty->Clone((hName+"qqHTT120"+hNameSuffix).c_str());
        if(!hqqHiggs125) hqqHiggs125 = (TH1F*)hEmpty->Clone((hName+"qqHTT125"+hNameSuffix).c_str());
        if(!hqqHiggs130) hqqHiggs130 = (TH1F*)hEmpty->Clone((hName+"qqHTT130"+hNameSuffix).c_str());

        if(!hZHiggs120) hZHiggs120 = (TH1F*)hEmpty->Clone((hName+"ZHTT120"+hNameSuffix).c_str());
        if(!hZHiggs125) hZHiggs125 = (TH1F*)hEmpty->Clone((hName+"ZHTT125"+hNameSuffix).c_str());
        if(!hZHiggs130) hZHiggs130 = (TH1F*)hEmpty->Clone((hName+"ZHTT130"+hNameSuffix).c_str());

        if(!hWplusHiggs120) hWplusHiggs120 = (TH1F*)hEmpty->Clone((hName+"WplusHTT120"+hNameSuffix).c_str());
        if(!hWplusHiggs125) hWplusHiggs125 = (TH1F*)hEmpty->Clone((hName+"WplusHTT125"+hNameSuffix).c_str());
        if(!hWplusHiggs130) hWplusHiggs130 = (TH1F*)hEmpty->Clone((hName+"WplusHTT130"+hNameSuffix).c_str());

        if(!hWminusHiggs120) hWminusHiggs120 = (TH1F*)hEmpty->Clone((hName+"WminusHTT120"+hNameSuffix).c_str());
        if(!hWminusHiggs125) hWminusHiggs125 = (TH1F*)hEmpty->Clone((hName+"WminusHTT125"+hNameSuffix).c_str());
        if(!hWminusHiggs130) hWminusHiggs130 = (TH1F*)hEmpty->Clone((hName+"WminusHTT130"+hNameSuffix).c_str());

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

        float lumi = HTTAnalysis::getLumi();
        float weight = 0;
        float scale = 0;

        std::string sampleName = "WJets";
        pair<float,float> dataToMCScale = getWNormalisation(iCategory, iSystEffect);
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi*dataToMCScale.first;
        hWJets->Scale(scale);

        ///EWK 2Jets scaled for preselection during stiching step
        sampleName = "EWK2Jets";
        scale = lumi;
        hEWK2Jets->Scale(scale);

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

        sampleName = "ggHTT120";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hggHiggs120->Scale(scale);

        sampleName = "qqHTT120";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hqqHiggs120->Scale(scale);

        sampleName = "ggHTT125";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hggHiggs125->Scale(scale);

        sampleName = "qqHTT125";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hqqHiggs125->Scale(scale);

        sampleName = "ggHTT130";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hggHiggs130->Scale(scale);

        sampleName = "qqHTT130";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hqqHiggs130->Scale(scale);

        sampleName = "ZHTT120";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hZHiggs120->Scale(scale);

        sampleName = "ZHTT125";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hZHiggs125->Scale(scale);

        sampleName = "ZHTT130";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hZHiggs130->Scale(scale);

        sampleName = "WminusHTT120";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWminusHiggs120->Scale(scale);

        sampleName = "WminusHTT125";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWminusHiggs125->Scale(scale);

        sampleName = "WminusHTT130";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWminusHiggs130->Scale(scale);

        sampleName = "WplusHTT120";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWplusHiggs120->Scale(scale);

        sampleName = "WplusHTT125";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWplusHiggs125->Scale(scale);

        sampleName = "WplusHTT130";
        weight = getSampleNormalisation(sampleName);
        scale = weight*lumi;
        hWplusHiggs130->Scale(scale);

        //Join Wplus and Wminus processes////////////////////
        TH1F *hWHiggs120 = hWplusHiggs120;
        hWHiggs120->Add(hWminusHiggs120);
        hWHiggs120->SetName((hName+"WHTT120"+hNameSuffix).c_str());
        TH1F *hWHiggs125 = hWplusHiggs125;
        hWHiggs125->Add(hWminusHiggs125);
        hWHiggs125->SetName((hName+"WHTT125"+hNameSuffix).c_str());
        TH1F *hWHiggs130 = hWplusHiggs130;
        hWHiggs130->Add(hWminusHiggs130);
        hWHiggs130->SetName((hName+"WHTT130"+hNameSuffix).c_str());
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
        hDYJetsLowM->SetFillColor(kGreen+4);
        hQCD->SetFillColor(kMagenta-10);
        hHiggs->SetFillColor(kCyan+2);

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

        outputStream<<"Event count summary for selecion name: "<<hNameSuffix<<std::endl;
        outputStream<<"         Data: "<<hSoup->Integral(0,hSoup->GetNbinsX()+1)<<std::endl;
        outputStream<<"MC (no Higgs): "<<hMCSum->Integral(0,hMCSum->GetNbinsX()+1)<<std::endl;
        outputStream<<"       W->l: "<<hWJets->Integral(0,hWJets->GetNbinsX()+1)<<std::endl;
        outputStream<<"     TTbarJ: "<<hTTbarJ->Integral(0,hTTbarJ->GetNbinsX()+1)<<std::endl;
        outputStream<<"     TTbarT: "<<hTTbarT->Integral(0,hTTbarT->GetNbinsX()+1)<<std::endl;
        outputStream<<"   single T: "<<hST->Integral(0,hST->GetNbinsX()+1)<<std::endl;
        outputStream<<"   DiBosonJ: "<<hVVJ->Integral(0,hVVJ->GetNbinsX()+1)<<std::endl;
        outputStream<<"   DiBosonT: "<<hVVT->Integral(0,hVVT->GetNbinsX()+1)<<std::endl;
        outputStream<<"        ZTT: "<<hDYJetsZTT->Integral(0,hDYJetsZTT->GetNbinsX()+1)<<std::endl;
        outputStream<<"         ZL: "<<hDYJetsZL->Integral(0,hDYJetsZL->GetNbinsX()+1)<<std::endl;
        outputStream<<"         ZJ: "<<hDYJetsZJ->Integral(0,hDYJetsZJ->GetNbinsX()+1)<<std::endl;
        outputStream<<"Z->ll(m<50): "<<hDYJetsLowM->Integral(0,hDYJetsLowM->GetNbinsX()+1)<<std::endl;
        outputStream<<"  EWK 2Jets: "<<hEWK2Jets->Integral(0,hEWK2Jets->GetNbinsX()+1)<<std::endl;
        outputStream<<"       H(125)->tau tau: "<<hHiggs->Integral(0,hHiggs->GetNbinsX()+1)<<std::endl;
        outputStream<<"     qqH(125)->tau tau: "<<hqqHiggs125->Integral(0,hqqHiggs125->GetNbinsX()+1)<<std::endl;
        outputStream<<"     ggH(125)->tau tau: "<<hggHiggs125->Integral(0,hggHiggs125->GetNbinsX()+1)<<std::endl;
        outputStream<<"      ZH(125)->tau tau: "<<hZHiggs125->Integral(0,hZHiggs125->GetNbinsX()+1)<<std::endl;
        outputStream<<"  WplusH(125)->tau tau: "<<hWplusHiggs125->Integral(0,hWplusHiggs125->GetNbinsX()+1)<<std::endl;
        outputStream<<" WminusH(125)->tau tau: "<<hWminusHiggs125->Integral(0,hWminusHiggs125->GetNbinsX()+1)<<std::endl;
        outputStream<<"QCD: "<<hQCD->Integral(0,hQCD->GetNbinsX()+1)<<std::endl;
        outputStream<<"QCD_MC: "<<hQCD_MC->Integral(0,hQCD_MC->GetNbinsX()+1)<<std::endl;
        outputStream<<"Correction factors:"<<std::endl;
        outputStream<<"QCD control to signal: "<<qcdOStoSS.first<<" +- "<<qcdOStoSS.second<<std::endl;
        outputStream<<"         W MC to DATA: "<<dataToMCScale.first<<" +- "<<dataToMCScale.second<<std::endl;
        outputStream<<"----------------------------------------"<<std::endl;

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

        if(varName.find("Phi-")!=std::string::npos) lowEnd = -M_PI;

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

        std::string plotName;
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

        return hs;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getQCDControlToSignal(unsigned int iCategory,
                                                            unsigned int iSystEffect){

        //return std::pair<float,float>(1.0/3.95319,0.0); //TEST

        ///Return fixed values according to SMH2016 TWiki:
        ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMTauTau2016#QCD_background_estimation_in_Lta
        std::pair<float,float> result(1,0);
        if(myCategoryRejester[iCategory]->name().find("0jet")!=std::string::npos) result = std::pair<float,float>(1.0,0.15);
        else if(myCategoryRejester[iCategory]->name().find("boosted")!=std::string::npos) result = std::pair<float,float>(1.15,0.15*1.15);
        else if(myCategoryRejester[iCategory]->name().find("vbf")) result = std::pair<float,float>(1.2,0.30*1.2);
        else result = std::pair<float,float>(1.0,0.15);

        if(iSystEffect==(unsigned int)HTTAnalysis::QCDSFUp) result.first+=result.second;
        if(iSystEffect==(unsigned int)HTTAnalysis::QCDSFDown) result.first-=result.second;

        std::string varName = "MassVis";

        unsigned int iCategoryNum = myCategoryRejester[iCategory]->qcdSFDenominator()->id();
        TH1F *hSoupLoose = get1DHistogram(iCategoryNum, varName+"Data", iSystEffect);
        if(!hSoupLoose) return result; //MuTau has fixed QCD control to signal transfer factors. In MuTau this category is not defined

        TH1F *hMCSumLoose = getMCSum(iCategoryNum, varName, iSystEffect);
        if(hMCSumLoose) hSoupLoose->Add(hMCSumLoose, -1);

        unsigned int iCategoryDenom = myCategoryRejester[iCategory]->qcdSFNumerator()->id();
        TH1F *hSoupTight = get1DHistogram(iCategoryDenom, varName+"Data", iSystEffect);
        TH1F *hMCSumTight = getMCSum(iCategoryDenom, varName, iSystEffect);
        if(hMCSumTight && hSoupTight) hSoupTight->Add(hMCSumTight, -1);

        double sumTightErr = 0;
        double sumTight = hSoupTight->IntegralAndError(0,hSoupTight->GetNbinsX()+2,sumTightErr);

        double sumLooseErr = 0;
        double sumLoose = hSoupLoose->IntegralAndError(0,hSoupLoose->GetNbinsX()+2,sumLooseErr);

        hSoupTight->Divide(hSoupLoose);
        double ratioErr = (sumTight/sumLoose)*sqrt( (sumTightErr*sumTightErr)/(sumTight*sumTight) + (sumLooseErr*sumLooseErr)/(sumLoose*sumLoose));

        //funtion fitting
        TF1 *line=new TF1("line","[0]",0,350);
        line->SetParameter(0,sumTight/sumLoose);
        TCanvas* c = new TCanvas("QCD_LooseToTight","QCD_LooseToTight",460,500);
        hSoupLoose->SetLineWidth(3);
        hSoupLoose->GetYaxis()->SetTitleOffset(1.4);
        hSoupLoose->GetYaxis()->SetTitle("Tight/Loose");
        hSoupLoose->GetXaxis()->SetTitle("mass");
        gStyle->SetOptStat(11);
        gStyle->SetOptFit(11);
        hSoupTight->Draw();
        //hSoupTight->Fit("line","","",100,250);

        std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
        std::string hNameSuffix = "_"+categoryName;
        std::string plotName = varName+"_"+hNameSuffix+"_"+systEffectName;
        c->Print(TString::Format("fig_png/%s.png",plotName.c_str()).Data());
        c->Print(TString::Format("fig_C/%s.C",plotName.c_str()).Data());

        float param, dparam;
        param=line->GetParameter(0);
        dparam=line->GetParError(0);
/*
        std::cout<<"Category: "<<categoryName<<std::endl
                 <<"\tQCD Tight/Loose` ratio: "<<std::endl
                 <<"\tRatio: "<<sumTight/sumLoose<<" +- "<<ratioErr<<std::endl
                 <<"\tFit: "<<param<<" +- "<<dparam<<std::endl;
 */
        result = std::make_pair(sumTight/sumLoose,ratioErr);
        return result;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F* HTTHistograms::getQCDbackground(unsigned int iCategory,
                                      std::string varName,
                                      unsigned int iSystEffect){

        std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string hName = "h1D" + varName+"QCDEstimate_"+categoryName;
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
        hName+=systEffectName;

        float qcdScale = getQCDControlToSignal(iCategory, iSystEffect).first;
        iCategory = myCategoryRejester[iCategory]->qcdEstimate()->id();

        TH1F *hMCSum = getMCSum(iCategory, varName, iSystEffect);
        TH1F *hSoup = get1DHistogram(iCategory, varName+"Data", iSystEffect);
        if(!hSoup) return 0;
        if(hMCSum) hSoup->Add(hMCSum,-1);

        hSoup->SetName(hName.c_str());
        hSoup->Scale(qcdScale);
        return hSoup;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
std::pair<float,float> HTTHistograms::getWNormalisation(unsigned int iCategory, unsigned int iSystEffect){

        iCategory = myCategoryRejester[iCategory]->wSF()->id();
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
	std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string hNameSuffix =  "_"+categoryName+systEffectName;

        std::string varName = "MassTrans";
        std::string hName = "h1D" + varName;

        TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix));
        if(!hWJets) return std::make_pair(1.0, 0.0);
        TH1F *hMCSum = getMCSum(iCategory, varName, iSystEffect);
        TH1F *hQCD = (TH1F*)getQCDbackground(iCategory,varName, iSystEffect);
        TH1F *hSoup = get1DHistogram(iCategory,varName+"Data", iSystEffect);
        float lumi = HTTAnalysis::getLumi();

        ///Normalise MC histograms according to cross sections
        std::string sampleName = "WJets";
        float weight = getSampleNormalisation(sampleName);
        float scale = weight*lumi;
        hWJets->Scale(scale);

        if(hMCSum) hSoup->Add(hMCSum,-1);
        if(hQCD) hSoup->Add(hQCD,-1);
        hSoup->Add(hWJets);

        float inthWJets=hWJets->Integral(0,hWJets->GetNbinsX()+1);
        float intdata=hSoup->Integral(0,hSoup->GetNbinsX()+1);
        weight=intdata/inthWJets;

        float WSFUncertainty = 0.0;
        if(categoryName.find("0jet")!=std::string::npos) WSFUncertainty = 0.1;
        if(categoryName.find("boosted")!=std::string::npos) WSFUncertainty = 0.1;
        if(categoryName.find("vbf")!=std::string::npos) WSFUncertainty = 0.3;

        if(iSystEffect==(unsigned int)HTTAnalysis::WSFUp) weight*=(1+WSFUncertainty);
        if(iSystEffect==(unsigned int)HTTAnalysis::WSFDown) weight*=(1-WSFUncertainty);

	outputStream<<"WJets scale with uncertainty: "<<weight<<" WSFUncertainty "<<WSFUncertainty
        <<" iSystEffect: "<<iSystEffect
	<<" HTTAnalysis::WSFUp "<<HTTAnalysis::WSFUp
	<<" HTTAnalysis::WSFDown "<<HTTAnalysis::WSFDown
	<<endl;

        return std::make_pair(weight, WSFUncertainty);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::get1DHistogram(unsigned int iCategory, std::string varName,
                                    unsigned int iSystEffect){

        std::string hName = "h1D" + varName;
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
	std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string hNameSuffix =  "_"+categoryName+systEffectName;

        return get1DHistogram(hName+hNameSuffix);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
TH1F *HTTHistograms::getMCSum(unsigned int iCategory, std::string varName,
                              unsigned int iSystEffect,
                              bool sumForW){

        std::string hName = "h1D" + varName;
        std::string systEffectName = HTTAnalysis::systEffectName(iCategory, iSystEffect, myCategoryRejester);
	std::string categoryName = myCategoryRejester[iCategory]->name();
        std::string hNameSuffix =  "_"+categoryName+systEffectName;

        TH1F *hWJets = get1D_WJet_Histogram((hName+"WJets"+hNameSuffix));
        TH1F *hDYJetsLowM = get1D_DYJet_Histogram((hName+"DYLowM"+hNameSuffix));
        TH1F *hDYJets = get1D_DYJet_Histogram((hName+"DYJets"+hNameSuffix));
        TH1F *hTTbar = get1D_TT_Histogram((hName+"TTbar"+hNameSuffix));
        TH1F *hST = get1D_ST_Histogram((hName+"ST"+hNameSuffix));
        TH1F *hVV = get1D_VV_Histogram((hName+"DiBoson"+hNameSuffix));
        TH1F *hEWK2Jets = get1D_EWK2JetsSum(hName+"EWK2Jets"+hNameSuffix);

        ///Protection against null pointers
        ///Null pointers happen when sample was not read, or there were no
        ///events passing particular selection.
        TH1F *hEmpty = 0;
        if(hDYJets) hEmpty = (TH1F*)hDYJets->Clone("hEmpty");
        else if(hDYJetsLowM) hEmpty = (TH1F*)hDYJetsLowM->Clone("hEmpty");
        else if(hWJets) hEmpty = (TH1F*)hWJets->Clone("hEmpty");
        else if(hTTbar) hEmpty = (TH1F*)hTTbar->Clone("hEmpty");
        else if(hST) hEmpty = (TH1F*)hST->Clone("hEmpty");
        else if(hVV) hEmpty = (TH1F*)hVV->Clone("hEmpty");
        else if(hEWK2Jets) hEmpty = (TH1F*)hEWK2Jets->Clone("hEmpty");
        else return 0;
        hEmpty->Reset();

        if(!hWJets) hWJets = (TH1F*)hEmpty->Clone((hName+"WJets"+hNameSuffix).c_str());
        if(!hDYJetsLowM) hDYJetsLowM = (TH1F*)hEmpty->Clone((hName+"DYLowM"+hNameSuffix).c_str());
        if(!hDYJets) hDYJets = (TH1F*)hEmpty->Clone((hName+"DYJets"+hNameSuffix).c_str());
        if(!hTTbar) hTTbar = (TH1F*)hEmpty->Clone((hName+"TTbar"+hNameSuffix).c_str());
        if(!hST) hST = (TH1F*)hEmpty->Clone((hName+"ST"+hNameSuffix).c_str());
        if(!hVV) hVV = (TH1F*)hEmpty->Clone((hName+"DiBoson"+hNameSuffix).c_str());
        if(!hEWK2Jets) hEWK2Jets = (TH1F*)hEmpty->Clone((hName+"EWK2Jets"+hNameSuffix).c_str());
        //////////////////////////////////////////////////////////////////////
        float lumi = HTTAnalysis::getLumi();
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
        std::pair<float, float> wselCorrection(1,1);
        //FIXME if(!sumForW) wselCorrection = getWNormalisation(iCategory, iSystEffect);
        scale = getSampleNormalisation(sampleName)*lumi*wselCorrection.first;
        hWJets->Scale(scale);

        sampleName = "TTbar";
        scale = lumi;
        hTTbar->Scale(scale);

        sampleName = "ST";
        scale = lumi;
        hST->Scale(scale);

        sampleName = "DiBoson";
        scale = lumi;
        hVV->Scale(scale);

        ///EWK 2Jets scaled for preselection during stiching step
        sampleName = "EWK2Jets";
        scale = lumi;
        hEWK2Jets->Scale(scale);

        hName = "h1D" + varName+"MCSum_"+categoryName;
        hName+=systEffectName;

        hEmpty->SetName(hName.c_str());
        hEmpty->Add(hWJets);
        hEmpty->Add(hDYJetsLowM);
        hEmpty->Add(hDYJets);
        hEmpty->Add(hTTbar);
        hEmpty->Add(hST);
        hEmpty->Add(hVV);
        hEmpty->Add(hEWK2Jets);

        return hEmpty;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
