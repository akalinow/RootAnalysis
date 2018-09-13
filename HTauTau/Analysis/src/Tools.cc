#include <string>
#include <vector>
#include "Tools.h"
#include "AnalysisEnums.h"

namespace HTTAnalysis {

float getLumi(){
  /*
        //./.local/bin/brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i processedLumis_SingleMuon.json       
	float run2016 = 35.87*1E3; //Updated Run2016 luminosity
        return run2016; //pb-1 data for NTUPLES_05_12_2016
  */
  
	/* Run2017X-31Mar2018-v1/MINIAOD */
        float run2017 = 41.501*1E3; //pb-1 data for NTUPLES_12_07_2018
		        
	return run2017;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float getCrossSection(const std::string & sampleName){

        float crossSection = 0;

       ///Cross sections taken from
        if(sampleName=="DYLowM") {
              //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
                crossSection = 18610;
        }
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
        if(sampleName.find("DYJetsMatch")!=std::string::npos || sampleName=="DYJets") {
                //xsection for 3xZ->mu mu M50 in [pb]
                crossSection = 3*1921.8;
        }
        if(sampleName.find("WJets")!=std::string::npos) {
                //xsection for 3xW->mu nu in [pb]
                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
                crossSection = 3*20508.9;
        }
        if(sampleName.find("TTbar")!=std::string::npos) {
                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
	  crossSection = 831.76;
        }
	if(sampleName.find("TTTo2L2Nu")!=std::string::npos) {
	  ///Inclusive TT*BR(WW->2L2Nu)
	  crossSection = 831.76*0.1086*0.1086;
	}
	if(sampleName.find("TTToHadronic")!=std::string::npos) {
	  ///Inclusive TT*BR(WW->hadronic)
	  crossSection = 831.76*(1-0.1086)*(1-0.1086);
	}
	if(sampleName.find("TTToSemiLeptonic")!=std::string::npos) {
	  ///Inclusive TT*BR(WW->SemiLeptonic)
	  crossSection = 831.76*(1-0.1086)*0.1086*2;
	}
        //https://twiki.cern.ch/twiki/pub/LHCPhysics/LHCHXSWG/Higgs_XSBR_YR4_update.xlsx
        //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2014#Higgs_2_fermions
        //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014
        //Xsection for mass!=125 are calculated using luminosity ratio, and cross section for 8 TeV

        if(sampleName=="ggHTT110") crossSection = 5.790E+01*7.91E-02; //Analysis Note
        if(sampleName=="ggHTT120") crossSection = 5.222E+01*6.98E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ggHTT125") crossSection = 4.858E+01*6.27E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ggHTT130") crossSection = 4.531E+01*5.41E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ggHTT140") crossSection = 3.600E+01*3.60E-02; //Analysis Note

        if(sampleName=="qqHTT110") crossSection = 4.434E+00*7.95E-02; //Analysis Note
        if(sampleName=="qqHTT120") crossSection = 3.935E+00*6.98E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="qqHTT125") crossSection = 3.782E+00*6.27E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="qqHTT130") crossSection = 3.637E+00*5.41E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="qqHTT140") crossSection = 3.492E+00*3.60E-02; //Analysis Note

        ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
        if(sampleName=="WplusHTT110") crossSection = 1.335*7.91E-02; //Analysis Note
        if(sampleName=="WplusHTT120") crossSection = 0.9558*0.0698;
        if(sampleName=="WplusHTT125") crossSection = 0.8400*0.0627;
        if(sampleName=="WplusHTT130") crossSection = 0.7414*0.0541;
        if(sampleName=="WplusHTT140") crossSection = 0.6308*3.60E-02; //Analysis Note

        if(sampleName=="WminusHTT110") crossSection = 0.8587*7.91E-02; //Analysis Note
        if(sampleName=="WminusHTT120") crossSection = 0.6092*0.0698;
        if(sampleName=="WminusHTT125") crossSection = 0.5328*0.0627;
        if(sampleName=="WminusHTT130") crossSection = 0.4676*0.0541;
        if(sampleName=="WminusHTT140") crossSection = 0.3940*3.60E-02; //Analysis Note

        if(sampleName=="ZHTT110") crossSection = 1.309*7.15E-02; //Analysis Note
        if(sampleName=="ZHTT120") crossSection = 0.994*0.0698;
        if(sampleName=="ZHTT125") crossSection = 0.884*0.0627;
        if(sampleName=="ZHTT130") crossSection = 0.790*0.0541;
        if(sampleName=="ZHTT140") crossSection = 0.6514*3.60E-02; //Analysis Note

	if(sampleName=="ttHTT110") crossSection = 0E+00*7.95E-02; //FIXME
        if(sampleName=="ttHTT120") crossSection = 0.5697E+00*6.98E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ttHTT125") crossSection = 0.5071E+00*6.27E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ttHTT130") crossSection = 0.4539E+00*5.41E-02; //CERNYellowReportPageAt13TeV*CERNYellowReportPageBR
        if(sampleName=="ttHTT140") crossSection = 0E+00*3.60E-02; //FIXME

        ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
        if(sampleName.find("ZZTo2L2Q")!=std::string::npos) crossSection = 3.22;
        if(sampleName.find("ZZTo4L")!=std::string::npos) crossSection = 1.212;
        if(sampleName.find("WZTo1L3Nu")!=std::string::npos) crossSection = 3.05;
	if(sampleName.find("WZTo3LNu")!=std::string::npos) crossSection = 3.05/0.2*3*0.0337;//WZTo1L3Nu/BR(Z->nunu)*BR(Z->ll)	
        if(sampleName.find("WZJToLLLNu")!=std::string::npos) crossSection = 4.708;
        if(sampleName.find("WWTo1L1Nu2Q")!=std::string::npos) crossSection = 49.997;
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

        return crossSection;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string systEffectName(unsigned int iSystEffect){
        if(iSystEffect==(int)NOMINAL) return "";
        else if(iSystEffect==(int)TESUp) return "_CMS_scale_t_mt_13TeVUp";
        else if(iSystEffect==(int)TESDown) return "_CMS_scale_t_mt_13TeVDown";
        else if(iSystEffect==(int)JESUp) return "_CMS_scale_j_13TeVUp";
        else if(iSystEffect==(int)JESDown) return "_CMS_scale_j_13TeVDown";
        else if(iSystEffect==(int)M2TUp) return "_CMS_htt_ZLShape_mt_13TeVUp";
        else if(iSystEffect==(int)M2TDown) return "_CMS_htt_ZLShape_mt_13TeVDown";
        else if(iSystEffect==(int)E2TUp) return "_CMS_htt_ZEShape_et_13TeVUp";
        else if(iSystEffect==(int)E2TDown) return "_CMS_htt_ZEShape_et_13TeVDown";
        else if(iSystEffect==(int)J2TUp) return "_CMS_htt_jetToTauFake_13TeVUp";
        else if(iSystEffect==(int)J2TDown) return "_CMS_htt_jetToTauFake_13TeVDown";
        else if(iSystEffect==(int)ZPtUp) return "_CMS_htt_dyShape_13TeVUp";
        else if(iSystEffect==(int)ZPtDown) return "_CMS_htt_dyShape_13TeVDown";
        else if(iSystEffect==(int)TTUp) return "_CMS_htt_ttbarShape_13TeVUp";
        else if(iSystEffect==(int)TTDown) return "_CMS_htt_ttbarShape_13TeVDown";
        else if(iSystEffect==(int)QCDSFUp) return "_QCDSFUncert_mt_CAT_13TeVUp";
        else if(iSystEffect==(int)QCDSFDown) return "_QCDSFUncert_mt_CAT_13TeVDown";
        else if(iSystEffect==(int)WSFUp) return "_WSFUncert_mt_CAT_13TeVUp";
        else if(iSystEffect==(int)WSFDown) return "_WSFUncert_mt_CAT_13TeVDown";
        else if(iSystEffect==(int)ggUp) return "_CMS_scale_gg_13TeVUp";
        else if(iSystEffect==(int)ggDown) return "_CMS_scale_gg_13TeVDown";
        else if(iSystEffect==(int)ZmumuUp) return "_CMS_htt_zmumuShape_CAT_13TeVUp";
        else if(iSystEffect==(int)ZmumuDown) return "_CMS_htt_zmumuShape_CAT_13TeVDown";
        else return "_UnknownSystEffect";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string systEffectName(unsigned int iCategory, unsigned int iSystEffect,
                            const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester){

        std::string systEffectName = HTTAnalysis::systEffectName(iSystEffect);
        if(systEffectName.find("CAT")!=std::string::npos) {
                std::string categoryName = aCategoryRejester[iCategory]->name();
                systEffectName.replace(systEffectName.find("CAT"),3,categoryName);
        }
        return systEffectName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getSampleName(const EventProxyHTT & myEventProxy){

        return getSampleNameFromFileName(myEventProxy);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        std::string sampleName = "Unknown";

        if(fileName.find("W1JetsToLNu")!=std::string::npos) sampleName = "W1Jets";
        else if(fileName.find("W2JetsToLNu")!=std::string::npos) sampleName = "W2Jets";
        else if(fileName.find("W3JetsToLNu")!=std::string::npos) sampleName = "W3Jets";
        else if(fileName.find("W4JetsToLNu")!=std::string::npos) sampleName = "W4Jets";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) sampleName = "W0Jets";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) sampleName = "W1JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) sampleName = "W2JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) sampleName = "W3JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) sampleName = "W4JetsIncl";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) sampleName = "WAllJets";

        else if(fileName.find("SingleMuonRun201")!=std::string::npos) sampleName =  "Data";//FIX!
        else if(fileName.find("SUSYGluGluToHToTauTau")!=std::string::npos) sampleName =  "ATT";

        else if(fileName.find("GluGluHToTauTauM110")!=std::string::npos) sampleName =  "ggHTT110";
        else if(fileName.find("GluGluHToTauTauM120")!=std::string::npos) sampleName =  "ggHTT120";
        else if(fileName.find("GluGluHToTauTauM125")!=std::string::npos) sampleName =  "ggHTT125";
        else if(fileName.find("GluGluHToTauTauM130")!=std::string::npos) sampleName =  "ggHTT130";
        else if(fileName.find("GluGluHToTauTauM140")!=std::string::npos) sampleName =  "ggHTT140";

        else if(fileName.find("VBFHToTauTauM110")!=std::string::npos) sampleName =  "qqHTT110";
        else if(fileName.find("VBFHToTauTauM120")!=std::string::npos) sampleName =  "qqHTT120";
        else if(fileName.find("VBFHToTauTauM125")!=std::string::npos) sampleName =  "qqHTT125";
        else if(fileName.find("VBFHToTauTauM130")!=std::string::npos) sampleName =  "qqHTT130";
        else if(fileName.find("VBFHToTauTauM140")!=std::string::npos) sampleName =  "qqHTT140";

	else if(fileName.find("ttHToTauTauM110")!=std::string::npos) sampleName =  "ttHTT110";
        else if(fileName.find("ttHToTauTauM120")!=std::string::npos) sampleName =  "ttHTT120";
        else if(fileName.find("ttHToTauTauM125")!=std::string::npos) sampleName =  "ttHTT125";
        else if(fileName.find("ttHToTauTauM130")!=std::string::npos) sampleName =  "ttHTT130";
        else if(fileName.find("ttHToTauTauM140")!=std::string::npos) sampleName =  "ttHTT140";

        else if(fileName.find("WplusHToTauTauM110")!=std::string::npos) sampleName =  "WplusHTT110";
        else if(fileName.find("WplusHToTauTauM120")!=std::string::npos) sampleName =  "WplusHTT120";
        else if(fileName.find("WplusHToTauTauM125")!=std::string::npos) sampleName =  "WplusHTT125";
        else if(fileName.find("WplusHToTauTauM130")!=std::string::npos) sampleName =  "WplusHTT130";
        else if(fileName.find("WplusHToTauTauM140")!=std::string::npos) sampleName =  "WplusHTT140";

        else if(fileName.find("WminusHToTauTauM110")!=std::string::npos) sampleName =  "WminusHTT110";
        else if(fileName.find("WminusHToTauTauM120")!=std::string::npos) sampleName =  "WminusHTT120";
        else if(fileName.find("WminusHToTauTauM125")!=std::string::npos) sampleName =  "WminusHTT125";
        else if(fileName.find("WminusHToTauTauM130")!=std::string::npos) sampleName =  "WminusHTT130";
        else if(fileName.find("WminusHToTauTauM140")!=std::string::npos) sampleName =  "WminusHTT140";

        else if(fileName.find("ZHToTauTauM110")!=std::string::npos) sampleName =  "ZHTT110";
        else if(fileName.find("ZHToTauTauM120")!=std::string::npos) sampleName =  "ZHTT120";
        else if(fileName.find("ZHToTauTauM125")!=std::string::npos) sampleName =  "ZHTT125";
        else if(fileName.find("ZHToTauTauM130")!=std::string::npos) sampleName =  "ZHTT130";
        else if(fileName.find("ZHToTauTauM140")!=std::string::npos) sampleName =  "ZHTT140";

        else if(fileName.find("STtWantitop")!=std::string::npos) sampleName =  "Wantitop";
        else if(fileName.find("STtWtop")!=std::string::npos) sampleName =  "Wtop";
        else if(fileName.find("STtchannel__antitop")!=std::string::npos) sampleName =  "t-channel_antitop";
        else if(fileName.find("STtchannel__top")!=std::string::npos) sampleName =  "t-channel_top";
        else if(fileName.find("ZZTo2L2Q")!=std::string::npos) sampleName =  "ZZTo2L2Q";
        else if(fileName.find("ZZTo4L")!=std::string::npos) sampleName =  "ZZTo4L";
        else if(fileName.find("WZTo1L3Nu")!=std::string::npos) sampleName =  "WZTo1L3Nu";
	else if(fileName.find("WZTo3LNu")!=std::string::npos) sampleName =  "WZTo3LNu";	
        else if(fileName.find("WZJToLLLNu")!=std::string::npos) sampleName =  "WZJToLLLNu";
        else if(fileName.find("WWTo1L1Nu2Q")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";
        else if(fileName.find("WWToLNuQQ")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";	
	///else if(fileName.find("WWTo2L2Nu")!=std::string::npos) sampleName =  "WWTo2L2Nu";	
        else if(fileName.find("WZTo1L1Nu2Q")!=std::string::npos) sampleName =  "WZTo1L1Nu2Q";
        else if(fileName.find("VVTo2L2Nu")!=std::string::npos) sampleName =  "VVTo2L2Nu";
        else if(fileName.find("WZTo2L2Q")!=std::string::npos) sampleName =  "WZTo2L2Q";
        else if(fileName.find("EWKWMinus")!=std::string::npos) sampleName =  "EWKWMinus";
        else if(fileName.find("EWKWPlus")!=std::string::npos) sampleName =  "EWKWPlus";
        else if(fileName.find("EWKZ2JetsZToLL")!=std::string::npos) sampleName =  "EWKZ2JetsZToLL";
        else if(fileName.find("EWKZ2JetsZToNuNu")!=std::string::npos) sampleName =  "EWKZ2JetsZToNuNu";
        else if(fileName.find("QCD")!=std::string::npos) sampleName =  "QCD_MC";
        else if(fileName.find("DY")!=std::string::npos) sampleName =  getDYSampleName(myEventProxy);
	else if(fileName.find("TTTo2L2Nu")!=std::string::npos) sampleName =  "TTTo2L2Nu";
	else if(fileName.find("TTToHadronic")!=std::string::npos) sampleName =  "TTToHadronic";
	else if(fileName.find("TTToSemiLeptonic")!=std::string::npos) sampleName =  "TTToSemiLeptonic";
	else if(fileName.find("TTTune")!=std::string::npos) sampleName =  "TTbar";
	
        std::string matchingMode = getMatchingName(myEventProxy);

        if(sampleName=="Unknown") std::cout<<"Unkwown sample type. "<<fileName<<std::endl;

        sampleName += matchingMode;
 
        return sampleName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getDYSampleName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        std::string jetsName = "";
        if(fileName.find("DY1JetsToLLM50")!=std::string::npos) jetsName ="1Jets";
        else if(fileName.find("DY2JetsToLLM50")!=std::string::npos) jetsName = "2Jets";
        else if(fileName.find("DY3JetsToLLM50")!=std::string::npos) jetsName = "3Jets";
        else if(fileName.find("DY4JetsToLLM50")!=std::string::npos) jetsName = "4Jets";
        else if(fileName.find("DYJetsToLLM10to50")!=std::string::npos) jetsName =  "LowM";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) jetsName =  "0Jets";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) jetsName =  "1JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) jetsName =  "2JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) jetsName =  "3JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) jetsName =  "4JetsIncl";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) jetsName =  "AllJets";

        //int decayModeBoson = myEventProxy.event->getDecayModeBoson();
        int leg1MCMatch = 6, leg2MCMatch = 6;
        if(myEventProxy.pairs->size()) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg1 = aPair.getLeg1();
                HTTParticle aLeg2 = aPair.getLeg2();
                leg1MCMatch = aLeg1.getProperty(PropertyEnum::mc_match);
                leg2MCMatch = aLeg2.getProperty(PropertyEnum::mc_match);
        }
        std::string decayName = "Unknown";
        if(fileName.find("MT_")!=std::string::npos) {
                if(leg2MCMatch<5) decayName = "L";
                else if(leg2MCMatch==5) decayName = "T";
                else decayName = "J";
        }
        if(fileName.find("TT_")!=std::string::npos) {
                if(leg1MCMatch==5 && leg2MCMatch==5) decayName = "T";
                else if(leg1MCMatch<6 && leg2MCMatch<6) decayName = "L";
                else decayName = "J";
        }
        return "DY"+jetsName+"Match"+decayName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getMatchingName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        std::vector<std::string> sampleNames = {"TTTune", "TTTo2L2Nu", "TTToHadronic", "TTToSemiLeptonic", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZTo3LNu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};

        bool sampleToAnalyze = false;
        for(auto sampleName : sampleNames) {
                if(fileName.find(sampleName)!=std::string::npos) sampleToAnalyze = true;
        }
        if(!sampleToAnalyze) return "";
        if(!(fileName.find("TT_")!=std::string::npos || fileName.find("MT_")!=std::string::npos)) return "";

        int tauMCMatch_1 = 6, tauMCMatch_2 = 6;
        if(myEventProxy.pairs->size() && fileName.find("MT_")!=std::string::npos) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg2 = aPair.getTau();
                tauMCMatch_1 = aLeg2.getProperty(PropertyEnum::mc_match);
                if(tauMCMatch_1 == 5) return "MatchT"; else return "MatchJ";
        }

        if(myEventProxy.pairs->size() && fileName.find("TT_")!=std::string::npos) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg1 = aPair.getLeg1(), aLeg2 = aPair.getLeg2();
                tauMCMatch_1 = aLeg1.getProperty(PropertyEnum::mc_match);
                tauMCMatch_2 = aLeg2.getProperty(PropertyEnum::mc_match);
                if(tauMCMatch_1 == 5 && tauMCMatch_2 == 5) return "MatchT"; else return "MatchJ";
        }

        return "";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<std::string> getTauDecayName(int decModeMinus, int decModePlus){

        std::vector<std::string> types;

        if(decModeMinus==HTTAnalysis::tauDecay1ChargedPion0PiZero && decModePlus==HTTAnalysis::tauDecay1ChargedPion0PiZero) types.push_back("PiPi0Pi0");

        if(HTTAnalysis::isOneProng(decModeMinus) && HTTAnalysis::isOneProng(decModePlus) ) types.push_back("1Prong1Prong");

        if( (decModeMinus==HTTAnalysis::tauDecay1ChargedPion0PiZero && isLepton(decModePlus) ) ||
            (HTTAnalysis::isLepton(decModeMinus) && decModePlus==HTTAnalysis::tauDecay1ChargedPion0PiZero)) types.push_back("Lepton1Prong0Pi0");

        if( (HTTAnalysis::isOneProng(decModeMinus) && HTTAnalysis::isLepton(decModePlus) ) ||
            (HTTAnalysis::isLepton(decModeMinus) && HTTAnalysis::isOneProng(decModePlus) ) ) types.push_back("Lepton1Prong");

        if(decModeMinus==HTTAnalysis::tauDecay1ChargedPion1PiZero && decModePlus==HTTAnalysis::tauDecay1ChargedPion1PiZero ) types.push_back("PiPlusPiMinus2Pi0");


        if(HTTAnalysis::isOneProng(decModeMinus) && decModeMinus!=HTTAnalysis::tauDecay1ChargedPion0PiZero &&
           HTTAnalysis::isOneProng(decModePlus) && decModePlus!=HTTAnalysis::tauDecay1ChargedPion0PiZero ) types.push_back("1Prong1ProngXPi0");

        if(HTTAnalysis::isLepton(decModePlus) && HTTAnalysis::isLepton(decModeMinus)) types.push_back("LeptonLepton");

        return types;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool isOneProng(int decMode){
        if(decMode==HTTAnalysis::tauDecay1ChargedPion0PiZero ||
           decMode==HTTAnalysis::tauDecay1ChargedPion1PiZero ||
           decMode==HTTAnalysis::tauDecay1ChargedPion2PiZero ||
           decMode==HTTAnalysis::tauDecay1ChargedPion3PiZero ) return true;
        else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool isLepton(int decMode){
        if(decMode==HTTAnalysis::tauDecaysElectron || decMode==HTTAnalysis::tauDecayMuon) return true;
        else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float getGenWeight(const EventProxyHTT & myEventProxy){

        ///MC weights can be quite large, but are always +-const.
        ///to avoid counter overflow we keep only sign.
        return myEventProxy.event->getMCWeight()/fabs(myEventProxy.event->getMCWeight());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
					  const HTTParticle & aLeg1,
					  const HTTParticle & aLeg2, 
					  float deltaR){

        std::vector<HTTParticle> separatedJets;

        for(auto aJet : *myEventProxy.jets) {
                float dRLeg2 = aJet.getP4().DeltaR(aLeg2.getP4());
                float dRLeg1 = aJet.getP4().DeltaR(aLeg1.getP4());
                bool loosePFJetID = aJet.getProperty(PropertyEnum::PFjetID)>=1;
                bool jetEtaCut = std::abs(aJet.getP4().Eta())<4.7;
                if(dRLeg1>deltaR && dRLeg2>deltaR && loosePFJetID && jetEtaCut) separatedJets.push_back(aJet);
        }
        return separatedJets;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////  

}
