#include <string>
#include <vector>
#include "Tools.h"
#include "AnalysisEnums.h"

namespace HTTAnalysis {

float getLumi(){
  /*
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
	run2016 = 35.87*1E3*1E6; //Updated Run2016 luminosity
  
	float run2017BPromptv1 = 4.0*1E3*1E6;
	float run2017BPromptv2 = 0.7*1E3*1E6;
        float run2017CPromptv1 = 1.3*1E3*1E6;
        float run2017CPromptv2 = 4.1*1E3*1E6;
	*/
        float run2017CPromptv3 = 4.3*1E3*1E6;
        float run2017DPromptv1 = 4.3*1E3*1E6;
        float run2017EPromptv1 = 9.1*1E3*1E6;

        float run2017B12SepReReco = 2.1*1E3*1E6;
        float run2017C12SepReReco = 5.4*1E3*1E6;
	/*run2017BPromptv1+run2017BPromptv2*/
        float run2017 = run2017B12SepReReco
                      + run2017C12SepReReco
                      + run2017CPromptv3
                      + run2017DPromptv1
                      + run2017EPromptv1;//data for NTUPLES_08_10_2017
	        
        //return run2016*1E-6; //pb-1 data for NTUPLES_05_12_2016
	return run2017*1E-6; //pb-1 data for NTUPLES_08_10_2017
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
                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/KlubTwikiRun2
                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeVInclusive
                crossSection = 831.76;
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

        ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_and_data_samples
        if(sampleName.find("ZZTo2L2Q")!=std::string::npos) crossSection = 3.22;
        if(sampleName.find("ZZTo4L")!=std::string::npos) crossSection = 1.212;
        if(sampleName.find("WZTo1L3Nu")!=std::string::npos) crossSection = 3.05;
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

}
