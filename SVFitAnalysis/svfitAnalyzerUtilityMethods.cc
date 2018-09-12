#include <sstream>
#include <bitset>

#include "svfitAnalyzer.h"
#include "svfitHistograms.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string svfitAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

        return getSampleNameFromFileName(myEventProxy);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string svfitAnalyzer::getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        std::string sampleName = "Unknown";

        if(fileName.find("W1JetsToLNu")!=std::string::npos) sampleName = "W1Jets";
        else if(fileName.find("W2JetsToLNu")!=std::string::npos) sampleName = "W2Jets";
        else if(fileName.find("W3JetsToLNu")!=std::string::npos) sampleName = "W3Jets";
        else if(fileName.find("W4JetsToLNu")!=std::string::npos) sampleName = "W4Jets";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) sampleName = "W0Jets";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) sampleName = "W1JetsIncl";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) sampleName = "W2JetsIncl";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) sampleName = "W3JetsIncl";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) sampleName = "W4JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>=0) sampleName = "WAllJets";

        else if(fileName.find("Run201")!=std::string::npos) sampleName =  "Data";
        else if(fileName.find("SUSYGluGluToHToTauTau")!=std::string::npos) sampleName =  "ATT";

        else if(fileName.find("GluGluHToTauTauM110")!=std::string::npos) sampleName =  "ggHTT110";
        else if(fileName.find("GluGluHToTauTauM120")!=std::string::npos) sampleName =  "ggHTT120";
        else if(fileName.find("GluGluHToTauTauM125")!=std::string::npos) sampleName =  "ggHTT125";
        else if(fileName.find("GluGluHToTauTauM130")!=std::string::npos) sampleName =  "ggHTT130";
        else if(fileName.find("GluGluHToTauTauM140")!=std::string::npos) sampleName =  "ggHTT140";

        else if(fileName.find("SUSYGluGluToBBHToTauTauM140")!=std::string::npos) sampleName =  "ggHTT140";
	else if(fileName.find("SUSYGluGluToBBHToTauTauM200")!=std::string::npos) sampleName =  "ggHTT200";
        else if(fileName.find("SUSYGluGluToBBHToTauTauM250")!=std::string::npos) sampleName =  "ggHTT250";

        else if(fileName.find("VBFHToTauTauM110")!=std::string::npos) sampleName =  "qqHTT110";
        else if(fileName.find("VBFHToTauTauM120")!=std::string::npos) sampleName =  "qqHTT120";
        else if(fileName.find("VBFHToTauTauM125")!=std::string::npos) sampleName =  "qqHTT125";
        else if(fileName.find("VBFHToTauTauM130")!=std::string::npos) sampleName =  "qqHTT130";
        else if(fileName.find("VBFHToTauTauM140")!=std::string::npos) sampleName =  "qqHTT140";

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
        else if(fileName.find("WZJToLLLNu")!=std::string::npos) sampleName =  "WZJToLLLNu";
        else if(fileName.find("WWTo1L1Nu2Q")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";
        else if(fileName.find("WWToLNuQQ")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";
        else if(fileName.find("WZTo1L1Nu2Q")!=std::string::npos) sampleName =  "WZTo1L1Nu2Q";
        else if(fileName.find("VVTo2L2Nu")!=std::string::npos) sampleName =  "VVTo2L2Nu";
        else if(fileName.find("WZTo2L2Q")!=std::string::npos) sampleName =  "WZTo2L2Q";
        else if(fileName.find("EWKWMinus")!=std::string::npos) sampleName =  "EWKWMinus";
        else if(fileName.find("EWKWPlus")!=std::string::npos) sampleName =  "EWKWPlus";
        else if(fileName.find("EWKZ2JetsZToLL")!=std::string::npos) sampleName =  "EWKZ2JetsZToLL";
        else if(fileName.find("EWKZ2JetsZToNuNu")!=std::string::npos) sampleName =  "EWKZ2JetsZToNuNu";
        else if(fileName.find("QCD")!=std::string::npos) sampleName =  "QCD_MC";
        else if(fileName.find("DY")!=std::string::npos) sampleName =  getDYSampleName(myEventProxy);
        else if(fileName.find("TTTune")!=std::string::npos) sampleName =  "TTbar";
        std::string matchingMode = getMatchingName(myEventProxy);

        if(sampleName=="Unknown") std::cout<<"Unkwown sample type. "<<fileName<<std::endl;

        sampleName += matchingMode;

        return sampleName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string svfitAnalyzer::getDYSampleName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        std::string jetsName = "";
        if(fileName.find("DY1JetsToLLM50")!=std::string::npos) jetsName ="1Jets";
        else if(fileName.find("DY2JetsToLLM50")!=std::string::npos) jetsName = "2Jets";
        else if(fileName.find("DY3JetsToLLM50")!=std::string::npos) jetsName = "3Jets";
        else if(fileName.find("DY4JetsToLLM50")!=std::string::npos) jetsName = "4Jets";
        else if(fileName.find("DYJetsToLLM10to50")!=std::string::npos) jetsName =  "LowM";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) jetsName =  "0Jets";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) jetsName =  "1JetsIncl";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) jetsName =  "2JetsIncl";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) jetsName =  "3JetsIncl";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) jetsName =  "4JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos) jetsName =  "AllJets";

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
	if(fileName.find("MM_")!=std::string::npos) {
                if(leg1MCMatch==2 && leg2MCMatch==2) decayName = "L";
        }
        return "DY"+jetsName+"Match"+decayName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string svfitAnalyzer::getMatchingName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        std::vector<std::string> sampleNames = {"TTTune", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};

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
