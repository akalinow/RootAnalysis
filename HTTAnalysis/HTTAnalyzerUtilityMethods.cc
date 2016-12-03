#include <sstream>
#include <bitset>

#include "RooRealVar.h"

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::getPreselectionEff(const EventProxyHTT & myEventProxy){

  if(true || ntupleFile_!=myEventProxy.getTTree()->GetCurrentFile()){
    ntupleFile_ = myEventProxy.getTTree()->GetCurrentFile();
    if(hStatsFromFile) delete hStatsFromFile;
    hStatsFromFile = (TH1F*)ntupleFile_->Get("hStats");
    
    std::string hName = "h1DStats"+getSampleName(myEventProxy);
    TH1F *hStats = myHistos_->get1DHistogram(hName.c_str(),true);
    
    float genWeight = getGenWeight(myEventProxy);
    
    hStats->SetBinContent(2,std::abs(hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1))*genWeight)); 
    hStats->SetBinContent(3,std::abs(hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3))*genWeight));
  }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){
  
  return getSampleNameFromFileName(myEventProxy);  
  
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

  std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
  std::string sampleName = "Unknown";

  if(fileName.find("W1JetsToLNu")!=std::string::npos) sampleName = "W1Jets";
  else if(fileName.find("W2JetsToLNu")!=std::string::npos) sampleName = "W2Jets";
  else if(fileName.find("W3JetsToLNu")!=std::string::npos) sampleName = "W3Jets";
  else if(fileName.find("W4JetsToLNu")!=std::string::npos) sampleName = "W4Jets";
  else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) sampleName = "W0Jets";
  else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) sampleName = "WAllJets";
  //TEST else if(fileName.find("WJetsToLNu")!=std::string::npos) sampleName = "WJets";
  
  else if(fileName.find("Run201")!=std::string::npos) sampleName = "Data";
  else if(fileName.find("SUSYGluGluToHToTauTau")!=std::string::npos) sampleName = "A";

  else if(fileName.find("GluGluHToTauTauM120")!=std::string::npos) sampleName = "ggH120";
  else if(fileName.find("GluGluHToTauTauM125")!=std::string::npos) sampleName = "ggH125";
  else if(fileName.find("GluGluHToTauTauM130")!=std::string::npos) sampleName = "ggH130";

  else if(fileName.find("VBFHToTauTauM120")!=std::string::npos) sampleName = "qqH120";  
  else if(fileName.find("VBFHToTauTauM125")!=std::string::npos) sampleName = "qqH125";
  else if(fileName.find("VBFHToTauTauM130")!=std::string::npos) sampleName = "qqH130";

  else if(fileName.find("WplusHToTauTau_M120")!=std::string::npos) sampleName = "WplusHM120";  
  else if(fileName.find("WplusHToTauTau_M125")!=std::string::npos) sampleName = "WplusHM125";  
  else if(fileName.find("WplusHToTauTau_M130")!=std::string::npos) sampleName = "WplusHM130";  

  else if(fileName.find("WminusHToTauTau_M120")!=std::string::npos) sampleName = "WminusHM120";  
  else if(fileName.find("WminusHToTauTau_M125")!=std::string::npos) sampleName = "WminusHM125";  
  else if(fileName.find("WminusHToTauTau_M130")!=std::string::npos) sampleName = "WminusHM130";  

  else if(fileName.find("ZHM120")!=std::string::npos) sampleName = "ZHM120";  
  else if(fileName.find("ZHM125")!=std::string::npos) sampleName = "ZHM125";  
  else if(fileName.find("ZHM130")!=std::string::npos) sampleName = "ZHM130";  

  else if(fileName.find("STtWantitop")!=std::string::npos) sampleName = "Wantitop";
  else if(fileName.find("STtWtop")!=std::string::npos) sampleName = "Wtop";
  else if(fileName.find("STtchannel__antitop")!=std::string::npos) sampleName = "t-channel_antitop";
  else if(fileName.find("STtchannel__top")!=std::string::npos) sampleName = "t-channel_top";   
  else if(fileName.find("ZZTo2L2Q")!=std::string::npos) sampleName = "ZZTo2L2Q";
  else if(fileName.find("ZZTo4L")!=std::string::npos) sampleName = "ZZTo4L";
  else if(fileName.find("WZTo1L3Nu")!=std::string::npos) sampleName = "WZTo1L3Nu";
  else if(fileName.find("WZJToLLLNu")!=std::string::npos) sampleName = "WZJToLLLNu";
  else if(fileName.find("WWTo1L1Nu2Q")!=std::string::npos) sampleName = "WWTo1L1Nu2Q";
  else if(fileName.find("WZTo1L1Nu2Q")!=std::string::npos) sampleName = "WZTo1L1Nu2Q";
  else if(fileName.find("VVTo2L2Nu")!=std::string::npos) sampleName = "VVTo2L2Nu";
  else if(fileName.find("WZTo2L2Q")!=std::string::npos) sampleName = "WZTo2L2Q";
  else if(fileName.find("EWKWMinus")!=std::string::npos) sampleName = "EWKWMinus";
  else if(fileName.find("EWKWPlus")!=std::string::npos) sampleName = "EWKWPlus";
  else if(fileName.find("EWKZ2JetsZToLL")!=std::string::npos) sampleName = "EWKZ2JetsZToLL";
  else if(fileName.find("EWKZ2JetsZToNuNu")!=std::string::npos) sampleName = "EWKZ2JetsZToNuNu";
  else if(fileName.find("QCD")!=std::string::npos) sampleName = "QCD_MC";
  else if(fileName.find("DY")!=std::string::npos) sampleName = getDYSampleName(myEventProxy);
  else if(fileName.find("TTTune")!=std::string::npos) sampleName = "TTbar";
  std::string matchingMode = getMatchingName(myEventProxy);
  
  return sampleName + matchingMode;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getDYSampleName(const EventProxyHTT & myEventProxy){

  std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

  std::string jetsName = "";    
  if(fileName.find("DY1JetsToLL")!=std::string::npos) jetsName ="1Jets";
  else if(fileName.find("DY2JetsToLL")!=std::string::npos) jetsName = "2Jets";
  else if(fileName.find("DY3JetsToLL")!=std::string::npos) jetsName = "3Jets";
  else if(fileName.find("DY4JetsToLL")!=std::string::npos) jetsName = "4Jets";
  else if(fileName.find("DYJetsToLLM10to50")!=std::string::npos) jetsName =  "LowM";
  //TEST else if(fileName.find("DYJetsToLLM50")!=std::string::npos) jetsName =  "Jets";  
  else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) jetsName =  "0Jets";
  else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) jetsName =  "AllJets";  
  
  int decayModeBoson = myEventProxy.event->getDecayModeBoson();
  int tauMCMatch = 6;
  if(myEventProxy.pairs->size()){
    HTTPair aPair = (*myEventProxy.pairs)[0];
    HTTParticle aTau = aPair.getTau();
    tauMCMatch = aTau.getProperty(PropertyEnum::mc_match);
  }    
  std::string decayName = "Unknown";
  if(tauMCMatch<5) decayName = "L";
  else if(tauMCMatch==5) decayName = "T";
  else decayName = "J";

  return "DY"+jetsName+"Match"+decayName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getMatchingName(const EventProxyHTT & myEventProxy){

  std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
  std::vector<std::string> sampleNames = {"TTTune", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};

  bool sampleToAnalyze = false;
  for(auto sampleName:sampleNames){
  	if(fileName.find(sampleName)!=std::string::npos) sampleToAnalyze = true;
  	}
  if(!sampleToAnalyze) return "";
  if(!(fileName.find("TT_")!=std::string::npos || fileName.find("MT_")!=std::string::npos)) return "";
  
  int tauMCMatch_1 = 6, tauMCMatch_2 = 6;
  if(myEventProxy.pairs->size() && fileName.find("MT_")!=std::string::npos){
    HTTPair aPair = (*myEventProxy.pairs)[0];
    HTTParticle aTau = aPair.getTau();
    tauMCMatch_1 = aTau.getProperty(PropertyEnum::mc_match);
    if(tauMCMatch_1 == 5) return "MatchT"; else return "MatchJ";
  }
  
  if(myEventProxy.pairs->size() && fileName.find("TT_")!=std::string::npos){
    HTTPair aPair = (*myEventProxy.pairs)[0];
    HTTParticle aTau_1 = aPair.getLeg1(), aTau_2 = aPair.getLeg2();
    tauMCMatch_1 = aTau_1.getProperty(PropertyEnum::mc_match);
    tauMCMatch_2 = aTau_2.getProperty(PropertyEnum::mc_match);
    if(tauMCMatch_1 == 5 && tauMCMatch_2 == 5) return "MatchT"; else return "MatchJ";
  }

  return "";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getPUWeight(const EventProxyHTT & myEventProxy){

  if(getSampleName(myEventProxy)=="Data") return 1.0;

  if(!puDataFile_ || !puMCFile_ ||
     puDataFile_->IsZombie() ||
     puMCFile_->IsZombie()) { return 1.0; }

  if(!hPUVec_.size())  hPUVec_.resize(1);

  if(!hPUVec_[0]){
    std::string hName = "pileup";
    TH1F *hPUData = (TH1F*)puDataFile_->Get(hName.c_str());
    TH1F *hPUSample = (TH1F*)puMCFile_->Get(hName.c_str());
    ///Normalise both histograms.
    hPUData->Scale(1.0/hPUData->Integral(0,hPUData->GetNbinsX()+1));
    hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
    ///
    hPUData->SetDirectory(0);
    hPUSample->SetDirectory(0);
    hPUData->Divide(hPUSample);
    hPUData->SetName(("h1DPUWeight"+getSampleName(myEventProxy)).c_str());
    hPUVec_[0] =  hPUData;
  }

  int iBinPU = hPUVec_[0]->FindBin(myEventProxy.event->getNPU());  
  return hPUVec_[0]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getGenWeight(const EventProxyHTT & myEventProxy){

  ///MC weights cab be quite large, but are always +-const.
  ///to avoid counter overflow we keep only sign.
  return myEventProxy.event->getMCWeight()/fabs(myEventProxy.event->getMCWeight());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::vector<std::string> HTTAnalyzer::getTauDecayName(int decModeMinus, int decModePlus){

  std::vector<std::string> types;

  if(decModeMinus==tauDecay1ChargedPion0PiZero && decModePlus==tauDecay1ChargedPion0PiZero) types.push_back("PiPi0Pi0");

  if(isOneProng(decModeMinus) && isOneProng(decModePlus) ) types.push_back("1Prong1Prong");

  if( (decModeMinus==tauDecay1ChargedPion0PiZero && isLepton(decModePlus) ) ||
      (isLepton(decModeMinus) && decModePlus==tauDecay1ChargedPion0PiZero)) types.push_back("Lepton1Prong0Pi0");
    
  if( (isOneProng(decModeMinus) && isLepton(decModePlus) ) ||
      ( isLepton(decModeMinus) && isOneProng(decModePlus) ) ) types.push_back("Lepton1Prong");

  if(decModeMinus==tauDecay1ChargedPion1PiZero && decModePlus==tauDecay1ChargedPion1PiZero ) types.push_back("PiPlusPiMinus2Pi0");


  if( isOneProng(decModeMinus) && decModeMinus!=tauDecay1ChargedPion0PiZero && 
      isOneProng(decModePlus) && decModePlus!=tauDecay1ChargedPion0PiZero )   types.push_back("1Prong1ProngXPi0");

  if(isLepton(decModePlus) && isLepton(decModeMinus)) types.push_back("LeptonLepton");

  return types;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::isOneProng(int decMode){
  if(decMode==tauDecay1ChargedPion0PiZero ||
     decMode==tauDecay1ChargedPion1PiZero ||
     decMode==tauDecay1ChargedPion2PiZero ||
     decMode==tauDecay1ChargedPion3PiZero ) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::isLepton(int decMode){
  if(decMode==tauDecaysElectron || decMode==tauDecayMuon) return true;
  else return false;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getLeptonCorrection(float eta, float pt, hadronicTauDecayModes tauDecayMode){

  if(sampleName.find("Data")!=std::string::npos) return 1.0;
  
  if(tauDecayMode == tauDecayMuon){
    scaleWorkspace->var("m_pt")->setVal(pt);
    scaleWorkspace->var("m_eta")->setVal(eta);
    float muon_id_scalefactor = scaleWorkspace->function("m_id_ratio")->getVal();
    float muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio")->getVal();
    float muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data")->getVal();//OR of the HLT_IsoMu22 and HLT_IsoTkMu22
    return  muon_id_scalefactor*muon_iso_scalefactor*muon_trg_efficiency;
  }
  else if(tauDecayMode == tauDecaysElectron) return 1.0;  
  else{
    if(sampleName.find("H")==std::string::npos &&
       !(sampleName.find("DY")!=std::string::npos && sampleName.find("MatchT")!=std::string::npos) //test: works fine
       ) return 1.0;
    scaleWorkspace->var("t_pt")->setVal(pt);
    scaleWorkspace->var("t_eta")->setVal(eta);
    scaleWorkspace->var("t_dm")->setVal(tauDecayMode);
    float tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_m_pt30_sf")->getVal();
    //float tau_id_scalefactor = 0.9;//according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
    return tau_id_scalefactor;
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
