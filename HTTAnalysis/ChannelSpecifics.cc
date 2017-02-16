#include "ChannelSpecifics.h"
#include "HTTAnalyzer.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::ChannelSpecifics(HTTAnalyzer *aAnalyzer){

        initializeCorrections();

        myAnalyzer = aAnalyzer;
        h2DMuonIdCorrections = 0;
        h2DMuonIsoCorrections = 0;
        h2DMuonTrgCorrections = 0;
        h3DTauCorrections = 0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::~ChannelSpecifics(){

        if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
        if(h2DMuonIsoCorrections) delete h2DMuonIsoCorrections;
        if(h2DMuonTrgCorrections) delete h2DMuonTrgCorrections;
        if(h3DTauCorrections) delete h3DTauCorrections;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializeCorrections(){

#pragma omp critical(ROOFIT_INITIALIZATION)
        {
        std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v5.root";
        TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");
        RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");

        RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
        RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio");
        RooAbsReal *muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data");         //OR of the HLT_IsoMu22 and HLT_IsoTkMu22
        //RooAbsReal *tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_m_pt30_sf");
        RooAbsReal *tau_trgOS_efficiency = scaleWorkspace->function("t_trgTightIso_data");
        RooAbsReal *tau_trgSS_efficiency = scaleWorkspace->function("t_trgTightIsoSS_data");

        h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
                                                                           *scaleWorkspace->var("m_pt"),RooFit::Binning(400,0,200),
                                                                           RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                           RooFit::Scaling(kFALSE));

        h2DMuonIsoCorrections = (TH2F*)muon_iso_scalefactor->createHistogram("h2DMuonIsoCorrections",
                                                                             *scaleWorkspace->var("m_pt"),RooFit::Binning(400,0,200),
                                                                             RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                             RooFit::Scaling(kFALSE));

        h2DMuonTrgCorrections = (TH2F*)muon_trg_efficiency->createHistogram("h2DMuonTrgCorrections",
                                                                            *scaleWorkspace->var("m_pt"),RooFit::Binning(400,0,200),
                                                                            RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                            RooFit::Scaling(kFALSE));
        ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
        h3DTauTrgOSCorrections = (TH3F*)tau_trgOS_efficiency->createHistogram("h3DTauTrgOSCorrections",
                                                                              *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                              RooFit::YVar(*scaleWorkspace->var("t_eta"),RooFit::Binning(1,-2.4,2.4)),
                                                                              RooFit::ZVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(1,-0.5,10.5)),
                                                                              RooFit::Scaling(kFALSE));

        h3DTauTrgSSCorrections = (TH3F*)tau_trgSS_efficiency->createHistogram("h3DTauTrgSSCorrections",
                                                                              *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                              RooFit::YVar(*scaleWorkspace->var("t_eta"),RooFit::Binning(1,-2.4,2.4)),
                                                                              RooFit::ZVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(1,-0.5,10.5)),
                                                                              RooFit::Scaling(kFALSE));

        delete aFile;
                                                                            }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getLeptonCorrection(float eta, float pt,
                                            HTTAnalysis::hadronicTauDecayModes tauDecayMode,
                                            bool useTauTrigger){

        if(myAnalyzer->sampleName.find("Data")!=std::string::npos) return 1.0;

        if(!h2DMuonIdCorrections) initializeCorrections();

        if(tauDecayMode == HTTAnalysis::tauDecayMuon) {
                int iBin = h2DMuonIdCorrections->FindBin(pt, eta);
                float muon_id_scalefactor = h2DMuonIdCorrections->GetBinContent(iBin);
                float muon_iso_scalefactor = h2DMuonIsoCorrections->GetBinContent(iBin);
                float muon_trg_efficiency = h2DMuonTrgCorrections->GetBinContent(iBin);
                return muon_id_scalefactor*muon_iso_scalefactor*muon_trg_efficiency;
        }
        else if(tauDecayMode == HTTAnalysis::tauDecaysElectron) return 1.0;
        else{
                //according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
                float tau_id_scalefactor = 0.95; //h3DTauIdCorrections->GetBinContent(iBin);
                float tau_trg_efficiency = 1.0;
                if(useTauTrigger) {
                        int iBin = h3DTauTrgOSCorrections->FindBin(pt, eta, (int)tauDecayMode);
                        tau_trg_efficiency = h3DTauTrgOSCorrections->GetBinContent(iBin);
                        if(myAnalyzer->sampleName.find("HTT")==std::string::npos &&
                           myAnalyzer->sampleName.find("ATT")==std::string::npos &&
                           myAnalyzer->sampleName.find("MatchT")==std::string::npos) {
                                tau_id_scalefactor = 1.0;
                                int iBin = h3DTauTrgSSCorrections->FindBin(pt, eta, (int)tauDecayMode);
                                tau_trg_efficiency = h3DTauTrgSSCorrections->GetBinContent(iBin);
                        }
                }
                return tau_id_scalefactor*tau_trg_efficiency;
        }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
