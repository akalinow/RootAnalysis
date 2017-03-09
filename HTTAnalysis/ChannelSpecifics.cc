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
        defineCategories();

        myAnalyzer = aAnalyzer;
        h2DMuonIdCorrections = 0;
        h2DMuonIsoCorrections = 0;
        h2DMuonTrgCorrections = 0;
        h2DTauTrgGenuineCorrections = 0;
        h2DTauTrgFakeCorrections = 0;
        h3DTauCorrections = 0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::~ChannelSpecifics(){

        if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
        if(h2DMuonIsoCorrections) delete h2DMuonIsoCorrections;
        if(h2DMuonTrgCorrections) delete h2DMuonTrgCorrections;

        if(h2DTauTrgGenuineCorrections) delete h2DTauTrgGenuineCorrections;
        if(h2DTauTrgFakeCorrections) delete h2DTauTrgFakeCorrections;
        if(h3DTauCorrections) delete h3DTauCorrections;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::defineCategories(){

        jet0 = new HTTAnalysis::eventCategory("jet0", categoryRejester);
        boosted = new HTTAnalysis::eventCategory("boosted", categoryRejester);
        vbf = new HTTAnalysis::eventCategory("vbf", categoryRejester);

        antiIso_jet0 = new HTTAnalysis::eventCategory("antiIso_jet0", categoryRejester);
        antiIso_boosted = new HTTAnalysis::eventCategory("antiIso_boosted", categoryRejester);
        antiIso_vbf = new HTTAnalysis::eventCategory("antiIso_vbf", categoryRejester);

//zmumuShape_jet0 = new HTTAnalysis::eventCategory("zmumuShape_jet0", categoryRejester);
//zmumuShape_boosted = new HTTAnalysis::eventCategory("zmumuShape_boosted", categoryRejester);
//zmumuShape_vbf = new HTTAnalysis::eventCategory("zmumuShape_vbf", categoryRejester);

        mu_pi = new HTTAnalysis::eventCategory("mu_pi", categoryRejester);
        mu_rho = new HTTAnalysis::eventCategory("mu_rho", categoryRejester);
        pi_pi = new HTTAnalysis::eventCategory("pi_pi", categoryRejester);
        pi_rho = new HTTAnalysis::eventCategory("pi_pi", categoryRejester);
        rho_rho = new HTTAnalysis::eventCategory("pi_pi", categoryRejester);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializeCorrections(){

#pragma omp critical(ROOFIT_INITIALIZATION)
        {
                std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v16_3.root";
                TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");
                RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");

                RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
                RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio");
                RooAbsReal *muon_trg_scalefactor = scaleWorkspace->function("m_trgIsoMu24orTkIsoMu24_desy_ratio");              
                RooAbsReal *tau_trg_genuine_efficiency = scaleWorkspace->function("t_genuine_TightIso_tt_data");
                RooAbsReal *tau_trg_fake_efficiency = scaleWorkspace->function("t_fake_TightIso_tt_data");

                h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
                                                                                   *scaleWorkspace->var("m_pt"),RooFit::Binning(400,10,200),
                                                                                   RooFit::YVar(*scaleWorkspace->var("m_abs_eta"),RooFit::Binning(10,0,2.4)),
                                                                                   RooFit::Extended(kFALSE),
                                                                                   RooFit::Scaling(kFALSE));

                h2DMuonIsoCorrections = (TH2F*)muon_iso_scalefactor->createHistogram("h2DMuonIsoCorrections",
                                                                                     *scaleWorkspace->var("m_pt"),RooFit::Binning(400,10,200),
                                                                                     RooFit::YVar(*scaleWorkspace->var("m_abs_eta"),RooFit::Binning(10,0,2.4)),
                                                                                     RooFit::Extended(kFALSE),
                                                                                     RooFit::Scaling(kFALSE));

                h2DMuonTrgCorrections = (TH2F*)muon_trg_scalefactor->createHistogram("h2DMuonTrgCorrections",
                                                                                     *scaleWorkspace->var("m_pt"),RooFit::Binning(400,10,200),
                                                                                     RooFit::YVar(*scaleWorkspace->var("m_abs_eta"),RooFit::Binning(10,0,2.4)),
                                                                                     RooFit::Extended(kFALSE),
                                                                                     RooFit::Scaling(kFALSE));
                ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
                RooArgSet dependentVars(*scaleWorkspace->var("t_pt"),*scaleWorkspace->var("t_dm"));
                RooArgSet projectedVars;

                const RooAbsReal * tau_trg_genuine_efficiency_proj = tau_trg_genuine_efficiency->createPlotProjection(dependentVars,projectedVars);
                const RooAbsReal * tau_trg_fake_efficiency_proj = tau_trg_fake_efficiency->createPlotProjection(dependentVars,projectedVars);

                h2DTauTrgGenuineCorrections = (TH2F*)tau_trg_genuine_efficiency_proj->createHistogram("h3DTauTrgGenuineCorrections",
                                                                                                      *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                                                      RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(3,-0.5,10.5)),
                                                                                                      RooFit::Extended(kFALSE),
                                                                                                      RooFit::Scaling(kFALSE));

                h2DTauTrgFakeCorrections = (TH2F*)tau_trg_fake_efficiency_proj->createHistogram("h3DTauTrgFakeCorrections",
                                                                                                *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                                                RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(3,-0.5,10.5)),
                                                                                                RooFit::Extended(kFALSE),
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
                float muon_trg_scalefactor = h2DMuonTrgCorrections->GetBinContent(iBin);
                return muon_id_scalefactor*muon_iso_scalefactor*muon_trg_scalefactor;
        }
        else if(tauDecayMode == HTTAnalysis::tauDecaysElectron) return 1.0;
        else{
                //according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
                float tau_id_scalefactor = 0.95; //h3DTauIdCorrections->GetBinContent(iBin);
                float tau_trg_efficiency = 1.0;
                if(useTauTrigger) {
                        int iBin = h2DTauTrgGenuineCorrections->FindBin(pt, (int)tauDecayMode);
                        tau_trg_efficiency = h2DTauTrgGenuineCorrections->GetBinContent(iBin);
                        if(myAnalyzer->sampleName.find("HTT")==std::string::npos &&
                           myAnalyzer->sampleName.find("ATT")==std::string::npos &&
                           myAnalyzer->sampleName.find("MatchT")==std::string::npos) {
                                tau_id_scalefactor = 1.0;
                                int iBin = h2DTauTrgFakeCorrections->FindBin(pt, (int)tauDecayMode);
                                tau_trg_efficiency = h2DTauTrgFakeCorrections->GetBinContent(iBin);
                        }
                }
                return tau_id_scalefactor*tau_trg_efficiency;
        }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
