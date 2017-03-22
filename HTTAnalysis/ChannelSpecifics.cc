#include "ChannelSpecifics.h"
#include "HTTAnalyzer.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooBinning.h"
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
        h3DMuonIsoCorrections = 0;
        h3DMuonTrgCorrections = 0;
        h1DMuonTrkCorrections = 0;
        h2DTauTrgGenuineCorrections = 0;
        h2DTauTrgFakeCorrections = 0;
        h3DTauCorrections = 0;

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::~ChannelSpecifics(){

        if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
        if(h3DMuonIsoCorrections) delete h3DMuonIsoCorrections;
        if(h3DMuonTrgCorrections) delete h3DMuonTrgCorrections;
        if(h1DMuonTrkCorrections) delete h1DMuonTrkCorrections;
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
        rho_rho = new HTTAnalysis::eventCategory("rho_rho", categoryRejester);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializeCorrections(){

#pragma omp critical(ROOFIT_INITIALIZATION)
        {

                std::cout<<"Initializing corrections for thread "<<omp_get_thread_num()<<std::endl;

                std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v16_4.root";
                TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");

                RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");

                RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
                RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_binned_ratio");
                RooAbsReal *muon_trg_scalefactor = scaleWorkspace->function("m_trgOR4_binned_ratio");//MB 24->22
                RooAbsReal *muon_trk_scalefactor = scaleWorkspace->function("m_trk_ratio");//MB not in HTTAnalysis
                RooAbsReal *tau_trg_genuine_efficiency = scaleWorkspace->function("t_genuine_TightIso_tt_ratio");//MB data->ratio
                RooAbsReal *tau_trg_fake_efficiency = scaleWorkspace->function("t_fake_TightIso_tt_ratio");//MB data->ratio

                h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
                                                                                   *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
                                                                                   RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
                                                                                   RooFit::Extended(kFALSE),
                                                                                   RooFit::Scaling(kFALSE));
                //
                RooArgSet dependentVarsForMu(*scaleWorkspace->var("m_pt"),*scaleWorkspace->var("m_eta"),*scaleWorkspace->var("m_iso"));
                RooArgSet projectedVarsForMu;
                const RooAbsReal *muon_iso_scalefactor_proj = muon_iso_scalefactor->createPlotProjection(dependentVarsForMu,projectedVarsForMu);
                h3DMuonIsoCorrections = (TH3F*)muon_iso_scalefactor_proj->createHistogram("h3DMuonIsoCorrections",
                                                                                          *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
                                                                                          RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),
                                                                                          RooFit::ZVar(*scaleWorkspace->var("m_iso"),RooFit::Binning(12,-0.05,0.55)),
                                                                                          RooFit::Extended(kFALSE),
                                                                                          RooFit::Scaling(kFALSE));
                const RooAbsReal *muon_trg_scalefactor_proj = muon_trg_scalefactor->createPlotProjection(dependentVarsForMu,projectedVarsForMu);
                h3DMuonTrgCorrections = (TH3F*)muon_trg_scalefactor_proj->createHistogram("h3DMuonTrgCorrections",
                                                                                          *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
                                                                                          RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),
                                                                                          RooFit::ZVar(*scaleWorkspace->var("m_iso"),RooFit::Binning(12,-0.05,0.55)),
                                                                                          RooFit::Extended(kFALSE),
                                                                                          RooFit::Scaling(kFALSE));
                h1DMuonTrkCorrections = (TH1F*)muon_trk_scalefactor->createHistogram("h1DMuonTrkCorrections",
                                                                                     *scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4),//MB m_abs_eta->m_eta
                                                                                     RooFit::Extended(kFALSE),
                                                                                     RooFit::Scaling(kFALSE));

                ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
                RooArgSet dependentVars(*scaleWorkspace->var("t_pt"),*scaleWorkspace->var("t_dm"));
                RooArgSet projectedVars;

                const RooAbsReal * tau_trg_genuine_efficiency_proj = tau_trg_genuine_efficiency->createPlotProjection(dependentVars,projectedVars);
                const RooAbsReal * tau_trg_fake_efficiency_proj = tau_trg_fake_efficiency->createPlotProjection(dependentVars,projectedVars);

		/*
                RooBinning binsForTauTrg(0,1000);
                binsForTauTrg.addUniform(1, 0, 30);
                binsForTauTrg.addUniform(2000, 30, 130);
                binsForTauTrg.addUniform(1000, 130, 330);
                binsForTauTrg.addUniform(1340, 330, 1000);
		*/
                h2DTauTrgGenuineCorrections = (TH2F*)tau_trg_genuine_efficiency_proj->createHistogram("h2DTauTrgGenuineCorrections",
                                                                                                      *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                                                      //*scaleWorkspace->var("t_pt"),RooFit::Binning(binsForTauTrg),
                                                                                                      RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(11,-0.5,10.5)),//MB proper binning (3->11) for DMs==0,1(2),10
                                                                                                      RooFit::Extended(kFALSE),
                                                                                                      RooFit::Scaling(kFALSE));

                h2DTauTrgFakeCorrections = (TH2F*)tau_trg_fake_efficiency_proj->createHistogram("h2DTauTrgFakeCorrections",
                                                                                                *scaleWorkspace->var("t_pt"),RooFit::Binning(5000,0,1000),
                                                                                                //*scaleWorkspace->var("t_pt"),RooFit::Binning(binsForTauTrg),
                                                                                                RooFit::YVar(*scaleWorkspace->var("t_dm"),RooFit::Binning(11,-0.5,10.5)),//MB proper binning (3->11) for DMs==0,1(2),10
                                                                                                RooFit::Extended(kFALSE),
                                                                                                RooFit::Scaling(kFALSE));

                delete aFile;
        }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getLeptonCorrection(float eta, float pt, float iso,
                                            HTTAnalysis::hadronicTauDecayModes tauDecayMode,
                                            bool useTauTrigger){

        if(myAnalyzer->sampleName.find("Data")!=std::string::npos) return 1.0;

        if(!h2DMuonIdCorrections) initializeCorrections();

        if(tauDecayMode == HTTAnalysis::tauDecayMuon) {
                int iBin = h3DMuonTrgCorrections->FindBin(std::min(pt,(Float_t)999.9), eta, std::min(iso,(Float_t)0.499));
                float trigweight = h3DMuonTrgCorrections->GetBinContent(iBin);

                iBin = h2DMuonIdCorrections->FindBin(std::min(pt,(Float_t)999.9), eta);
                float idweight = h2DMuonIdCorrections->GetBinContent(iBin);

                iBin = h3DMuonIsoCorrections->FindBin(std::min(pt,(Float_t)999.9), eta, std::min(iso,(Float_t)0.499));
                float isoweight = h3DMuonIsoCorrections->GetBinContent(iBin);

                iBin = h1DMuonTrkCorrections->FindBin(eta);
                float trackingweight = h1DMuonTrkCorrections->GetBinContent(iBin);

                return trigweight*idweight*isoweight*trackingweight;
        }
        else if(tauDecayMode == HTTAnalysis::tauDecaysElectron) return 1.0;
        else{
                bool fakeTau = myAnalyzer->sampleName.find("HTT")==std::string::npos &&
                               myAnalyzer->sampleName.find("ATT")==std::string::npos &&
                               myAnalyzer->sampleName.find("MatchT")==std::string::npos;

                //according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
                float tau_id_scalefactor = 1.0;
                float tau_trg_efficiency = 1.0;

                if(!fakeTau) tau_id_scalefactor = 0.95;
                if(useTauTrigger) {
                        int tauDecayModeFix = (int)tauDecayMode;
                        if(tauDecayModeFix==2) tauDecayModeFix=1;//DM=2 is not covered by parametrisation.
                        int iBin = h2DTauTrgGenuineCorrections->FindBin(std::min(pt,(Float_t)999.9), tauDecayModeFix);
                        tau_trg_efficiency = h2DTauTrgGenuineCorrections->GetBinContent(iBin);
                        if(fakeTau) {
                                iBin = h2DTauTrgFakeCorrections->FindBin(std::min(pt,(Float_t)999.9), tauDecayModeFix);
                                tau_trg_efficiency = h2DTauTrgFakeCorrections->GetBinContent(iBin);
                        }
                }
                return tau_id_scalefactor*tau_trg_efficiency;
        }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
