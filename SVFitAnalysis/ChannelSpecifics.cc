#include "ChannelSpecifics.h"
#include "svfitAnalyzer.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooBinning.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TRandom3.h"
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::ChannelSpecifics(svfitAnalyzer *aAnalyzer){

        myAnalyzer = aAnalyzer;

        initializeLeptonCorrections();
        h2DMuonIdCorrections = 0;
        h3DMuonIsoCorrections = 0;
        h3DMuonTrgCorrections = 0;
        h3DMuonXTrgCorrections = 0;
        h1DMuonTrkCorrections = 0;
        h2DTauTrgGenuineCorrections = 0;
        h2DTauTrgFakeCorrections = 0;
        h2DTauXTrgGenuineCorrections = 0;
        h2DTauXTrgFakeCorrections = 0;
        h3DTauCorrections = 0;

        //initializeLeptonCorrections();
        defineCategories();

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
ChannelSpecifics::~ChannelSpecifics(){

        if(h2DMuonIdCorrections) delete h2DMuonIdCorrections;
        if(h3DMuonIsoCorrections) delete h3DMuonIsoCorrections;
        if(h3DMuonTrgCorrections) delete h3DMuonTrgCorrections;
        if(h3DMuonXTrgCorrections) delete h3DMuonXTrgCorrections;
        if(h1DMuonTrkCorrections) delete h1DMuonTrkCorrections;
        if(h2DTauTrgGenuineCorrections) delete h2DTauTrgGenuineCorrections;
        if(h2DTauTrgFakeCorrections) delete h2DTauTrgFakeCorrections;
        if(h2DTauXTrgGenuineCorrections) delete h2DTauXTrgGenuineCorrections;
        if(h2DTauXTrgFakeCorrections) delete h2DTauXTrgFakeCorrections;
        if(h3DTauCorrections) delete h3DTauCorrections;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::defineCategories(){

        jet0 = new HTTAnalysis::eventCategory("0jet", categoryRejester);
        boosted = new HTTAnalysis::eventCategory("boosted", categoryRejester);
        vbf = new HTTAnalysis::eventCategory("vbf", categoryRejester);

        antiIso_jet0 = new HTTAnalysis::eventCategory("antiIso_0jet", categoryRejester);
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

        inclusive = new HTTAnalysis::eventCategory("inclusive", categoryRejester);
        antiIso_inclusive = new HTTAnalysis::eventCategory("antiIso_inclusive", categoryRejester);

        btag = new HTTAnalysis::eventCategory("btag", categoryRejester);
        nobtag = new HTTAnalysis::eventCategory("nobtag", categoryRejester);
        antiIso_btag = new HTTAnalysis::eventCategory("antiIso_btag", categoryRejester);
        antiIso_nobtag = new HTTAnalysis::eventCategory("antiIso_nobtag", categoryRejester);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ChannelSpecifics::initializeLeptonCorrections(){

  return;//TEST

#pragma omp critical(ROOFIT_INITIALIZATION)
        {

                std::cout<<"Initializing corrections for thread "<<omp_get_thread_num()<<std::endl;
                TFile::SetCacheFileDir("/tmp/");
                std::string correctionFileName = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_v16_4.root";
                TFile *aFile = TFile::Open(correctionFileName.c_str(),"CACHEREAD");

                RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile->Get("w");

                RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
                RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_binned_ratio");
                RooAbsReal *muon_trg_scalefactor = scaleWorkspace->function("m_trgOR4_binned_ratio");//MB 24->22
                RooAbsReal *muon_xtrg_scalefactor_iso  = scaleWorkspace->function("m_trgMu19leg_eta2p1_desy_ratio");//mu-iso<0.15?
                //RooAbsReal *muon_xtrg_scalefactor_aiso = scaleWorkspace->function("m_trgMu19leg_eta2p1_aiso0p15to0p3_desy_ratio");//0.15<mu-iso<0.3?
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
                //Uwaga: ponizsze korekcje sa 2d, do wlozenia w 3d
                TH2F *h2DMuonXTrgCorrections_iso = (TH2F*)muon_xtrg_scalefactor_iso->createHistogram("h2DMuonXTrgCorrections_iso",
									                       *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
									                       RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
									                       RooFit::Extended(kFALSE),
									                       RooFit::Scaling(kFALSE));
                TH2F *h2DMuonXTrgCorrections_aiso = (TH2F*)muon_xtrg_scalefactor_iso->createHistogram("h2DMuonXTrgCorrections_aiso",
									                       *scaleWorkspace->var("m_pt"),RooFit::Binning(1980,10,1000),
									                       RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(48,-2.4,2.4)),//MB m_abs_eta->m_eta
									                       RooFit::Extended(kFALSE),
									                       RooFit::Scaling(kFALSE));
                h3DMuonXTrgCorrections = new TH3F("h3DMuonXTrgCorrections","",1980,10,1000,48,-2.4,2.4,8,-0.05,0.35);//less iso-bins than for id/single-trg due to different range
                for(int iX=0; iX<1980; ++iX){
                  for(int iY=0; iY<48; ++iY){
                    float binContent = h2DMuonXTrgCorrections_iso->GetBinContent(iX+1,iY+1);
                    for(int iZ=0; iZ<4; ++iZ){//bins 1-4 correspond with iso<0.15
                      h3DMuonXTrgCorrections->SetBinContent(iX+1,iY+1,iZ+1,binContent);
		    }
		    binContent = h2DMuonXTrgCorrections_aiso->GetBinContent(iX+1,iY+1);
		    for(int iZ=4; iZ<8; ++iZ){//bins 4-8 correspond with 0.15<iso<0.3
                      h3DMuonXTrgCorrections->SetBinContent(iX+1,iY+1,iZ+1,binContent);
		    }
                  }
                }

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
                //if(h2DMuonXTrgCorrections_iso) delete h2DMuonXTrgCorrections_iso; FIXME
                //if(h2DMuonXTrgCorrections_aiso) delete h2DMuonXTrgCorrections_aiso; FIXME

                ////
                //Open different RooWorkspace for missig SFs - ugly:(

                std::string correctionFileName_b = "http://akalinow.web.cern.ch/akalinow/htt_scalefactors_sm_moriond_v2.root";
                TFile *bFile = TFile::Open(correctionFileName_b.c_str(),"CACHEREAD");

                RooWorkspace *scaleWorkspace_b = (RooWorkspace*)bFile->Get("w");

                RooAbsReal *tau_xtrg_genuine_efficiency = scaleWorkspace_b->function("t_genuine_TightIso_mt_ratio");//MB data->ratio
                RooAbsReal *tau_xtrg_fake_efficiency = scaleWorkspace_b->function("t_fake_TightIso_mt_ratio");//MB data->ratio

                ///WARNING: t_eta and t_dm not used, so histograms have only one bin in this directions
                RooArgSet dependentVarsForMT(*scaleWorkspace_b->var("t_pt"),*scaleWorkspace_b->var("t_eta"));
                RooArgSet projectedVarsForMT;

                const RooAbsReal * tau_xtrg_genuine_efficiency_proj = tau_xtrg_genuine_efficiency->createPlotProjection(dependentVarsForMT,projectedVarsForMT);
                const RooAbsReal * tau_xtrg_fake_efficiency_proj = tau_xtrg_fake_efficiency->createPlotProjection(dependentVarsForMT,projectedVarsForMT);

                h2DTauXTrgGenuineCorrections = (TH2F*)tau_xtrg_genuine_efficiency_proj->createHistogram("h2DTauXTrgGenuineCorrections",
											                *scaleWorkspace_b->var("t_pt"),RooFit::Binning(5000,0,1000),
											                RooFit::YVar(*scaleWorkspace_b->var("t_eta"),RooFit::Binning(4,-2.5,2.5)),//MB binning in eta: barrel/endcaps, i.e. <1.5
											                RooFit::Extended(kFALSE),
											                RooFit::Scaling(kFALSE));

                h2DTauXTrgFakeCorrections = (TH2F*)tau_xtrg_fake_efficiency_proj->createHistogram("h2DTauXTrgFakeCorrections",
										                  *scaleWorkspace_b->var("t_pt"),RooFit::Binning(5000,0,1000),
										                  RooFit::YVar(*scaleWorkspace_b->var("t_eta"),RooFit::Binning(4,-2.5,2.5)),//MB binning in eta: barrel/endcaps, i.e. <1.5
										                  RooFit::Extended(kFALSE),
										                  RooFit::Scaling(kFALSE));

                delete bFile;
        }
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getLeptonCorrection(float eta, float pt, float iso,
                                            HTTAnalysis::hadronicTauDecayModes tauDecayMode,
                                            bool useTauTrigger, int mc_match, bool useXTrigger){

        if(myAnalyzer->sampleName_.find("Data")!=std::string::npos) return 1.0;

        if(!h2DMuonIdCorrections) initializeLeptonCorrections();

        if(tauDecayMode == HTTAnalysis::tauDecayMuon) {
                int iBin = h3DMuonTrgCorrections->FindBin(std::min(pt,(Float_t)999.9), eta, std::min(iso,(Float_t)0.499));
                float trigweight = h3DMuonTrgCorrections->GetBinContent(iBin);
                if(pt<23){//xtrigger mu-tau
                  iBin = h3DMuonXTrgCorrections->FindBin(std::min(pt,(Float_t)999.9), eta, std::min(iso,(Float_t)0.299));
                  trigweight = h3DMuonXTrgCorrections->GetBinContent(iBin);
                }

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
                bool fakeTau = mc_match!=5;

                //according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
                float tau_id_scalefactor = getTauIDSF(eta, mc_match);
                float tau_trg_efficiency = 1.0;

                //if(!fakeTau) tau_id_scalefactor = 0.95;
                if(useTauTrigger) {
                        int tauDecayModeFix = (int)tauDecayMode;
                        if(tauDecayModeFix==2) tauDecayModeFix=1; //DM=2 is not covered by parametrisation.
                        int iBin = h2DTauTrgGenuineCorrections->FindBin(std::min(pt,(Float_t)999.9), tauDecayModeFix);
                        tau_trg_efficiency = h2DTauTrgGenuineCorrections->GetBinContent(iBin);
                        if(fakeTau) {
                                iBin = h2DTauTrgFakeCorrections->FindBin(std::min(pt,(Float_t)999.9), tauDecayModeFix);
                                tau_trg_efficiency = h2DTauTrgFakeCorrections->GetBinContent(iBin);
                        }
                        if(useXTrigger){//xtrigger mu-tau
                          if(!fakeTau){ //genuine tau
	                          iBin = h2DTauXTrgGenuineCorrections->FindBin(std::min(pt,(Float_t)999.9),eta);
	                          tau_trg_efficiency = h2DTauXTrgGenuineCorrections->GetBinContent(iBin);
                          }
                          else{ //fake tau
	                          iBin = h2DTauXTrgFakeCorrections->FindBin(std::min(pt,(Float_t)999.9), eta);
	                          tau_trg_efficiency = h2DTauXTrgFakeCorrections->GetBinContent(iBin);
                          }
                        }
                }
                return tau_id_scalefactor*tau_trg_efficiency;
        }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getDYReweight(const std::string & categoryName, const HTTAnalysis::sysEffects & aSystEffect){

        float weight = 1.0;
///Valuie read from Fig. 14 of AN2016_355_v9
        if(categoryName.find("0jet")!=std::string::npos) { weight = 1.02;}
        else if(categoryName.find("boosted")!=std::string::npos) {

                float higgsPt =  (myAnalyzer->leg1_.getP4(aSystEffect) + myAnalyzer->leg2_.getP4(aSystEffect) + myAnalyzer->MET_.getP4(aSystEffect)).Pt();

                if(higgsPt<100) weight = 1.02;
                else if (higgsPt>100 && higgsPt<200) weight = 1.0;
                else if (higgsPt>200 && higgsPt<250) weight = 1.02;
                else if (higgsPt>250) weight = 0.98;
        }
        else if(categoryName.find("vbf")!=std::string::npos) {

                float jetsMass = (myAnalyzer->jet1_.getP4(aSystEffect)+myAnalyzer->jet2_.getP4(aSystEffect)).M();

                if(jetsMass<700) weight = 1.15;
                else if (jetsMass>700 && jetsMass<1100) weight = 1.0;
                else if (jetsMass>1100 && jetsMass<1500) weight = 0.92;
                else if (jetsMass>1500 && jetsMass<2000) weight = 0.9;
        }
        return weight;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getTauIDSF(float eta, int mc_match){

        if(mc_match == 5) return 0.95;
        if(mc_match == 2 || mc_match == 4){
          int iBin = tauID_FRSF_mu->FindBin(std::abs(eta));
          return tauID_FRSF_mu->GetBinContent(iBin);
          }
        if(mc_match == 1 || mc_match == 3){
          int iBin = tauID_FRSF_ele->FindBin(std::abs(eta));
          return tauID_FRSF_ele->GetBinContent(iBin);
          }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
