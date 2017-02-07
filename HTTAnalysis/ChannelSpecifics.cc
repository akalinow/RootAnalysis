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

#pragma omp critical
        {
                std::string filePath = "htt_scalefactors_v5.root";
                TFile aFile(filePath.c_str());
                RooWorkspace *scaleWorkspace = (RooWorkspace*)aFile.Get("w");

                RooAbsReal *muon_id_scalefactor = scaleWorkspace->function("m_id_ratio");
                RooAbsReal *muon_iso_scalefactor = scaleWorkspace->function("m_iso_ratio");
                RooAbsReal *muon_trg_efficiency = scaleWorkspace->function("m_trgOR_data");//OR of the HLT_IsoMu22 and HLT_IsoTkMu22
                //RooAbsReal *tau_id_scalefactor = scaleWorkspace->function("t_iso_mva_m_pt30_sf");

                h2DMuonIdCorrections = (TH2F*)muon_id_scalefactor->createHistogram("h2DMuonIdCorrections",
                                                                                   *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                   RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                   RooFit::Scaling(kFALSE));

                h2DMuonIsoCorrections = (TH2F*)muon_iso_scalefactor->createHistogram("h2DMuonIsoCorrections",
                                                                                     *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                     RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                     RooFit::Scaling(kFALSE));

                h2DMuonTrgCorrections = (TH2F*)muon_trg_efficiency->createHistogram("h2DMuonTrgCorrections",
                                                                                    *scaleWorkspace->var("m_pt"),RooFit::Binning(300,0,300),
                                                                                    RooFit::YVar(*scaleWorkspace->var("m_eta"),RooFit::Binning(10,-2.1,2.1)),
                                                                                    RooFit::Scaling(kFALSE));
        }

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
float ChannelSpecifics::getLeptonCorrection(float eta, float pt, HTTAnalysis::hadronicTauDecayModes tauDecayMode){

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
                if(myAnalyzer->sampleName.find("H")==std::string::npos &&
                   !(myAnalyzer->sampleName.find("A")!=std::string::npos &&
                     myAnalyzer->sampleName.find("All")==std::string::npos) && //pseudoscalar
                   myAnalyzer->sampleName.find("MatchT")==std::string::npos
                   ) return 1.0;
                float tau_id_scalefactor = 0.9;//according to https://twiki.cern.ch/twiki/bin/view/CMS/SMTauTau2016#MC_corrections
                return tau_id_scalefactor;
        }
        return 1.0;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
